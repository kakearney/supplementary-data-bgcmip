% Comparing ice edge location across hindcasts
%
% Comparing: H16, K20, new to SSMI
%
% Focusing on 1990-2019 b/c that's the overlap between all 4 datasets

% yr = 1990:2019;
% nyr = length(yr);

grdfile = fullfile(moxdir, 'ROMS_Datasets', 'grids', 'AlaskaGrids_Bering10K.nc');
Grd = ncstruct(grdfile);
[nxi, neta] = size(Grd.h);

%% SSMI historical and near-real-time data
%  Historical data should be redownloaded for recent years as it becomes
%  available.  We prefer that when available, but fall back on the
%  near-real-time dataset otherwise.

latlim = minmax(Grd.lat_rho, 'expand');
lonlim = minmax(Grd.lon_rho, 'expand');

% For coordinates, use pre-calculated values from arcticseaice

fid = fopen('psn25lons_v3.dat');
Ssmi.lon = fread(fid,[304 448],'long','ieee-le')/100000;
fclose(fid);
Ssmi.lon = wrapTo360(Ssmi.lon);

fid = fopen('psn25lats_v3.dat');
Ssmi.lat = fread(fid,[304 448],'long','ieee-le')/100000;
fclose(fid);

mask = Ssmi.lon >= lonlim(1) & Ssmi.lon <= lonlim(2) & Ssmi.lat >= latlim(1) & Ssmi.lat <= latlim(2);

fint = scatteredInterpolant(Ssmi.lon(mask), Ssmi.lat(mask), zeros(nnz(mask),1), 'natural', 'none');

% Interpolate historical SSMR data (daily) to ROMS grid

ssmifol =  '/Volumes/LaCie2023/SeaIce_SSMI';
fol1 = 'b10k_regrid';

if ~exist(fullfile(ssmifol, fol1), 'dir')
    mkdir(fullfile(ssmifol, fol1));
end

F = dir(fullfile(ssmifol, 'NSIDC0051_SEAICE_PS_N25km_*_v2.0.nc'));

tstr = strrep(strrep({F.name}, 'NSIDC0051_SEAICE_PS_N25km_', ''), '_v2.0.nc', '');
isdaily = cellfun(@length, tstr) == 8;
tstr = tstr(isdaily);
F = F(isdaily);

t = datetime(tstr, 'InputFormat', 'yyyyMMdd');

yrs = unique(year(t));
for iy = 1:length(yrs)
    disp(yrs(iy));
    newfile = fullfile(ssmifol, fol1, sprintf('ssmr_%d.mat', yrs(iy)));
    if ~exist(newfile, 'file')
    
        isin = year(t) == yrs(iy);

        nfile = nnz(isin);
        DataTmp = struct('t', t(isin), ...
                      'aice', nan(nxi, neta, nfile));

        ftmp = fullfile({F(isin).folder}, {F(isin).name});

        for ii = 1:length(ftmp)
            fprintf('  %d/%d\n', ii, nfile);
            
            I = ncinfo(ftmp{ii});
            isicecon = endsWith({I.Variables.Name}, 'ICECON');
            if any(isicecon)
                fld = I.Variables(isicecon).Name;
                Tmp = ncstruct(ftmp{ii}, fld);            
                Tmp.(fld)(Tmp.(fld) > 1) = NaN;
                fint.Values = Tmp.(fld)(mask);
                DataTmp.aice(:,:,ii) = fint(Grd.lon_rho, Grd.lat_rho);
            end
        end
        nodata = all(isnan(cube2rect(DataTmp.aice)), 2);
        DataTmp.t = DataTmp.t(~nodata);
        DataTmp.aice = DataTmp.aice(:,:,~nodata);
        
        save(newfile, '-struct', 'DataTmp');
    end    
end

% Same for near-real-time DMSP data (using F18 satellite data for now... no
% real reason, just seems to be the most recent one)

F = dir(fullfile(ssmifol, 'NSIDC0081_SEAICE_PS_N25km_*_v2.0.nc'));

tstr = strrep(strrep({F.name}, 'NSIDC0081_SEAICE_PS_N25km_', ''), '_v2.0.nc', '');
isdaily = cellfun(@length, tstr) == 8;
tstr = tstr(isdaily);
F = F(isdaily);

t = datetime(tstr, 'InputFormat', 'yyyyMMdd');

yrs = unique(year(t));

for iy = 1:length(yrs)
    disp(yrs(iy));
    newfile = fullfile(ssmifol, fol1, sprintf('dmsp_%d.mat', yrs(iy)));
    if ~exist(newfile, 'file')
    
        isin = year(t) == yrs(iy);

        nfile = nnz(isin);
        DataTmp = struct('t', t(isin), ...
                      'aice', nan(nxi, neta, nfile));

        ftmp = fullfile({F(isin).folder}, {F(isin).name});

        for ii = 1:length(ftmp)
            fprintf('  %d/%d\n', ii, nfile);
            Tmp = ncstruct(ftmp{ii});
            if isfield(Tmp, 'F18_ICECON')
                Tmp.F18_ICECON(Tmp.F18_ICECON > 1) = NaN;
                fint.Values = Tmp.F18_ICECON(mask);
                DataTmp.aice(:,:,ii) = fint(Grd.lon_rho, Grd.lat_rho);
            end
        end
        nodata = all(isnan(cube2rect(DataTmp.aice)), 2);
        DataTmp.t = DataTmp.t(~nodata);
        DataTmp.aice = DataTmp.aice(:,:,~nodata);
        
        save(newfile, '-struct', 'DataTmp');
    end    
end

%% Setup: Measure ice along 170W and total on shelf

linelon = -168.5;
[Line.lat, Line.lon] = interpm([54 68], ones(1,2)*linelon, km2deg(10));
Line.lon = wrapTo360(Line.lon);

% ROMS grid cells along line

f = scatteredInterpolant(Grd.lon_rho(:), Grd.lat_rho(:), (1:numel(Grd.h))', 'nearest', 'none');

Line.idx = f(Line.lon, Line.lat);
isn = isnan(Line.idx);
Line = structfun(@(x) x(~isn), Line, 'uni', 0);

% % SSMI grid cells along line
% 
% ssmifol = '/Volumes/KYA/SSMIice';
% Ssmi = load(fullfile(ssmifol, 'ssmiBeringCoords.mat'));
% 
% fs = scatteredInterpolant(Ssmi.lon(:), Ssmi.lat(:), (1:numel(Ssmi.lat))', 'nearest', 'none');
% 
% Ssmi.idx = unique(fs(Line.lon, Line.lat), 'stable');

% Shelf region mask (includes Gulf of Anadyr, NBS)

M = load(fullfile(moxdir, 'bering10k', 'output', 'forecasts', 'cmip6', 'analysis_kak', 'nbudget_masks.mat'));
mask = nan(size(Grd.h));
for ii = 1:size(M.masks,1)
    mask(M.masks{ii,2}) = ii;
    [ln, lt] = mask2poly(Grd.lon_psi, Grd.lat_psi, M.masks{ii,2}(2:end-1,2:end-1));
    pmask{ii} = polyshape(ln, lt);
end
pshelf = union([pmask{1:8}]);

w = Grd.area_feast;

grdmask = reshape(pshelf.isinterior(Grd.lon_rho(:), Grd.lat_rho(:)), size(Grd.h));
% obsmask = reshape(pshelf.isinterior(Ssmi.lon(:), Ssmi.lat(:)), size(Ssmi.lon));

% f = scatteredInterpolant(Ssmi.lon, Ssmi.lat(:), zeros(numel(Ssmi.lat),1), 'nearest', 'none');
% a = alphaShape(wrapTo360(Ssmi.lon(:)), Ssmi.lat(:));
% amask = inShape(a, Grd.lon_rho, Grd.lat_rho);

% maxyearlyice = nan([size(Grd.h), nyr, 4]);



%% ... SSMI coordinate data



% [~,Ssmi.lat,Ssmi.lon] = arcticseaice('noplot', datenum(2000,1,1)); % date doesn't matter
% [nrow, ncol] = size(Ssmi.lat);
% [nxi, neta] = size(Grd.h);
% 
% latlim = minmax(Grd.lat_rho, 'expand');
% lonlim = minmax(Grd.lon_rho, 'expand');
% isnear = Ssmi.lat >= latlim(1) & ...
%          Ssmi.lat <= latlim(2) & ...
%          wrapTo360(Ssmi.lon) >= lonlim(1) & ...
%          wrapTo360(Ssmi.lon) <= lonlim(2);
% [r,c] = ind2sub(size(Ssmi.lat), find(isnear));     
% rlim = minmax(r);
% clim = minmax(c);
% 
% lat = Ssmi.lat(rlim(1):rlim(2), clim(1):clim(2));
% lon = wrapTo360(Ssmi.lon(rlim(1):rlim(2), clim(1):clim(2)));
% 
% save ssmiBeringCoords lat lon;


%% ... Observations: SSMI

thresh = 0.1;
ncell = 5;

yrs = 1978:year(datetime('today'));
nyr = length(yrs);

Obs = struct;
Obs.t = cell(nyr,1);
Obs.aice = cell(nyr,1);
Obs.gt10 = cell(nyr,1);
Obs.ltedge = cell(nyr,1);

for iy = 1:length(yrs)
    ftmp1 = fullfile(ssmifol, fol1, sprintf('ssmr_%d.mat', yrs(iy)));
    ftmp2 = fullfile(ssmifol, fol1, sprintf('dmsp_%d.mat', yrs(iy)));
    
    if exist(ftmp1, 'file')
        A = load(ftmp1);
    elseif exist(ftmp2, 'file')
        A = load(ftmp2);
    else
        continue;
    end
    
    Obs.t{iy} = A.t';
    Obs.aice{iy} = local(A.aice,         grdmask, 'weight', Grd.area_feast, 'omitnan');
    Obs.gt10{iy} = local(A.aice>=thresh, grdmask, 'weight', Grd.area_feast, 'omitnan');
    Obs.ltedge{iy} = nan(size(Obs.aice{iy}));

    aiceline = reshape(A.aice, [], size(A.aice,3));
    aiceline = aiceline(Line.idx,:);
    for it = 1:size(aiceline,2)
        if any(aiceline(:,it) > thresh)
            [b,n,bi] = RunLength(aiceline(:,it)>= thresh);
            isice = b == 1 & n' > ncell; 
            if any(isice)
                idxedge = find(isice,1,'last');
                Obs.ltedge{iy}(it) = Grd.lat_rho(Line.idx(bi(idxedge)));
            end
        end
    end
end

% isemp = cellfun(@isempty, Obs.t);

Obs = structfun(@(x) cat(1, x{:}), Obs, 'uni', 0);
Data.Obs = table2timetable(struct2table(Obs));


% % Which coordinates are closest to our line
% 
% % ell = referenceEllipsoid('earth');
% % dfun = @(a,b) distance(a(:,1), a(:,2), b(:,1), b(:,2), ell);
% % [d,idx] = pdist2([lat(:) lon(:)], [Line.lat(:) Line.lon(:)], dfun, 'smallest', 1);
% 
% % Read data from previously-built files
% 
% Obs = struct;
% Obs.t = cell(nyr,1);
% Obs.aice = cell(nyr,1);
% Obs.gt10 = cell(nyr,1);
% Obs.ltedge = cell(nyr,1);
% 
% ssmifol = '/Volumes/KYA/SSMIice';
% 
% for ii = 1:nyr
%     A = load(fullfile(ssmifol, sprintf('ssmiBeringOct%d-Oct%d.mat', yr(ii), yr(ii)+1)));
%     
%     Obs.t{ii} = datetime(A.t', 'ConvertFrom', 'datenum');
%     Obs.aice{ii} = local(A.ssmifrac./100,         obsmask, 'omitnan');
%     Obs.gt10{ii} = local(A.ssmifrac/.100>=thresh, obsmask, 'omitnan');
%     
% %     f.Values = reshape(max(A.ssmifrac, [], 3), [], 1);
% %     tmp = f(Grd.lon_rho, Grd.lat_rho)./100;
% %     tmp(~amask) = NaN;
% %     maxyearlyice(:,:,ii,1) = tmp;
%     Obs.ltedge{ii} = nan(size(Obs.aice{ii}));
% 
%     aiceline = reshape(A.ssmifrac/100, [], size(A.ssmifrac,3));
%     aiceline = aiceline(Ssmi.idx,:);
%     for it = 1:size(aiceline,2)
%         if any(aiceline(:,it) > thresh)
%             [b,n,bi] = RunLength(aiceline(:,it)>= thresh);
%             isice = b == 1 & n' > ncell; 
%             if any(isice)
%                 idxedge = find(isice,1,'last');
%                 Obs.ltedge{ii}(it) = Ssmi.lat(Ssmi.idx(bi(idxedge)));
%             end
%         end
%     end
% end
% 
% Obs = structfun(@(x) cat(1, x{:}), Obs, 'uni', 0);
% 
% Data.Obs = table2timetable(struct2table(Obs));

%% Observations take 2: SSMI and DMSP separately

thresh = 0.1;
ncell = 5;

yrs = 1978:year(datetime('today'));
nyr = length(yrs);

Ssmi = struct;
Ssmi.t = cell(nyr,1);
Ssmi.aice = cell(nyr,1);
Ssmi.gt10 = cell(nyr,1);

Dmsp = Ssmi;

for iy = 1:length(yrs)
    ftmp1 = fullfile(ssmifol, fol1, sprintf('ssmr_%d.mat', yrs(iy)));
    ftmp2 = fullfile(ssmifol, fol1, sprintf('dmsp_%d.mat', yrs(iy)));
    
    if exist(ftmp1, 'file')
        fprintf('Loading %d: SSMI\n', yrs(iy));
        A = load(ftmp1);

        Ssmi.t{iy} = A.t';
        Ssmi.aice{iy} = local(A.aice,         grdmask, 'weight', Grd.area_feast, 'omitnan');
        Ssmi.gt10{iy} = local(A.aice>=thresh, grdmask, 'weight', Grd.area_feast, 'omitnan');
    end

    if exist(ftmp2, 'file')
        fprintf('Loading %d: DMSP\n', yrs(iy));
        A = load(ftmp2);

        Dmsp.t{iy} = A.t';
        Dmsp.aice{iy} = local(A.aice,         grdmask, 'weight', Grd.area_feast, 'omitnan');
        Dmsp.gt10{iy} = local(A.aice>=thresh, grdmask, 'weight', Grd.area_feast, 'omitnan');

    end
    % 
    % Obs.t{iy} = A.t';
    % Obs.aice{iy} = local(A.aice,         grdmask, 'weight', Grd.area_feast, 'omitnan');
    % Obs.gt10{iy} = local(A.aice>=thresh, grdmask, 'weight', Grd.area_feast, 'omitnan');
    % Obs.ltedge{iy} = nan(size(Obs.aice{iy}));
    % 
    % aiceline = reshape(A.aice, [], size(A.aice,3));
    % aiceline = aiceline(Line.idx,:);
    % for it = 1:size(aiceline,2)
    %     if any(aiceline(:,it) > thresh)
    %         [b,n,bi] = RunLength(aiceline(:,it)>= thresh);
    %         isice = b == 1 & n' > ncell; 
    %         if any(isice)
    %             idxedge = find(isice,1,'last');
    %             Obs.ltedge{iy}(it) = Grd.lat_rho(Line.idx(bi(idxedge)));
    %         end
    %     end
    % end
end

% isemp = cellfun(@isempty, Obs.t);

Ssmi = structfun(@(x) cat(1, x{:}), Ssmi, 'uni', 0);
Data.Ssmi = table2timetable(struct2table(Ssmi));

Dmsp = structfun(@(x) cat(1, x{:}), Dmsp, 'uni', 0);
Data.Dmsp = table2timetable(struct2table(Dmsp));


%% ... Models

sims = ["B10K-H16_CORECFS", "B10K-K20_CORECFS", "B10K-K20nobio_CORE", "B10K-K20nobio_CFS", "B10K-BGCMIPbanas_CFS"];
fld = strrep(sims, 'B10K-', '');    

for im = 1:length(sims)
    fprintf('%d/%d: %s\n', im, length(sims), sims(im));
    
    F = dir(fullfile(moxdir, 'roms_for_public', sims(im), 'Level1', '*average_aice.nc'));
    fname = fullfile({F.folder}, {F.name});

    Mod = struct;
    Mod.t = ncdateread(fname, 'ocean_time');

%     isin = Mod.t >= datetime(yr(1),1,1) & Mod.t < datetime(yr(end)+1,1,1);
%     Mod.t = Mod.t(isin);

    A = ncstruct(fname, 'aice');
%     A = ncstruct(fname, 'aice', struct('ocean_time', [find(isin,1) nnz(isin) 1]));

    Mod.aice = local(A.aice,         grdmask, 'weight', Grd.area_feast, 'omitnan');
    Mod.gt10 = local(A.aice>=thresh, grdmask, 'weight', Grd.area_feast, 'omitnan');
    Mod.ltedge = nan(size(Mod.aice));

    aiceline = reshape(A.aice, [], size(A.aice,3));
    aiceline = aiceline(Line.idx,:);
    for it = 1:size(aiceline,2)
        if any(aiceline(:,it) > thresh)
            [b,n,bi] = RunLength(aiceline(:,it)>= thresh);
            isice = b == 1 & n' > ncell; 
            if any(isice)
                idxedge = find(isice,1,'last');
                Mod.ltedge(it) = Grd.lat_rho(Line.idx(bi(idxedge)));
            end
        end
    end

    Data.(fld(im)) = table2timetable(struct2table(Mod));
    
end

save ice_comparison_data -struct Data;

return

%% Calculate climatologies

fld = fieldnames(Data);
nfld = length(fld);

for ii = 1:nfld
    isn = isnan(Data.(fld{ii}).ltedge);
    Data.(fld{ii}).ltedge(isn) = max(Line.lat);
end

nvar = size(Data.(fld{1}),2);

for ii = 1:nfld
    for iv = 1:nvar
        [Clim(ii,iv), Aligned(ii,iv)] = romsavgclimatology(Data.(fld{ii}){:,iv}, Data.(fld{ii}).t, 'realign', true, 'pivotmonth', 10);
    end
end

% For easier comparison, pad missing years in aligned arrays with NaNs

yr = unique(year(cat(1, Aligned.t)));
nyr = length(yr)-1;
nbin = size(Aligned(1).y,2);

for ii = 1:numel(Aligned)
    tmp = nan(nyr, nbin);
    [tf,loc] = ismember(year(Aligned(ii).t(:,1)), yr(1:end-1));
    tmp(loc,:) = Aligned(ii).y;
    Aligned(ii).y = tmp;
end


%% ... plot climatologies

h = plotgrid('size', [nvar 1]);
for iv = 1:nvar
    axes(h.ax(iv));
    
    [h.ln, h.p] = boundedline(Clim(1,iv).doy, cat(1, Clim(:,iv).mean)', permute(cat(1, Clim(:,iv).std), [2 3 1]), 'alpha');

end

h.leg = legendflex(h.ln, fld, 'interpreter', 'none');


%% Maps of max yearly ice

h = plotgrid('size', [2 2]);

% for ii = 1:4
%     axes(h.ax(ii));
%     plotromsrho(Grd, maxyearlyice(:,:,2,ii));
% end

iy = 10;

for ii = 1:4
    axes(h.ax(ii));
    plotromsrho(Grd,  mean(maxyearlyice(:,:,:,ii),3)- mean(maxyearlyice(:,:,:,1),3));
end

set(h.ax, 'clim', [-1 1], 'colormap', cmocean('balance'));

%% Hovmoller-type plots (obs w/ model diff from obs)

vars = Data.Obs.Properties.VariableNames;
vlong = {...
    'fraction EBS ice cover (aice mean)'
    'fraction EBS with aice>0.1'
    'latitude of main ice edge along 168.5W'};

yrall = unique(year(cat(1, Aligned.t)));
tmid = Aligned(1).t(1,1)+days(Clim(1).doy);

xday = days(tmid - tmid(1));
xtk = datetime(year(tmid(1)), month(tmid(1))+(0:12), 1);
xlbl = datestr(xtk, 'm');
xtk = days(xtk - xtk(1));

nyr = length(yrall);
nday = length(xday);

for iv = 1:nvar
    
    plt = cell(nfld,1);
    [plt{:}] = deal(nan(nyr,nday));
    
    for ii = 1:nfld
        [tf,loc] = ismember(year(Aligned(ii,iv).t(:,1)), yrall);
        plt{ii}(loc,:) = Aligned(ii,iv).y;
    end
    
    h = plotgrid('size', [2,3], 'sv', 0.1, 'mar', 0.05, 'mb', 0.18);
    
    for ii = 1:6
        axes(h.ax(ii));
        if ii == 1
            pcolorpad(xday, yrall, plt{ii});
        else
            pcolorpad(xday, yrall, plt{ii}-plt{1});
        end
        shading flat;
    end
    
    set(h.ax(2:end), 'clim', minmax([h.ax(2:end).CLim], 'center'), ...
                     'colormap', cmocean('balance'));
    
    set(h.ax, 'xlim', minmax(xtk), 'xtick', xtk, 'xticklabel', xlbl, ...
        'layer', 'top');
    multitextloc(h.ax(:), fld, 'northwestoutsideabove', 'interpreter', 'none');
    
    h.cb(1) = colorbar(h.ax(1), 'south');
    h.cb(2) = colorbar(h.ax(3), 'south');
    
    setpos(h.cb(1), '# 0.08 # #');
    setpos(h.cb(2), '# 0.08 # #');
    set(h.cb, 'axisloc', 'out');
    xlabel(h.cb(1), vlong{iv});

end
    
%% Skill stats

% Metric: southernmost yearly ice edge location

latmin = cell2mat(arrayfun(@(x) min(x.y,[],2), Aligned(:,3), 'uni', 0)');

% Date of advance and retreat

wknum = (1:52);

[wkadvance, wkretreat] = deal(nan(nyr,nfld));

for ii = 1:nfld
    for iy = 2:nyr-1
        if ~all(isnan(Aligned(ii,3).y(iy,:)))
            [xc,yc] = polyxpoly([wknum NaN], [Aligned(ii,3).y(iy,:) NaN], [0 53 NaN], [63 63 NaN]);
            wkadvance(iy,ii) = min(xc);
            wkretreat(iy,ii) = max(xc);
        end
    end
end

% Plot

marker = {'x', 'o', 's', '^', 'v', 'p'};
cmap = num2cell(cptcmap('Dark2_08'),2);
cmap = cmap(1:nfld);

metrics = {latmin wkadvance wkretreat};
metlabel = {'Southernmost latitude extent of ice', ...
            'Week of ice advance below 63N', ...
            'Week of ice retreat above 63N'};
nmet = length(metrics);

for im = 1:nmet
    S(im) = skillstats(metrics{im}(:,1), metrics{im}(:,2:end));
end


h = plotgrid('size', [nmet 2]);

for im = 1:nmet
    
    % Taylor
    
    h.tax(im) = tayloraxis(h.ax(im,1), 'npanel', 1, 'stdref', 1, 'stdmax', max(S(im).stdnorm));
    h.ln(im,:) = polarplot([real(acos(S(im).cor)); nan(1,nfld)], [S(im).stdnorm; nan(1,nfld)]);
    
    set(h.ln(im,:), {'marker', 'color'}, [marker' cmap], 'markersize', 10, 'linestyle', 'none');
%     set(h.tax(im).ax, 'rlim', [0 max(S(im).std)]);
    
    % Target
    
    h.tar(im,:) = plot(h.ax(im,2), [S(im).crmsd; nan(1,nfld)].*sign(S(im).std - S(im).std(1)), [S(im).bias; nan(1,nfld)]);
    set(h.tar(im,:), {'marker', 'color'}, [marker' cmap], 'markersize', 10, 'linestyle', 'none');
    ymax = max(abs(([S(im).crmsd S(im).bias])));
    
    set(h.ax(im,2), 'xlim', [-1 1]*ymax, 'ylim', [-1 1]*ymax);
    
end

multitextloc(h.ax(:,2), metlabel, 'northoutside');

set([h.tax.ax], 'clipping', 'off');
set(h.ax(:,2), 'xaxisloc', 'origin', 'yaxisloc', 'origin', 'dataaspectratio', [1 1 1], 'box', 'off');

legendflex(h.ln(1,:), fld, 'ref', h.fig, 'anchor', {'n','n'}, 'buffer', [0 0], ...
    'nrow', 2, 'interpreter', 'none');


%% Quick plot of recent observations

[xg1,yr1,tmid] = reshapetimeseries(Data.Ssmi.t, Data.Ssmi.aice);
[xg2,yr2] = reshapetimeseries(Data.Dmsp.t, Data.Dmsp.aice);

h1 = plot(tmid, xg1, 'color', rgb('light gray'));
hold on;
h2 = plot(tmid, xg2, 'color', rgb('light gray'));

set(h1(yr1==2023), 'color', rgb('orange'), 'linewidth', 2);
set(h2(yr2==2023), 'color', rgb('orange'), 'linewidth', 2);
set(h2(yr2==2024), 'color', rgb('red'), 'linewidth', 2);

set(h1(yr1==2018), 'color', rgb('green'), 'linewidth', 2);
set(h1(yr1==2019), 'color', rgb('blue'), 'linewidth', 2);








