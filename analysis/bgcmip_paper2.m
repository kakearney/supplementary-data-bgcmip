% bgcmip_paper2

% Simulation folder: [simfolder] = [moxdir]/kearney/BGC_hindcasts_workdir.
%
% Processing run prior to this script:
%
% [simfolder]/analysis/readbgcmetrics.m: saves metrics from specific
% station locations to tables, saved in nstation x nvariant structure Metrics with
% fields corresponding to the simulation names.
% variants are ["_CFS", "_CFS_loop", "bury_CFS", "bury_CFS_loop", "burynoinfauna_CFS", "burynoinfauna_CFS_loop"];
%
% [simfolder]/analysis/readnbudget.m:

% Setup

% simfolder = fullfile(moxdir, 'kearney', 'BGC_hindcasts_workdir');
simfolder = fullfile(mounteddir('klone'), 'GR011377_bgcmip', 'BGC_hindcasts_workdir');

run(fullfile(simfolder, 'analysis', 'bgcmip_setup'));

Metrics = load(fullfile(simfolder, 'analysis', 'bgcmip_metrics_at_stations.mat'));
Metrics = Metrics.Metrics;

bgcname = fieldnames(Metrics);
nbgc = length(bgcname);

%% ... colors

% GeoDataViz qualitative palette

colgeodataviz = ["ff1f5b" "00cd6c" "009ade" "af58ba" "ffc61e" "f28522" "a0b1ba" "a6761d"];

colors = [...
 "purple"   colgeodataviz(4)
 "green"    colgeodataviz(2)
 "yellow"   colgeodataviz(5) 
 "brown"    colgeodataviz(8) % actually brown
 "red"      colgeodataviz(1)
 "gray"     colgeodataviz(7)
 "dark"     "000000"       
 "orange"   colgeodataviz(6)
 "blue"     colgeodataviz(3) % actually blue
 ];

colors = table(colors(:,2), 'rownames', colors(:,1));
colors.Properties.VariableNames = "hex";

tmp = arrayfun(@hex2rgb, colors.hex, 'uni', 0);
colors.rgb = cat(1, tmp{:});

colmodel = colors{["green", "blue", "purple"], "rgb"};

adjustcol = @(col,x) interp1([-1 0 1], [0 0 0; col; 1 1 1], x);

%% ... station reordering

nstat = height(Station);

% statordergrd = [1  7 11   9 NaN; ...  % Alaska peninsula, outer to inner
%                 3 12 16  15  14; ...  % Outer shelf, south to north
%                 2 10  4   5   6; ...  % Middle shelf, south to north
%                 8 13 17 NaN NaN; ...  % Inner shelf, south to north;
%                21 18 19  20 NaN];

statordergrd = [21 18 19  20; ...
               14   5 17 NaN; ...
               16   4 13 NaN; ...
                3   2  8 NaN; ...
                1   7 11   9];


statorder = statordergrd';
statorder = statorder(~isnan(statorder));

statletter = string(char((1:nstat)+96)');

[tf,loc] = ismember(1:nstat, statorder);
Station.letter(tf) = statletter(loc(tf));

%% Map: Grid, stations, etc.

%% ... Borders

[blat{1}, blon{1}] = borders('alaska');
[blat{2}, blon{2}] = borders('canada');
[blat{3}, blon{3}] = borders('russia');

blon = cellfun(@wrapTo360, blon, 'uni', 0);

% Polygon (for shaded)

Border.poly = polyshape(blon, blat);

% Lines (for outline w/o 180 line)

is180 = abs(blon{3} - 180) < 0.01;
[b,n,bi] = RunLength(is180);
n = n(b);
bi = bi(b);
for ii = 1:length(bi)
    if n(ii)>2
        blon{3}(bi(ii)+(1:(n(ii)-2))) = NaN;
        blat{3}(bi(ii)+(1:(n(ii)-2))) = NaN;
    end
end
needextra = bi(n==2);
nper = diff([0 needextra length(blat{3})]);
blat(3+(0:length(needextra))) = mat2cell(blat{3}, 1, nper);
blon(3+(0:length(needextra))) = mat2cell(blon{3}, 1, nper);

Border.lon = [blon{:}];
Border.lat = [blat{:}];

%% ... data setup

latlim = minmax(Grd.lat_rho);
lonlim = minmax(Grd.lon_rho);

htmp = figure;
worldmap(latlim, lonlim);
proj = getm(gca);
close(htmp);

% Project stuff

[P.xgrd, P.ygrd] = projfwd(proj, Grd.lat_psi, Grd.lon_psi);

xlim = minmax(P.xgrd, 'expand');
ylim = minmax(P.ygrd, 'expand');

xybox = interparc(linspace(0,1,200), xlim([1 1 2 2 1]), ylim([1 2 2 1 1]), 'linear');
[ltbox, lnbox] = projinv(proj, xybox(:,1), xybox(:,2));

ltextract = minmax(ltbox);
lnextract = minmax(wrapTo360(lnbox));

% ETOPO_2022 data

tilename = 'https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/60s/60s_surface_elev_netcdf/ETOPO_2022_v1_60s_N90W180_surface.nc';
E = ncstruct(tilename, 'lat', 'lon');
latmask = E.lat >= ltextract(1) & E.lat <= ltextract(2);
lonmask = wrapTo360(E.lon) >= lnextract(1) & wrapTo360(E.lon) <= lnextract(2);

latscs = [find(latmask, 1, 'first') nnz(latmask) 1];
[b, n, bi] = RunLength(lonmask);
lonscs = [bi(b==1) n(b==1) ones(nnz(b),1)];

Etop(1) = ncstruct(tilename, struct('lat', latscs, 'lon', lonscs(1,:)));
Etop(2) = ncstruct(tilename, struct('lat', latscs, 'lon', lonscs(2,:)));

Etop(1).lon = wrapTo360(cat(1, Etop([2 1]).lon));
Etop(1).z   = cat(1, Etop([2 1]).z);
Etop = Etop(1);

[lntop, lttop] = ndgrid(Etop.lon, Etop.lat);
[P.xtop, P.ytop] = projfwd(proj, lttop, lntop);

% Borders, stations, polygons projected

[P.xbor, P.ybor] = projfwd(proj, Border.lat, Border.lon);

[P.xstat, P.ystat] = projfwd(proj, Station.lat, Station.lon);

Bsierp = shapeprjread('/Volumes/LaCie2023/AlaskaShapefiles/bsierp/BSIERP_regions_2012.shp');
[P.xbsierp, P.ybsierp] = projfwd(proj, [Bsierp.Lat], [Bsierp.Lon]);

Strata = catstruct(1, shapeprjread('/Volumes/LaCie2023/AlaskaShapefiles/gis_updated/EBS_NBS_2019.shp'), ...
                      shapeprjread('/Volumes/LaCie2023/AlaskaShapefiles/gis_updated/RUSSIAN_2019.shp'));
[P.xstrata, P.ystrata] = projfwd(proj, [Strata.Lat], [Strata.Lon]);
issebs = [Strata.STRATUM] <= 62;
pebs = polyshape([Strata(issebs).Lon], [Strata(issebs).Lat]);
[P.xebs, P.yebs] = projfwd(proj, pebs.Vertices(:,2), pebs.Vertices(:,1));

%% ... plot

% Create figure

h = plotgrid('size', [1 1], 'margin', 0.01);
setpos(h.fig, '# # 600 600');
set(h.fig, 'color', 'w');

% ETOPO shaded relief

h.ax(2) = copyobj(h.ax(1), h.fig);
axes(h.ax(1));

h.topo = pcolor(P.xtop', P.ytop', Etop.z');
colormap(adjustcol(colors{'gray','rgb'}, 0.5));
set(h.ax, 'xlim', minmax(P.xgrd), 'ylim', minmax(P.ygrd), 'dataaspectratio', [1 1 1]);
shadem([90 45], 4);
shading flat;

% B10K grid

axes(h.ax(2));

h.roms = pcolor(P.xgrd, P.ygrd, Grd.h(2:end,2:end));
shading flat;
set(h.roms, 'facealpha', 0.6);
cmapbathy = crameri('-oslo'); %cmocean('deep'); %adjustcol(colors{'purple','rgb'}, linspace(1,0,50));
set(h.ax(2), 'colormap', cmapbathy, 'clim', [0 7000]);
hold(h.ax(2), 'on');

[cc, h.cont] = contour(P.xtop', P.ytop', Etop.z', [-50 -100 -200], 'k');
set(h.cont, 'color', rgb('gray'));

% Borders

plot(P.xbor, P.ybor, 'color', 'k');

% Stations and BSIERP regions

dx = diff(minmax(P.xgrd))*0.01;

h.letter = text(P.xstat-dx, P.ystat, Station.letter, 'horiz', 'right', ...
    'fontname', 'menlo', 'fontsize', 12);
ismoor = startsWith(Station.name, 'BS-');
h.stat1 = plot(P.xstat(ismoor), P.ystat(ismoor), 'r^');
h.stat2 = plot(P.xstat(~ismoor), P.ystat(~ismoor), 'ro');
set([h.stat1 h.stat2], 'markerfacecolor', 'w', 'markersize', 5, 'markeredgecolor', 'k');

h.bsebs = plot(P.xebs, P.yebs, 'color', adjustcol(colors{'red','rgb'}, 0.5), 'linewidth', 1.5);
h.breg = plot(P.xbsierp, P.ybsierp, 'color', colors{'yellow','rgb'}, 'linestyle', '-');
h.bstrata = plot(P.xstrata, P.ystrata, 'color', colors{'red','rgb'}, 'linestyle', ':');


uistack(h.letter, 'top');

% Axes

set(h.ax(2), 'visible', 'off');
set(h.ax, 'xlim', minmax(P.xgrd), 'ylim', minmax(P.ygrd), 'dataaspectratio', [1 1 1], ...
    'layer', 'top', 'xtick', [], 'ytick', []);

% Colorbar

% drawnow;
% puase(1);
h.cb = colorbar('south');
h.cb.Position(3) = 0.3;
% setpos(h.cb, '0.1 0.2 0.3 #');
set(h.cb, 'ticks', [], 'edgecolor', colors{'dark','rgb'});
tprops = {'edgecolor', 'none', 'fontsize', 12, 'vert', 'middle'};
annotation('textbox', h.cb.Position, 'string', '0m', 'horiz', 'left', tprops{:});
annotation('textbox', h.cb.Position, 'string', '7000m', 'horiz', 'right', tprops{:}, 'color', 'w');

% Inset map

pos = plotboxpos(h.ax(2));
h.axin = axes('position', [pos(1)+pos(3)-0.2 pos(2)+pos(4)-0.2, 0.2, 0.2]);
globefill;
hold on;
C = load('coastlines');
globeplot(C.coastlat, C.coastlon, 'color', 'k');
globegraticule('linestyle', ':');
axis tight;
view(-90,50);
globeplot(ltbox, lnbox, 'color', 'k', 'linewidth', 1.5);


export_fig(h.fig, 'B10K_map', '-png', '-r150');


%% ... EPOC slide version

set(h.breg, 'visible', 'off');
set(h.stat1, 'visible', 'off');
set(h.stat2, 'visible', 'off');
set(h.letter, 'visible', 'off');

export_fig(h.fig, '~/Documents/Conferences/202409_EPOC/EPOC_Kearney_slides/B10K_domain', '-png', '-r150');


%% NO3/NH4 timeseries: How quickly does nutrient environment diverge between models?

nstat = height(Station);
% models = ["banas" "bestnpz" "cobalt"];
% nmodel = length(models);

h = plotgrid('size', [nstat 2], 'sh', 0.05, 'mar', 0.05, 'sv', 0, 'ml', 0.08);
setpos(h.fig, '# # 8.5in 11in');
arrayfun(@(x) hold(x,'on'), h.ax);
set(h.fig, 'color', 'w');


Mplt = Metrics(statorder,1);


for is = 1:nstat
    for im = 1:nbgc
        h.ln(is,im,1) = plot(h.ax(is,1), Mplt(is).(bgcname{im}), 'NO3b', 'color', colmodel(im,:));
        h.ln(is,im,2) = plot(h.ax(is,2), Mplt(is).(bgcname{im}), 'NH4b', 'color', colmodel(im,:));
    end
end

set(h.ax(:,1), 'ylim', [0 40], 'xlim', minmax(Mplt(1).banas.Time));
set(h.ax(:,2), 'ylim', [0 20], 'xlim', minmax(Mplt(1).banas.Time));
set(h.ax, 'clipping', 'off', 'box', 'off', 'fontsize', 8, 'tickdir', 'out');
set(h.ax(1:end-1,:), 'xcolor', 'none');
set(h.ax(1:2:end,:), 'color', ones(1,3)*0.95);
set(h.ax(1:2:end,:), 'yaxisloc', 'right');
arrayfun(@(x) set(get(x, 'ylabel'), 'string', ''), h.ax);

h.lbl = labelaxes(h.ax(:,1),statletter+")", 'westoutside', ...
    'hbuffer', 0.08, 'fontweight', 'b', 'fontname', 'menlo');

h.ttl = labelaxes(h.ax(1,:), {'Bottom NO_3 (mmol m^{-3})', 'Bottom NH_4 (mmol m^{-3})'}, 'northoutside', ...
    'fontweight', 'b');


%% Climatological metrics by station: Do they agree across models?

%% ... Climatological metrics by station, bgc model

mvars = Metrics(1).banas.Properties.VariableNames;
nmetric = length(mvars);

MetClim.avg = nan(52,nbgc+1,nmetric,nstat);
MetClim.all = nan(52,nbgc+1,nmetric,nstat,30);

% Primary variants (loop 2) of each BGC model

for ii = 1:nbgc
    for im = 1:nmetric
        for is = 1:nstat
            t = Metrics(is,2).(bgcname{ii}).Time;
            y = Metrics(is,2).(bgcname{ii}).(mvars{im});
            if ~all(isnan(y))
                [Clim, Aligned] = romsavgclimatology(y,t);
                MetClim.avg(:,ii,im,is) = Clim.mean;
                MetClim.all(:,ii,im,is,:) = permute(Aligned.y, [2 3 4 5 1]);
            end
        end
    end
end

MetClim.doy = Clim.doy;

% Add +burial -infauna bestnpz version

for im = 1:nmetric
    for is = 1:nstat
        t = Metrics(is,6).bestnpz.Time;
        y = Metrics(is,6).bestnpz.(mvars{im});
        if ~all(isnan(y))
            [Clim, Aligned] = romsavgclimatology(y,t);
            MetClim.avg(:,nbgc+1,im,is) = Clim.mean;
            MetClim.all(:,nbgc+1,im,is,:) = permute(Aligned.y, [2 3 4 5 1]);
        end
    end
end

%% ...... cluster stations with similar metrics

mclus = reshape(metricavg(:,:,8:9,:), [], nstat)';

z = linkage(mclus, 'ward');
ctmp = clusterdata(mclus, ...
        'linkage', 'ward', 'distance', 'euclidean', 'criterion', ...
        'distance', 'cutoff', 0.2 * max(z(:,3)));
cmapcluster = cptcmap('d3cat20');


figure;
[hh,t,p] = dendrogram(z,0, ...
            'orientation', 'right', 'colorthreshold', 0.2 * max(z(:,3)));

figure;

plot(wrapTo360([Bsierp.Lon]), [Bsierp.Lat], 'k');
hold on;
scatter(wrapTo360(Station.lon), Station.lat, 100, ctmp, 'filled');
text(wrapTo360(Station.lon), Station.lat, compose('%d--', 1:nstat), 'horiz', 'right')
text(wrapTo360(Station.lon), Station.lat, Station.name, 'horiz', 'left');
set(gca, 'clim', [0.5 20.5], 'colormap', cmapcluster);

%% ... one figure per metric

mvars = Metrics(1).banas.Properties.VariableNames;

pltmask = isnan(statordergrd);

lprops = {'linewidth', 1.5};

for im = 8:9 %1:length(mvars)

    h = plotgrid('size', size(statordergrd), 'sp', 0.01, 'mar', 0.05);
    
    for ii = 1:numel(statordergrd)
        if ~pltmask(ii)
            is = statordergrd(ii);

            for ibgc = 1:3
                plot(h.ax(ii), MetClim.doy, MetClim.avg(:,ibgc,im,is), 'color', colmodel(ibgc,:), lprops{:});
                hold(h.ax(ii), 'on');
            end

            % 
            % if all(isnan(Metrics(is,1).banas.(mvars{im})))
            %     plot(h.ax(ii), NaN, NaN, 'color', colmodel(1,:));
            % else
            %     Clim = romsavgclimatology(Metrics(is,1).banas.(mvars{im}), Metrics(is,1).banas.Time);
            %     plot(h.ax(ii), Clim.doy, Clim.mean, 'color', colmodel(1,:), lprops{:});
            % end
            % hold(h.ax(ii), 'on');
            % 
            % Clim = romsavgclimatology(Metrics(is,1).bestnpz.(mvars{im}), Metrics(is,1).banas.Time);
            % plot(h.ax(ii), Clim.doy, Clim.mean, 'color', colmodel(2,:), lprops{:});
            % 
            % Clim = romsavgclimatology(Metrics(is,1).cobalt.(mvars{im}), Metrics(is,1).banas.Time);
            % plot(h.ax(ii), Clim.doy, Clim.mean, 'color', colmodel(3,:), lprops{:});
        end
    end

    set(h.ax(pltmask), 'visible', 'off');
    set(h.ax, 'ylim', minmax(cat(1, h.ax(~pltmask).YLim)), 'box', 'off', 'xlim', [0 365]);
    set(h.ax(:,2:end), 'yticklabel', '');
    set(h.ax(1:end-1,:), 'xticklabel', '');
    set(h.fig, 'color', 'w', 'name', mvars{im});

    h.lbl = labelaxes(h.ax(~pltmask), Station.letter(statordergrd(~pltmask)), 'northwest', 'vbuffer', 0, 'hbuffer', 0.05);
    set(h.lbl, {'color'}, num2cell(cmapcluster(ctmp(statordergrd(~pltmask)),:),2), 'fontweight', 'b');

end

%% Metrics in north/south, cross-shelf axes, EPOC slide version

%% ... Axes

pltmask = isnan(statordergrd);
lprops = {'linewidth', 1.5};

h = plotgrid('size', size(statordergrd), 'sp', 0.01, 'mar', 0.05);
h.fig.Position(3:4) = [800 600];
h.lbl = labelaxes(h.ax(~pltmask), Station.letter(statordergrd(~pltmask)), 'northwest', 'vbuffer', 0, 'hbuffer', 0.05, 'fontname', 'menlo', 'fontsize', 18, 'fontweight', 'b');
xtk = datetime(2024,1:12,1);
set(h.ax, 'ylim', minmax(cat(1, h.ax(~pltmask).YLim)), 'box', 'off', 'xlim', [0 365], 'xtick', doy(xtk), 'xticklabel', datestr(xtk, 'm'));
set(h.ax, 'yticklabel', '');
set(h.ax(1:end-1,:), 'xticklabel', '');
set(h.fig, 'color', 'w');
set(h.ax(pltmask), 'visible', 'off');
set(h.ax(~pltmask), 'linewidth', 1.0);

set(h.lbl(ismember(string({h.lbl.String}), ["e","h","k"])), 'color', colors{'orange','rgb'});
set(h.lbl(ismember(string({h.lbl.String}), ["f","i","l"])), 'color', colors{'yellow','rgb'});
set(h.lbl(ismember(string({h.lbl.String}), ["g","j","m"])), 'color', colors{'red','rgb'});
set(h.lbl(ismember(string({h.lbl.String}), ["a","b","c","d"])), 'color', colors{'gray','rgb'});
set(h.lbl(ismember(string({h.lbl.String}), ["n","o","p","q"])), 'color', colors{'brown','rgb'});

export_fig(h.fig, '~/Documents/Conferences/202409_EPOC/EPOC_Kearney_slides/metrics_axes', '-png', '-r150', '-nocrop');
close(h.fig);

%% ... NPP

im = find(strcmp(mvars, 'npp'));

h = plotgrid('size', size(statordergrd), 'sp', 0.01, 'mar', 0.05);
h.fig.Position(3:4) = [800 600];

for ii = 1:numel(statordergrd)
    if ~pltmask(ii)
        is = statordergrd(ii);

        for ibgc = 1:3
            h.ln(ii,ibgc) = plot(h.ax(ii), MetClim.doy, MetClim.avg(:,ibgc,im,is), 'color', colmodel(ibgc,:), lprops{:});
            hold(h.ax(ii), 'on');
        end
        h.ln2(ii) = plot(h.ax(ii), MetClim.doy, MetClim.avg(:,nbgc+1,im,is), 'color', brighten(colmodel(2,:), -0.5), lprops{:});
    end
end

h.lbl = labelaxes(h.ax(~pltmask), Station.letter(statordergrd(~pltmask)), 'northwest', 'vbuffer', 0, 'hbuffer', 0.05, 'fontname', 'menlo', 'fontsize', 18, 'fontweight', 'b');
set(h.ax, 'ylim', minmax(cat(1, h.ax(~pltmask).YLim)), 'box', 'off', 'xlim', [0 365], 'xtick', doy(xtk), 'xticklabel', datestr(xtk, 'm'));
set(h.ax(:,2:end), 'yticklabel', '');
set(h.ax(1:end-1,:), 'xticklabel', '');
set(h.fig, 'color', 'w');
set(h.ax(pltmask), 'visible', 'off');
set(h.ax(~pltmask), 'linewidth', 1.0);

set(h.lbl(ismember(string({h.lbl.String}), ["e","h","k"])), 'color', colors{'orange','rgb'});
set(h.lbl(ismember(string({h.lbl.String}), ["f","i","l"])), 'color', colors{'yellow','rgb'});
set(h.lbl(ismember(string({h.lbl.String}), ["g","j","m"])), 'color', colors{'red','rgb'});
set(h.lbl(ismember(string({h.lbl.String}), ["a","b","c","d"])), 'color', colors{'gray','rgb'});
set(h.lbl(ismember(string({h.lbl.String}), ["n","o","p","q"])), 'color', colors{'brown','rgb'});

legprop = {'fontsize', 20, 'fontweight', 'b', 'hbuffer', 0.2, 'interpreter', 'none'};
h.leg(1) = labelaxes(h.ax(3,end), "Banas", 'northwest', 'color', colmodel(1,:), legprop{:});
h.leg(2) = labelaxes(h.ax(3,end), "BEST_NPZ",   'west', 'color', colmodel(2,:), legprop{:});
h.leg(3) = labelaxes(h.ax(3,end), "COBALT",'southwest', 'color', colmodel(3,:), legprop{:});

set(h.ln2(~pltmask), 'visible', 'off');
export_fig(h.fig, '~/Documents/Conferences/202409_EPOC/EPOC_Kearney_slides/metrics_npp', '-png', '-r150', '-nocrop');
set(h.ln(~pltmask,1), 'linewidth', 3);
export_fig(h.fig, '~/Documents/Conferences/202409_EPOC/EPOC_Kearney_slides/metrics_npp_banas', '-png', '-r150', '-nocrop');
set(h.ln(~pltmask,:), 'linewidth', 1.5);
set(h.ln(~pltmask,2), 'linewidth', 3);
export_fig(h.fig, '~/Documents/Conferences/202409_EPOC/EPOC_Kearney_slides/metrics_npp_bestnpz', '-png', '-r150', '-nocrop');
set(h.ln(~pltmask,:), 'linewidth', 1.5);
set(h.ln(~pltmask,3), 'linewidth', 3);
export_fig(h.fig, '~/Documents/Conferences/202409_EPOC/EPOC_Kearney_slides/metrics_npp_cobalt', '-png', '-r150', '-nocrop');

set(h.ln2(~pltmask), 'visible', 'on', 'linewidth', 3);
set(h.ln(~pltmask,:), 'linewidth', 1.5);
export_fig(h.fig, '~/Documents/Conferences/202409_EPOC/EPOC_Kearney_slides/metrics_npp_bestnpzinfauna', '-png', '-r150', '-nocrop');

close(h.fig);

%% ... NO3

im1 = find(strcmp(mvars, 'NO3s'));
im2 = find(strcmp(mvars, 'NO3b'));

h = plotgrid('size', size(statordergrd), 'sp', 0.01, 'mar', 0.05);
h.fig.Position(3:4) = [800 600];

for ii = 1:numel(statordergrd)
    if ~pltmask(ii)
        is = statordergrd(ii);

        for ibgc = 1:3
            h.ln(ii,ibgc,1) = plot(h.ax(ii), MetClim.doy, MetClim.avg(:,ibgc,im1,is), 'color', colmodel(ibgc,:), lprops{:});
            hold(h.ax(ii), 'on');
            h.ln(ii,ibgc,2) = plot(h.ax(ii), MetClim.doy, MetClim.avg(:,ibgc,im2,is), 'color', colmodel(ibgc,:), lprops{:}, 'linestyle', '-.');
        end
    end
end

h.lbl = labelaxes(h.ax(~pltmask), Station.letter(statordergrd(~pltmask)), 'northwest', 'vbuffer', 0, 'hbuffer', 0.05, 'fontname', 'menlo', 'fontsize', 18, 'fontweight', 'b');
set(h.ax, 'ylim', minmax(cat(1, h.ax(~pltmask).YLim)), 'box', 'off', 'xlim', [0 365], 'xtick', doy(xtk), 'xticklabel', datestr(xtk, 'm'));
set(h.ax(:,2:end), 'yticklabel', '');
set(h.ax(1:end-1,:), 'xticklabel', '');
set(h.fig, 'color', 'w');
set(h.ax(pltmask), 'visible', 'off');
set(h.ax(~pltmask), 'linewidth', 1.0);

set(h.lbl(ismember(string({h.lbl.String}), ["e","h","k"])), 'color', colors{'orange','rgb'});
set(h.lbl(ismember(string({h.lbl.String}), ["f","i","l"])), 'color', colors{'yellow','rgb'});
set(h.lbl(ismember(string({h.lbl.String}), ["g","j","m"])), 'color', colors{'red','rgb'});
set(h.lbl(ismember(string({h.lbl.String}), ["a","b","c","d"])), 'color', colors{'gray','rgb'});
set(h.lbl(ismember(string({h.lbl.String}), ["n","o","p","q"])), 'color', colors{'brown','rgb'});

legprop = {'fontsize', 20, 'fontweight', 'b', 'hbuffer', 0.2, 'interpreter', 'none'};
h.leg(1) = labelaxes(h.ax(3,end), "Banas", 'northwest', 'color', colmodel(1,:), legprop{:});
h.leg(2) = labelaxes(h.ax(3,end), "BEST_NPZ",   'west', 'color', colmodel(2,:), legprop{:});
h.leg(3) = labelaxes(h.ax(3,end), "COBALT",'southwest', 'color', colmodel(3,:), legprop{:});

export_fig(h.fig, '~/Documents/Conferences/202409_EPOC/EPOC_Kearney_slides/metrics_no3', '-png', '-r150', '-nocrop', '-painters');
close(h.fig);

%% ... NH4

im1 = find(strcmp(mvars, 'NH4s'));
im2 = find(strcmp(mvars, 'NH4b'));

h = plotgrid('size', size(statordergrd), 'sp', 0.01, 'mar', 0.05);
h.fig.Position(3:4) = [800 600];


for ii = 1:numel(statordergrd)
    if ~pltmask(ii)
        is = statordergrd(ii);

        for ibgc = 1:3
            h.ln(ii,ibgc,1) = plot(h.ax(ii), MetClim.doy, MetClim.avg(:,ibgc,im1,is), 'color', colmodel(ibgc,:), lprops{:});
            hold(h.ax(ii), 'on');
            h.ln(ii,ibgc,2) = plot(h.ax(ii), MetClim.doy, MetClim.avg(:,ibgc,im2,is), 'color', colmodel(ibgc,:), lprops{:}, 'linestyle', '-.');
        end
    end
end

h.lbl = labelaxes(h.ax(~pltmask), Station.letter(statordergrd(~pltmask)), 'northwest', 'vbuffer', 0, 'hbuffer', 0.05, 'fontname', 'menlo', 'fontsize', 18, 'fontweight', 'b');
set(h.ax, 'ylim', minmax(cat(1, h.ax(~pltmask).YLim)), 'box', 'off', 'xlim', [0 365], 'xtick', doy(xtk), 'xticklabel', datestr(xtk, 'm'));
set(h.ax(:,2:end), 'yticklabel', '');
set(h.ax(1:end-1,:), 'xticklabel', '');
set(h.fig, 'color', 'w');
set(h.ax(pltmask), 'visible', 'off');
set(h.ax(~pltmask), 'linewidth', 1.0);

set(h.lbl(ismember(string({h.lbl.String}), ["e","h","k"])), 'color', colors{'orange','rgb'});
set(h.lbl(ismember(string({h.lbl.String}), ["f","i","l"])), 'color', colors{'yellow','rgb'});
set(h.lbl(ismember(string({h.lbl.String}), ["g","j","m"])), 'color', colors{'red','rgb'});
set(h.lbl(ismember(string({h.lbl.String}), ["a","b","c","d"])), 'color', colors{'gray','rgb'});
set(h.lbl(ismember(string({h.lbl.String}), ["n","o","p","q"])), 'color', colors{'brown','rgb'});

legprop = {'fontsize', 20, 'fontweight', 'b', 'hbuffer', 0.2, 'interpreter', 'none'};
h.leg(1) = labelaxes(h.ax(3,end), "Banas", 'northwest', 'color', colmodel(1,:), legprop{:});
h.leg(2) = labelaxes(h.ax(3,end), "BEST_NPZ",   'west', 'color', colmodel(2,:), legprop{:});
h.leg(3) = labelaxes(h.ax(3,end), "COBALT",'southwest', 'color', colmodel(3,:), legprop{:});

export_fig(h.fig, '~/Documents/Conferences/202409_EPOC/EPOC_Kearney_slides/metrics_nh4', '-png', '-r150', '-nocrop', '-painters');
close(h.fig);


%% Metrics in north/south, cross-shelf axes, paper version

% Plots

pltvars = {...
     "npp",         'NPP (mmol N m^{-2})'   60  
    ["NO3s" "NO3b"] 'NO_3 (mmol N m^{-3})' NaN
    ["NH4s" "NH4b"] 'NH_4 (mmol N m^{-3})'      NaN
     "zprod_micro"  'Microzoo. prod. (mmol N m^{-2})' 40
     "zprod_meso"   'Mesozoo. prod. (mmol N m^{-2})' NaN
     "plfrac"       'Fraction large phyto' 1
     "zmicro"       'Microzoo. biomass (mmol N m^{-2})' NaN
     "zmeso"        'Mesozoo. biomass (mmol N m^{-2})' NaN
    };

pltvars = cell2table(pltvars, 'variablenames', {'vars', 'lname', 'maxval'});

% pltmask = isnan(statordergrd);
% lprops = {'linewidth', 1.5};

xtk = datetime(2024,1:12,1);


% Loop over variables to plot

lineprops = {'linewidth', 1.5};
legprop = {'fontsize', 12, 'fontweight', 'b', 'hbuffer', 0.2, 'interpreter', 'none'};
legprop = {'anchor', {'w','w'}, 'buffer', [0 0], 'interpreter', 'none', 'box', 'off', 'fontsize', 12};
labelprop = {'vbuffer', 0, 'hbuffer', 0.05, 'fontname', 'menlo', 'fontsize', 12};

colplt = [colmodel; colors{'yellow','rgb'}];

for iplt = 1:height(pltvars)

    [~,im] = ismember(pltvars.vars{iplt}, mvars);
    % im = find(strcmp(mvars, 'npp'));
    
    h = plotgrid('size', size(statordergrd), 'sp', 0.015, 'mar', 0.05);
    h.fig.Position(3:4) = [800 600];
    
    for ii = 1:numel(statordergrd)
        if ~pltmask(ii)
            is = statordergrd(ii);
    
            if length(im) == 1
                [~,h.ens(ii)] = ensemble2bnd(MetClim.doy, permute(MetClim.all(:,1:4,im,is,:), [1 2 5 3 4]), 'prc', [10 90], 'plot', 'boundedline', 'alpha', true, 'axis', h.ax(ii), 'cmap', colplt, 'tlim', [0.2 0.2]);
                set(h.ens(ii).ln, lineprops{:});
            elseif length(im) == 2
                [~,h.ens(ii,1)] = ensemble2bnd(MetClim.doy, permute(MetClim.all(:,1:4,im(1),is,:), [1 2 5 3 4]), 'prc', [10 90], 'plot', 'boundedline', 'alpha', true, 'axis', h.ax(ii), 'cmap', colplt, 'tlim', [0.2 0.2]);
                set(h.ens(ii,1).ln, lineprops{:});
                [~,h.ens(ii,2)] = ensemble2bnd(MetClim.doy, permute(MetClim.all(:,1:4,im(2),is,:), [1 2 5 3 4]), 'prc', [10 90], 'plot', 'boundedline', 'alpha', true, 'axis', h.ax(ii), 'cmap', colplt, 'tlim', [0.2 0.2]);
                set(h.ens(ii,2).ln, lineprops{:}, 'linestyle', '-.');
            end
        end
    end
    % set(h.ens, 'edgealpha', 0.1);
    
    h.lbl = labelaxes(h.ax(~pltmask), Station.letter(statordergrd(~pltmask)), 'northwest', labelprop{:});
    set(h.ax, 'ylim', minmax(cat(1, h.ax(~pltmask).YLim)), 'box', 'off', 'xlim', [0 365], 'xtick', doy(xtk), 'xticklabel', datestr(xtk, 'm'));
    if ~isnan(pltvars.maxval(iplt))
        set(h.ax, 'ylim', [0 pltvars.maxval(iplt)]);
    end
    set(h.ax(:,2:end), 'yticklabel', '');
    set(h.ax(1:end-1,:), 'xticklabel', '');
    set(h.fig, 'color', 'w');
    set(h.ax(pltmask), 'visible', 'off');
    
    if length(im) == 1
        h.leg = legendflex(h.ens(1).ln, {'Banas', 'BEST_NPZ', 'COBALT', 'BEST_NPZ no-ben'}, 'ref', h.ax(3,end), legprop{:});
    else
        h.leg = legendflex([h.ens(1,:).ln], {'Banas, surface', 'BEST_NPZ, surface', 'COBALT, surface', 'BEST_NPZ no-ben, surface', ...
                                             'Banas, bottom',  'BEST_NPZ, bottom',  'COBALT, bottom',  'BEST_NPZ no-ben, bottom'}, 'ref', h.ax(3,end), legprop{:});
    end

    ylabel(h.ax(1,1), pltvars.lname{iplt});
    
    % export_fig(h.fig, "metrics_"+pltvars.vars{iplt}(1), '-png', '-r150', '-painters');
    % close(h.fig);

end

%% ... same axes setup, but NO3,NH4 from initialization

h = plotgrid('size', [2 1], 'mar', 0.05, 'mr', 0.01);
setpos(h.fig, '# # 8.5in 11in');
set(h.fig, 'color', 'w');
sz = size(statordergrd);
h.ax = subgridaxes(h.ax, sz(1), sz(2));


for ii = 1:numel(statordergrd)
    if ~pltmask(ii)
        is = statordergrd(ii);

        for im = 1:nbgc
            h.ln(ii,im,1) = plot(h.ax(ii), Metrics(is).(bgcname{im}), 'NH4b', 'color', colmodel(im,:), 'linewidth', 0.1);
            hold(h.ax(ii), 'on');
            h.ln(ii,im,2) = plot(h.ax(ii+prod(sz)), Metrics(is).(bgcname{im}), 'NO3b', 'color', colmodel(im,:), 'linewidth', 0.1);
            hold(h.ax(ii+prod(sz)), 'on');
            % h.ln(ii,im,3) = plot(h.ax(ii), Metrics(is).(bgcname{im}), 'NH4s', 'color', adjustcol(colmodel(im,:), 0.5));
            % h.ln(ii,im,3) = plot(h.ax(ii), Metrics(is).(bgcname{im}), 'NH4b', 'color', adjustcol(colmodel(im,:), 0.5), 'linestyle', '-.');
        end
    end
end
set(h.ax(:,:,1), 'ylim', [0 20]); 
set(h.ax(:,:,2), 'ylim', [0 40]); 
set(h.ax, 'tickdir', 'out', 'fontsize', 8);
set(h.fig, 'color', 'w');

set(h.ax(cat(3, pltmask, pltmask)), 'visible', 'off');

for ii = 1:numel(h.ax)
    h.ax(ii).Position(3) = 0.9*h.ax(ii).Position(3);
    h.ax(ii).Position(4) = 0.9*h.ax(ii).Position(4);
    h.ax(ii).XLabel.String = '';
    h.ax(ii).YLabel.String = '';
end
ylabel(h.ax(1,1,1), 'Bottom NH_4 (mmol N m^{-3})');
ylabel(h.ax(1,1,2), 'Bottom NO_3 (mmol N m^{-3})');

tmp = h.ax(:,:,1);
h.lbl = labelaxes(tmp(~pltmask), Station.letter(statordergrd(~pltmask)), 'northwest', labelprop{:});
set(tmp(:,2:end), 'yticklabel', '');
set(tmp(1:end-1,1:end-1), 'xticklabel', '');
tmp = h.ax(:,:,2);
h.lbl = labelaxes(tmp(~pltmask), Station.letter(statordergrd(~pltmask)), 'northwest', labelprop{:});
set(tmp(:,2:end), 'yticklabel', '');
set(tmp(1:end-1,1:end-1), 'xticklabel', '');

export_fig(h.fig, 'bottom_nuts_timeseries', '-png', '-r150');

% for is = 1:nstat
%     for im = 1:nbgc
%         h.ln(is,im,1) = plot(h.ax(is,1), Mplt(is).(bgcname{im}), 'NO3b', 'color', colmodel(im,:));
%         h.ln(is,im,2) = plot(h.ax(is,2), Mplt(is).(bgcname{im}), 'NH4b', 'color', colmodel(im,:));
%     end
% end




%% Nbudget graphs

nbudfile = fullfile(simfolder, 'analysis/nbudgets_all_SEBSmiddle.mat');
Nbud = load(nbudfile);

%% ... Colors

types = ["N","P","Z","D","B","X"];

ncmap = table(ones(length(types),3), ...
    'variablenames', {'rgb'}, 'rownames', types');

ncmap.rgb = colors{["red", "green", "blue", "brown", "orange", "gray"],"rgb"};

edgegrp = {...
    ["nit" "Nit"]
     "dnit"
    ["npp" "gpp" "Gpp"]
    ["gra" "Gra"]
    ["mor" "Mor" "vir" "exu"]
     "agg"
    ["ege" "Ege" "Exc" "Res"]
    ["rem" "Rem"]
     "sed"
    ["Frz" "Twi" "hadv", "vadv", "hdif", "vdif"]
    ["Ver", "snk"]
    };
nper = cellfun(@length, edgegrp);
egrpnum = arrayfun(@(x,y) ones(y,1)*x, 1:length(edgegrp), nper', 'uni', 0);
egrpnum = cat(1, egrpnum{:});
edgegrp = cat(2, edgegrp{:});

% ecol = rgb(["orange"; "red orange"; "green"; "blue"; "black"; "olive green"; "gold"; "burnt orange"; "yellow orange"; "light gray"; "gray"]);
ecol = colors{["red","red","green","blue","brown","brown","yellow","red","red","gray","gray"],"rgb"};


ecmap = table(ecol(egrpnum,:), 'variablenames', {'rgb'}, 'rownames', edgegrp');


%% ... Graph plot setup: colors and positions

% Location by type

ngraph = length(Nbud.G);

tmp = array2table(zeros(length(types), ngraph));
tmp.Properties.RowNames = types;

for ii = 1:ngraph
    [g, gtype] = findgroups(Nbud.G{ii}.Nodes.type);
    tmp{gtype,ii} = splitapply(@length, g, g);
end

nmax = max(tmp{:,:},[],2);

th = linspace(0, 2*pi, sum(nmax)+1);           
idx = [1; cumsum(nmax)+1];

for ii = 1:ngraph
    N = Nbud.G{ii}.Nodes;
    N.theta = nan(height(N),1);
    for it = 1:length(types)
        N.theta(strcmp(N.type, types{it})) = th(idx(it)+(1:tmp{it,ii})-1);
    end
    
    Nbud.G{ii}.Nodes.theta = N.theta;
    
end

%% ... gplotpoly version, total over one year

[~,graphplt] = ismember(["../bgcmip_loop_nbudget_banas/Out", ...
                         "../bgcmip_loop_nbudget_bestnpzbury/Out", ...
                         "../bgcmip_loop_nbudget_bestnpz/Out", ...
                         "../bgcmip_loop_nbudget_bestnpzburynoinfauna/Out", ...
                         "../bgcmip_loop_nbudget_cobalt/Out"], Nbud.outpath);
nplt = length(graphplt);
lbl = {'Banas', 'BEST_NPZ +burial', 'BEST_NPZ', 'BEST_NPZ +burial -infauna', 'COBALT'};

h = plotgrid('size', [2 3], 'sp', 0, 'mar', 0, 'sv', 0.03);
setpos(h.fig, '# # 18in 12in');

for ii = 1:nplt

    is = graphplt(ii);

    isin = year(Nbud.tg{is}) == 1993;
    Gplt = Nbud.G{is};

    % Node values: mean yearly (gC m^-2)

    ndval = max(cellfun(@(x) mean(x(isin)), Gplt.Nodes.val), 0);

    % Edge values: sum over year (gC m^-2)

    New = Gplt.Edges;
    New.EndNodes = New.EndNodes(:,[2 1]);
    New.val = cellfun(@(x) -x, New.val, 'uni', 0);
    Gplt = addedge(Gplt, New);

    edval = max(cellfun(@(x) sum(x(isin)),  Gplt.Edges.val), 0);    
    is0 = edval == 0;
    Gplt = rmedge(Gplt, find(is0));
    edval = edval(~is0);

    % ndval = cellfun(@(x) max(max(x),0), G.(sims{is}).Nodes.val);
    % edval = cellfun(@(x) max(max(x),0), G.(sims{is}).Edges.val);
    
    % First, need to decompose graphs multiple adjacency matrices to
    % account to repeated edges between nodes
    
    [~,loc] = ismember(Gplt.Edges.EndNodes, Gplt.Nodes.Name);
    [~,cloc] = ismember(Gplt.Edges.type, ecmap.Properties.RowNames);
    [~,ncloc] = ismember(Gplt.Nodes.type, ncmap.Properties.RowNames);
    
    [~,~,unq] = unique(loc, 'rows');
    cnt = diag(cumsum(unq == unq', 1));
    
    na = max(cnt);

    nnode = numnodes(Gplt);
    
    adj = zeros(nnode, nnode, na);  
    cval = zeros(nnode, nnode, na);
   
    idx = sub2ind([nnode nnode na], loc(:,1), loc(:,2), cnt);
    
    adj(idx) = sqrt(edval);
    cval(idx) = cloc;
        
    % Scale data so node area is proportional to area and and flux/day
    % matches diameter (so a flux representing x mass/day has the same 
    % width as the diameter of a node representing x mass. 
    
    axes(h.ax(ii));
    
    for ia = 1:na
        ctmp = cval(:,:,ia);
        ctmp = ctmp(adj(:,:,ia)>0);
        
        rpos = ones(nnode,1);
        rpos(ismember(Gplt.Nodes.type, ["X"])) = 1.5;
        rpos(contains(Gplt.Nodes.Name, "Ben")) = 1.5;
        xy = [rpos.*cos(Gplt.Nodes.theta), ...
              rpos.*sin(Gplt.Nodes.theta)];

        h.gp{ii}(ia) = gplotpoly(adj(:,:,ia), ...
            xy, ...
            sqrt(ndval), ...
            'edgecurve', 0.08*ia, ...
            'wlim', [0 0.25], ...
            'elim', [0 35], ...
            'rlim', [0.03 0.35], ...
            'nlim', [0 35], ...
            'ecdata', ctmp, ...
            'ncdata', ncloc+height(ecmap));
        if ia > 1
            set(h.gp{ii}(ia).node, 'visible', 'off');
        end
        
    end 
end


set(h.ax, 'dataaspectratio', [1 1 1], 'xlim', [-1 1]*1.5, 'ylim', [-1.75 1.35], ...
    'clim', [0.5 height(ecmap)+height(ncmap)+0.5], 'colormap', [ecmap.rgb; ncmap.rgb], ...
    'clipping', 'off');

h.gp = cat(2, h.gp{:});
set(cat(1, h.gp.arrow), 'edgecolor', 'none', 'facealpha', 0.5);
set(cat(1, h.gp.node), 'edgecolor', 'none', 'facealpha', 0.8);


h.lbl = labelaxes(h.ax(1:nplt), lbl, 'northwest');
set(h.lbl, 'interpreter', 'none', 'fontsize', 20, 'fontweight', 'b');
% set(h.lbl(1), 'color', colors{'green','rgb'});
% set(h.lbl(end), 'color', colors{'red','rgb'});
% multitextloc(h.ax(:,1), {Budget(bidx).name}, 'west');
set(h.ax, 'visible', 'off');
set(h.fig, 'color', 'w');
% set(h.fig, 'color', 'none');

%% ... EPOC slide versions

pltopt{1} = ["../bgcmip_loop_nbudget_banas/Out", ...
          "../bgcmip_loop_nbudget_bestnpz/Out", ...
          "../bgcmip_loop_nbudget_cobalt/Out"];
pltopt{2} = ["../bgcmip_loop_nbudget_bestnpz/Out", ...
             "../bgcmip_loop_nbudget_bestnpzbury/Out", ...
             "../bgcmip_loop_nbudget_bestnpzburynoinfauna/Out"];

pltlbl{1} = ["Banas", "BEST_NPZ", 'COBALT'];
pltlbl{2} = ["BEST_NPZ", "BEST_NPZ +burial", "BEST_NPZ +burial -infauna"]; 

iplt = 1;
popout = true;

[~,graphplt] = ismember(pltopt{iplt}, Nbud.outpath);
nplt = length(graphplt);
% lbl = {'Banas', 'BEST_NPZ +burial', 'BEST_NPZ', 'BEST_NPZ +burial -infauna', 'COBALT'};

h = plotgrid('size', [1 3], 'sp', 0, 'mar', 0);
h.fig.Position(3:4) = [1500 650];

for ii = 1:nplt

    is = graphplt(ii);

    isin = year(Nbud.tg{is}) == 1993;
    Gplt = Nbud.G{is};

    % Node values: mean yearly (gC m^-2)

    ndval = max(cellfun(@(x) mean(x(isin)), Gplt.Nodes.val), 0);

    % Edge values: sum over year (gC m^-2)

    New = Gplt.Edges;
    New.EndNodes = New.EndNodes(:,[2 1]);
    New.val = cellfun(@(x) -x, New.val, 'uni', 0);
    Gplt = addedge(Gplt, New);

    edval = max(cellfun(@(x) sum(x(isin)),  Gplt.Edges.val), 0);    
    is0 = edval == 0;
    Gplt = rmedge(Gplt, find(is0));
    edval = edval(~is0);

    % ndval = cellfun(@(x) max(max(x),0), G.(sims{is}).Nodes.val);
    % edval = cellfun(@(x) max(max(x),0), G.(sims{is}).Edges.val);
    
    % First, need to decompose graphs multiple adjacency matrices to
    % account to repeated edges between nodes
    
    [~,loc] = ismember(Gplt.Edges.EndNodes, Gplt.Nodes.Name);
    [~,cloc] = ismember(Gplt.Edges.type, ecmap.Properties.RowNames);
    [~,ncloc] = ismember(Gplt.Nodes.type, ncmap.Properties.RowNames);
    
    [~,~,unq] = unique(loc, 'rows');
    cnt = diag(cumsum(unq == unq', 1));
    
    na = max(cnt);

    nnode = numnodes(Gplt);
    
    adj = zeros(nnode, nnode, na);  
    cval = zeros(nnode, nnode, na);
   
    idx = sub2ind([nnode nnode na], loc(:,1), loc(:,2), cnt);
    
    adj(idx) = sqrt(edval);
    cval(idx) = cloc;
        
    % Scale data so node area is proportional to area and and flux/day
    % matches diameter (so a flux representing x mass/day has the same 
    % width as the diameter of a node representing x mass. 
    
    axes(h.ax(ii));
    
    for ia = 1:na
        ctmp = cval(:,:,ia);
        ctmp = ctmp(adj(:,:,ia)>0);
        
        rpos = ones(nnode,1);
        if popout
            rpos(ismember(Gplt.Nodes.type, ["X"])) = 1.5;
            rpos(contains(Gplt.Nodes.Name, "Ben")) = 1.5;
        end
        xy = [rpos.*cos(Gplt.Nodes.theta), ...
              rpos.*sin(Gplt.Nodes.theta)];

        h.gp{ii}(ia) = gplotpoly(adj(:,:,ia), ...
            xy, ...
            sqrt(ndval), ...
            'edgecurve', 0.08*ia, ...
            'wlim', [0 0.25], ...
            'elim', [0 35], ...
            'rlim', [0.03 0.35], ...
            'nlim', [0 35], ...
            'ecdata', ctmp, ...
            'ncdata', ncloc+height(ecmap));
        if ia > 1
            set(h.gp{ii}(ia).node, 'visible', 'off');
        end

        isno3 = adj(:,:,ia)>0 & strcmpi(Gplt.Nodes.Name','no3') | ...
            (strcmpi(Gplt.Nodes.Name,'no3') & strcmpi(Gplt.Nodes.Name','out'));
        isno3 = isno3(adj(:,:,ia)>0);
        set(h.gp{ii}(ia).arrow(isno3), 'tag', 'no3');

        isnh4 = adj(:,:,ia)>0 & strcmpi(Gplt.Nodes.Name','nh4');
        isnh4 = isnh4(adj(:,:,ia)>0);
        set(h.gp{ii}(ia).arrow(isnh4), 'tag', 'nh4');

        
    end 
end


set(h.ax, 'dataaspectratio', [1 1 1], 'xlim', [-1.1 1.5], 'ylim', [-1.75 1.35], ...
    'clim', [0.5 height(ecmap)+height(ncmap)+0.5], 'colormap', [ecmap.rgb; ncmap.rgb], ...
    'clipping', 'off');

h.gp = cat(2, h.gp{:});
set(cat(1, h.gp.arrow), 'edgecolor', 'none', 'facealpha', 0.5);
set(cat(1, h.gp.node), 'edgecolor', 'none', 'facealpha', 0.8);


h.lbl = labelaxes(h.ax(1:nplt), pltlbl{iplt}, 'northwest');
set(h.lbl, 'interpreter', 'none', 'fontsize', 20, 'fontweight', 'b');

set(h.ax, 'visible', 'off');
set(h.fig, 'color', 'w');


if iplt == 1 && ~popout
    export_fig(h.fig, '~/Documents/Conferences/202409_EPOC/EPOC_Kearney_slides/nbudget1_nopop', '-r150', '-nocrop');
    set(cat(1,h.gp.arrow), 'visible', 'off');
    export_fig(h.fig, '~/Documents/Conferences/202409_EPOC/EPOC_Kearney_slides/nbudget1_nopop_nodes', '-r150', '-nocrop');
end
if iplt == 1 && popout

    export_fig(h.fig, '~/Documents/Conferences/202409_EPOC/EPOC_Kearney_slides/nbudget1', '-r150', '-nocrop');

    set(cat(1, h.gp.arrow), 'facealpha', 0.1);
    set(findall(cat(1, h.gp.arrow), 'tag', 'no3'), 'facealpha', 0.5);
    set(findall(cat(1, h.gp.arrow), 'tag', 'nh4'), 'facealpha', 0.5);

    export_fig(h.fig, '~/Documents/Conferences/202409_EPOC/EPOC_Kearney_slides/nbudget1_nuts_highlight', '-r150', '-nocrop');

end
if iplt == 2
    export_fig(h.fig, '~/Documents/Conferences/202409_EPOC/EPOC_Kearney_slides/nbudget2', '-r150', '-nocrop');
end





%% Total fluxes between pelagic and benthic/out

flxname = @(nd,tp) sprintf('%s_%s_%s', tp, nd{1}, nd{2});

% Banas

[~,is] = ismember("../bgcmip_loop_nbudget_banas/Out", Nbud.outpath);
% tmask = year(Nbud.tg{is}) == 1993;
tmask = year(Nbud.tg{is}) >= 1990;
Etmp = Nbud.G{is}.Edges;

isflux = Etmp.type=="snk";

Flx.banas = array2timetable(cat(2, Etmp.val{isflux}), ...
    'rowtimes', Nbud.tg{is}, ...
    'variablenames', rowfun(flxname, Etmp(isflux, {'EndNodes', 'type'}), 'outputformat', 'cell'));
Flx.banas = Flx.banas(tmask,:);

% BEST_NPZ

[~,is] = ismember("../bgcmip_loop_nbudget_bestnpz/Out", Nbud.outpath);
% tmask = year(Nbud.tg{is}) == 1993;
tmask = year(Nbud.tg{is}) >= 1990;
Etmp = Nbud.G{is}.Edges;

isflux = (Etmp.type == "Ver" & ismember(Etmp.EndNodes(:,2), {'Out', 'DetBen'})) | ...
    (any(ismember(Etmp.EndNodes, {'Ben','DetBen'}),2) & ~all(ismember(Etmp.EndNodes, {'Ben','DetBen'}),2));

Flx.bestnpz = array2timetable(cat(2, Etmp.val{isflux}), ...
    'rowtimes', Nbud.tg{is}, ...
    'variablenames', rowfun(flxname, Etmp(isflux, {'EndNodes', 'type'}), 'outputformat', 'cell'));
Flx.bestnpz = Flx.bestnpz(tmask,:);


% COBALT

[~,is] = ismember("../bgcmip_loop_nbudget_cobalt/Out", Nbud.outpath);
% tmask = year(Nbud.tg{is}) == 1993;
tmask = year(Nbud.tg{is}) >= 1990;
Etmp = Nbud.G{is}.Edges;

isflux = any(ismember(Etmp.EndNodes, {'out'}),2) & ~endsWith(Etmp.type, "adv") & ~endsWith(Etmp.type, "dif");

Flx.cobalt = array2timetable(cat(2, Etmp.val{isflux}), ...
    'rowtimes', Nbud.tg{is}, ...
    'variablenames', rowfun(flxname, Etmp(isflux, {'EndNodes', 'type'}), 'outputformat', 'cell'));
Flx.cobalt = Flx.cobalt(tmask,:);

% BEST_NPZ minimal-benthos variant

[~,is] = ismember("../bgcmip_loop_nbudget_bestnpzburynoinfauna/Out", Nbud.outpath);
tmask = year(Nbud.tg{is}) >= 1990;
Etmp = Nbud.G{is}.Edges;

isflux = (Etmp.type == "Ver" & ismember(Etmp.EndNodes(:,2), {'Out', 'DetBen'})) | ...
    (any(ismember(Etmp.EndNodes, {'Ben','DetBen'}),2) & ~all(ismember(Etmp.EndNodes, {'Ben','DetBen'}),2));

Flx.bestnpz2 = array2timetable(cat(2, Etmp.val{isflux}), ...
    'rowtimes', Nbud.tg{is}, ...
    'variablenames', rowfun(flxname, Etmp(isflux, {'EndNodes', 'type'}), 'outputformat', 'cell'));
Flx.bestnpz2 = Flx.bestnpz2(tmask,:);

% Means

FlxAvg = structfun(@(x) varfun(@mean,x), Flx, 'uni', 0);

% Stacked bar with labels

tmpname = FlxAvg.bestnpz.Properties.VariableNames;
flipsign = endsWith(tmpname, "Ben") | ...
           endsWith(tmpname, "Out");
FlxAvg.bestnpz{1,flipsign} = FlxAvg.bestnpz{1,flipsign}*-1;
FlxAvg.bestnpz2{1,flipsign} = FlxAvg.bestnpz2{1,flipsign}*-1;

tmpnameparts = split(string(tmpname'), "_");
[tmpnamesrt, isrt] = sortrows(table(sign(FlxAvg.bestnpz{1,:})'==1, tmpnameparts(:,2), tmpnameparts(:,4), tmpnameparts(:,3)));
% return
% 
% [tmpnamesrt,isrt] = sortrows(table(sign(FlxAvg.bestnpz{1,:})', tmpname'));
FlxAvg.bestnpz = FlxAvg.bestnpz(:,isrt);
FlxAvg.bestnpz2 = FlxAvg.bestnpz2(:,isrt);

FlxAvg.cobalt.mean_snk_ndet_out = FlxAvg.cobalt.mean_snk_ndet_out*-1;


%% ... plot bar graph (appended below benthic processes schematic)

dx = 0.2;
xax = 0.25;

barlabels = [...
    "mean_snk_out_det_small"        "Sinking, small detritus"
    "mean_snk_out_det_large"        "Sinking, large detritus"
    "mean_Gra_DetF_Ben"             "Infaunal grazing, DetF"
    "mean_Gra_Det_Ben"              "Infaunal grazing, Det"
    "mean_Gra_PhL_Ben"              "Infaunal grazing, PhL"
    "mean_Gra_PhS_Ben"              "Infaunal grazing, PhS"
    "mean_Ver_DetF_DetBen"          "Sinking to DetBen, DetF"
    "mean_Ver_DetF_Out"             "Burial and denitrification, DetF"
    "mean_Ver_Det_DetBen"           "Sinking to DetBen, Det"
    "mean_Ver_Det_Out"              "Burial and denitrification, Det"
    "mean_Ver_NCaO_DetBen"          "Sinking to DetBen, NCaO"
    "mean_Ver_PhL_DetBen"           "Sinking to DetBen, PhL"
    "mean_Ver_PhL_Out"              "Burial and denitrification, PhL"
    "mean_Ver_PhS_DetBen"           "Sinking to DetBen, PhS"
    "mean_Ver_PhS_Out"              "Burial and denitrification, PhS"
    "mean_Ver_NCaS_DetBen"          "Sinking to DetBen, NCaS"
    "mean_Exc_Ben_NH4"              "Infauna excretion"
    "mean_Rem_DetBen_NH4"           "Benthic remineralization"
    "mean_Res_Ben_NH4"              "Infauna respiration"
    "mean_snk_ndet_out"             "Sinking, N detritus"
    "mean_sed_out_nh4"              "Sedimentary remineralization"
    "mean_sed_out_no3"              "Sedimentary denitrification"
    ];

% Plot

h = plotgrid('size', [1 1], 'mr', 0, 'mt', 0.05, 'mb', 0.05, 'ml', xax/12.33);
setpos(h.fig, '5in 5in 12.33in 2in');

n = size(FlxAvg.banas,2)+1;
x = 1.5-xax+((1:n)-mean([1 n]))*dx;
h.b(1,1) = bar(x(1:end-1), FlxAvg.banas{1,:});
hold on;
h.b(1,2) = bar(x(end-1:end), [NaN sum(FlxAvg.banas{1,:})]);
% [~,loc] = ismember(FlxAvg.banas.Properties.VariableNames, barlabels(:,1));
% h.txt{1} = text(x, FlxAvg.banas{1,:}, barlabels(loc,2), 'rotation', 90, 'color', colmodel(1,:), 'fontsize', 8);

n = size(FlxAvg.bestnpz,2)+1;
x = 5.5-xax+((1:n)-mean([1 n]))*dx;
h.b(2,1) = bar(x(1:end-1), FlxAvg.bestnpz{1,:});
h.b(2,2) = bar(x(end-1:end), [NaN sum(FlxAvg.bestnpz{1,:})]);

n = size(FlxAvg.bestnpz,2)+1;
x = 5.5-xax+((1:n)-mean([1 n]))*dx;
h.b(4,1) = bar(x(1:end-1), FlxAvg.bestnpz2{1,:});
h.b(4,2) = bar(x(end-1:end), [NaN sum(FlxAvg.bestnpz2{1,:})]);


n = size(FlxAvg.cobalt,2)+1;
x = 10.5-xax+((1:n)-mean([1 n]))*dx;
h.b(3,1) = bar(x(1:end-1), FlxAvg.cobalt{1,:});
% h.b(3,2) = bar(x(end), sum(FlxAvg.cobalt{1,:}));
h.b(3,2) = bar(x(end-1:end), [NaN sum(FlxAvg.cobalt{1,:})]);

set(h.ax, 'xlim', [0 12.33], 'box', 'off', 'xcolor', 'none', 'fontsize', 8, 'ylim', [-2 1.5]);
labelaxes(h.ax, "Boundary flux (mmol N m^{-2} d^{-1})", 'northwest', 'fontsize', 10);

set(h.b(1:3,1), {'facecolor'}, num2cell(colmodel,2),{'edgecolor'}, num2cell(colmodel,2), 'linewidth', 1.0);
set(h.b(1:3,2), {'edgecolor'}, num2cell(colmodel,2), 'facecolor', 'none', 'linewidth', 1.5);
set(h.b(4,:), 'facecolor', 'none', 'edgecolor', colors{'yellow','rgb'}, 'linestyle', '-', 'linewidth', 2)
arrayfun(@(x) set(get(x, 'BaseLine'), 'color', rgb('light gray')), h.b);
uistack(h.b(1:3,2), 'top');
set(h.fig, 'color', 'w');

annotation('rectangle', [0 0 1 1], 'edgecolor', 'w');

return
exportgraphics(h.fig, 'bottom_flux_bars.pdf');

% TDOD Add labels on fluxes?  Words?  Schematics (type, x -> y w/
% pictographs for x,y?)

%% Is each node a source or sink?

[~,is] = ismember("../bgcmip_loop_nbudget_banas/Out", Nbud.outpath);
% tmask = year(Nbud.tg{is}) == 1993;
tmask = year(Nbud.tg{is}) >= 1990;
Etmp = Nbud.G{is}.Edges;
Ntmp = Nbud.G{is}.Nodes;

val = cellfun(@mean, Etmp.val);
for ii = 1:height(Ntmp)
    issrc = strcmp(Etmp.EndNodes(:,1), Ntmp.Name{ii});
    issnk = strcmp(Etmp.EndNodes(:,2), Ntmp.Name{ii});
    Ntmp.netval(ii) = sum(val(issrc))*-1 + sum(val(issnk));
end


[~,is] = ismember("../bgcmip_nbudget_bestnpz/Out", Nbud.outpath);
% tmask = year(Nbud.tg{is}) == 1993;
tmask = year(Nbud.tg{is}) == 1990;
Etmp = Nbud.G{is}.Edges;
Ntmp = Nbud.G{is}.Nodes;

val = cellfun(@mean, Etmp.val);
for ii = 1:height(Ntmp)
    issrc = strcmp(Etmp.EndNodes(:,1), Ntmp.Name{ii});
    issnk = strcmp(Etmp.EndNodes(:,2), Ntmp.Name{ii});
    Ntmp.netval(ii) = sum(val(issrc))*-1 + sum(val(issnk));
end

%% Total SEBS budget

nbudfile = fullfile(simfolder, 'analysis/nbudgets_all_SEBS.mat');
Nbud = load(nbudfile);

maskarea = sum(Grd.area_feast(Budget(5).mask)) * 1e6; % m^2


pltopt{1} = ["../bgcmip_loop_nbudget_banas/Out", ...
          "../bgcmip_loop_nbudget_bestnpz/Out", ...
          "../bgcmip_loop_nbudget_cobalt/Out"];
pltopt{2} = ["../bgcmip_loop_nbudget_bestnpz/Out", ...
             "../bgcmip_loop_nbudget_bestnpzbury/Out", ...
             "../bgcmip_loop_nbudget_bestnpzburynoinfauna/Out"];

pltlbl{1} = ["Banas", "BEST_NPZ", 'COBALT'];
pltlbl{2} = ["BEST_NPZ", "BEST_NPZ +burial", "BEST_NPZ +burial -infauna"]; 

iplt = 2;
popout = true;

[~,graphplt] = ismember(pltopt{iplt}, Nbud.outpath);
nplt = length(graphplt);

% Node color brightness by number-of-type

% ngraph = length(Nbud.G);
% 
% tmp = array2table(zeros(length(types), ngraph));
% tmp.Properties.RowNames = types;
% 
% for ii = 1:ngraph
%     [g, gtype] = findgroups(Nbud.G{ii}.Nodes.type);
%     tmp{gtype,ii} = splitapply(@length, g, g);
% end
% 
% nmax = max(tmp{:,:},[],2);
% 
% th = linspace(0, 2*pi, sum(nmax)+1);           
% idx = [1; cumsum(nmax)+1];

% cfac = [0:0.15:1.0; -0.15:-0.15:-1.0 NaN];
% 
% for ii = 1:ngraph
%     N = Nbud.G{ii}.Nodes;
%     N.colfac = nan(height(N),1);
%     for it = 1:length(types)
%         N.colfac(strcmp(N.type, types{it})) = cfac(1:tmp{it,ii});
%     end
% 
%     Nbud.G{ii}.Nodes.colfac = N.colfac;
% 
% end

% Classify state variables and assign type brightness based on that

vtype = {...
    "N, NH4"         ["nh4" "NH4"]                                         ncmap{'N','rgb'}
    "N, NO3"         ["no3", "NO3"]                              adjustcol(ncmap{'N','rgb'}, -0.2)
    "N, other"       ["IceNH4", "IceNO3", "n2"]                  adjustcol(ncmap{'N','rgb'},  0.4)
    "P"              "phyto"                                               ncmap{'P','rgb'}
    "P, small"       ["PhS" "nsm"]                               adjustcol(ncmap{'P','rgb'}, -0.2)
    "P, large"       ["PhL" "nlg"]                               adjustcol(ncmap{'P','rgb'},  0.2)
    "P, other"       ["IcePhL", "ndi"]                           adjustcol(ncmap{'P','rgb'},  0.6)
    "Z, micro"       ["microzoo" "MZL" "nsmz"]                             ncmap{'Z','rgb'}
    "Z, meso small"  ["Cop" "nmdz"]                              adjustcol(ncmap{'Z','rgb'}, -0.2)
    "Z, meso large"  ["NCaS" "NCaO" "EupS" "EupO" "nlgz"]        adjustcol(ncmap{'Z','rgb'}, 0.2)
    "Z, other"       ["Jel" "hip"]                               adjustcol(ncmap{'Z','rgb'}, -0.4)
    "Z, benthic"      "Ben"                                      adjustcol(ncmap{'Z','rgb'},  0.4)
    "D, small/dissolved" ["det_small", "Det" "ldon" "sldon" "srdon"]       ncmap{'D','rgb'}
    "D, particulate" ["det_large", "DetF", "ndet"]               adjustcol(ncmap{'D','rgb'}, -0.4)
    "D, benthic"     "DetBen"                                    adjustcol(ncmap{'D','rgb'},  0.4)
    "B"              "nbact"                                               ncmap{'B','rgb'}
    "X"              "out"                                                 ncmap{'X','rgb'}
    };
vtype = cell2table(vtype, 'variablenames', ["type", "vars", "rgb"]);

%% ... Mean total amount in each box

% Read total in each state variable

ngpermol = 14.0067;

Ntmp = cell(nplt,1);
for ii = 1:nplt

    Ntmp{ii} = Nbud.G{graphplt(ii)}.Nodes;
    Ntmp{ii}.meanval = cellfun(@mean, Ntmp{ii}.val) * maskarea * ngpermol * 1e-3 * 1e-12; % Tg N


    tf = cellfun(@(x) ismember(Ntmp{ii}.Name, x), vtype.vars, 'uni', 0);
    [rr,cc] = find(cat(2,tf{:}));
    Ntmp{ii}.subtype(rr) = vtype.type(cc);
    Ntmp{ii}.color(rr,:) = vtype.rgb(cc,:);

end

% Relative size of each total

ntot = cellfun(@(x) sum(x.meanval), Ntmp);
nrel = sqrt(ntot./max(ntot));

% Plot

h = plotgrid('size', [1 1], 'mar', 0.05, 'mb', 0.4);
hold(h.ax, 'on');
h.fig.Position(3:4) = [800 600];

leganchor = {{'sw','nw'}, {'s','n'}, {'se','ne'}};

for ii = 1:nplt
    [g, gtype] = findgroups(Ntmp{ii}.type);
    meanvaltot = splitapply(@sum, Ntmp{ii}.meanval, g);

    % Plot treemaps

    r = treemap(meanvaltot, nrel(ii), nrel(ii));
    rsub = nan(4,height(Ntmp{ii}));

    for itp = 1:max(g)
        rNew = treemap(Ntmp{ii}.meanval(g==itp),r(3,itp),r(4,itp));
        rNew(1,:) = rNew(1,:) + r(1,itp);
        rNew(2,:) = rNew(2,:) + r(2,itp);
        rsub(:,g==itp) = rNew;
    end

   
    % coltmp = ncmap{Ntmp{ii}.type,'rgb'};
    % for ic = 1:size(coltmp,1)
    %     coltmp(ic,:) = adjustcol(coltmp(ic,:), Ntmp{ii}.colfac(ic));
    % end

    rplt = rsub;
    rplt(1,:) = rplt(1,:) + sum(nrel(1:ii-1))+(ii-1)*0.1;

    xplt = [rplt(1,:); rplt(1,:)+rplt(3,:)];
    yplt = [rplt(2,:); rplt(2,:)+rplt(4,:)];

    h.sub(ii) = patch(xplt([1 1 2 2 1],:), yplt([1 2 2 1 1],:), 'r');
    set(h.sub(ii), 'cdata', permute(Ntmp{ii}.color, [1 3 2]), 'facecolor', 'flat', 'edgecolor', 'none');

    rplt = r;
    rplt(1,:) = rplt(1,:) + sum(nrel(1:ii-1))+(ii-1)*0.1;

    xplt = [rplt(1,:); rplt(1,:)+rplt(3,:)];
    yplt = [rplt(2,:); rplt(2,:)+rplt(4,:)];

    h.type(ii) = patch(xplt([1 1 2 2 1],:), yplt([1 2 2 1 1],:), 'w');
    set(h.type(ii), 'facecolor', 'none', 'edgecolor', 'w');

    % Invisible patches for labels

    for ic = 1:height(Ntmp{ii})-1
        h.ptmp(ic) = patch([0 1 1 0], [0 0 1 1]+100, Ntmp{ii}.color(ic,:));
    end
    set(h.ptmp, 'visible', 'off');
    % legend(h.ptmp, Ntmp{ii}.Name);

    lblstr = compose("%9s: ", string(Ntmp{ii}.Name)) + ...
             compose("%5.3f Tg N ", Ntmp{ii}.meanval) + ...
             compose("(%5.2f%%)", 100*Ntmp{ii}.meanval./sum(Ntmp{ii}.meanval));

    h.leg(ii) = legendflex(h.ptmp, cellstr(lblstr(1:end-1)), 'ref', h.ax, 'anchor', leganchor{ii}, 'buffer', [0 0], ...
        'interpreter', 'none', 'fontsize', 8, 'box', 'off', 'xscale', 0.5, 'fontname', 'menlo');
end
set(h.fig, 'color', 'w');
axis tight equal;
set(h.ax, 'visible', 'off');


if iplt == 1

    title(h.leg(1), sprintf('Banas: %.2f Tg N', sum(Ntmp{1}.meanval)));
    title(h.leg(2), sprintf('BEST\\_NPZ: %.2f Tg N', sum(Ntmp{2}.meanval)));
    title(h.leg(3), sprintf('COBALT: %.2f Tg N', sum(Ntmp{3}.meanval)));
    
    export_fig('nbudgets', h.fig, '-png', '-r150');
elseif iplt == 2

    title(h.leg(1), sprintf('default: %.2f Tg N', sum(Ntmp{1}.meanval)));
    title(h.leg(2), sprintf('all sinking buried: %.2f Tg N', sum(Ntmp{2}.meanval)));
    title(h.leg(3), sprintf('no infauna grazing: %.2f Tg N', sum(Ntmp{3}.meanval)));
    
    export_fig('nbudgets_bestnpz', h.fig, '-png', '-r150');
end

%% ... simplified graphs with treemap breakdown

% Run Total SEBS budget block first

% Read total in each state variable

ngpermol = 14.0067;

r1 = 5.0; % radius of PZD nodes
r2 = 3.0; % distance to X nodes

% Node and edge position: N in center, PZD around that, X entending out

ndtype = ["N", "P", "Z", "D", "NX", "PX", "ZX", "DX"];
nnodeplt = length(ndtype);

xnode = [0 r1 r1*cosd(120) r1*cosd(240) r2*cosd(90)  r1+r2 (r1+r2)*cosd(120) (r1+r2)*cosd(240)];
ynode = [0  0 r1*sind(120) r1*sind(240) r2*sind(90)      0 (r1+r2)*sind(120) (r1+r2)*sind(240)];

ncoord = table(ndtype', xnode', ynode', 'variablenames', {'type', 'x', 'y'}, 'rownames', ndtype);

ecoord = {...
  "N" "P"  r1*2/3           0 
  "N" "Z"  r1*2/3*cosd(120) r1*2/3*sind(120)
  "N" "D"  r1*2/3*cosd(240) r1*2/3*sind(240)
  "P" "Z"  r1.*cosd( 60) r1.*sind( 60)
  "P" "D"  r1.*cosd(300) r1.*sind(300)
  "Z" "D"  r1.*cosd(180) r1.*sind(180)
  "N" "NX" (r2*2/3)*cosd(90)  (r2*2/3)*sind(90)
  "P" "PX" (r1+r2*2/3)*cosd(  0)  (r1+r2*2/3)*sind(  0)
  "Z" "ZX" (r1+r2*2/3)*cosd(120)  (r1+r2*2/3)*sind(120)
  "D" "DX" (r1+r2*2/3)*cosd(240)  (r1+r2*2/3)*sind(240)
  "N" "N"  ncoord{"N","x"} ncoord{"N","y"}
  "P" "P"  ncoord{"P","x"} ncoord{"P","y"}
  "Z" "Z"  ncoord{"Z","x"} ncoord{"Z","y"}
  "D" "D"  ncoord{"D","x"} ncoord{"D","y"}
};
ecoord = cell2table(ecoord, 'variablenames', {'src','snk','x','y'});
nedgeplt = height(ecoord);

% Shared flux path lines

circs = {...
    "P" "Z" r1.*cosd(linspace(  0,120,30)) r1.*sind(linspace(  0,120,30))
    "Z" "D" r1.*cosd(linspace(120,240,30)) r1.*sind(linspace(120,240,30))
    "P" "D" fliplr(r1.*cosd(linspace(240,360,30))) fliplr(r1.*sind(linspace(240,360,30)))
};

[~,loc] = ismember(string(circs(:,1:2)), ecoord{:,1:2}, 'rows');
ecoord.xseg(loc) = circs(:,3);
ecoord.yseg(loc) = circs(:,4);


% th = linspace(0,2*pi,21);
% th2 = linspace(0,2*pi,100);
% xcirc = r1*cos(th2);
% ycirc = r1*sin(th2);

straight = [...
    "N" "P"  
    "N" "Z"  
    "N" "D"  
    "N" "NX" 
    "P" "PX" 
    "Z" "ZX" 
    "D" "DX"];
[~,segs] = ismember(straight', ncoord.type);

[~,loc] = ismember(straight, ecoord{:,1:2}, 'rows');
ecoord.xseg(loc) = num2cell(ncoord.x(segs)', 2);
ecoord.yseg(loc) = num2cell(ncoord.y(segs)', 2);

% Figure setup

h = plotgrid('size', [1 length(graphplt)], 'mar', 0.02, 'mb', 0.5, 'sp', 0.03);
h.fig.Position(3:4) = [1000 800];

leganchor = {{'sw','nw'}, {'s','n'}, {'se','ne'}};

for ig = 1:length(graphplt)

    % Full graph

    Gplt = Nbud.G{graphplt(ig)};

    Gplt.Nodes.meanval = cellfun(@mean, Gplt.Nodes.val) * maskarea * ngpermol * 1e-3 * 1e-12; % Tg N
    Gplt.Edges.meanval = cellfun(@mean, Gplt.Edges.val) * maskarea * ngpermol * 1e-3 * 1e-12 * 365; % Tg N/yr

    % Simplify nodes by type

    ndtypetmp = ["N","P","Z","D","X"];
   
    [~,nodeloc] = ismember(Gplt.Nodes.type, ndtypetmp);
    nodeloc(Gplt.Nodes.type=="B") = 3; % Treat bacteria as consumer for plotting

    nodesum = splitapply(@sum, Gplt.Nodes.meanval, nodeloc);

    % Calculte treemap subdivisions

    wnd = sqrt(nodesum(1:4));

    r = [ncoord.x(1:4)-wnd/2 ncoord.y(1:4)-wnd/2 wnd wnd]';
    Gplt.Nodes.rsub = nan(height(Gplt.Nodes),4);

    for itp = 1:4
        rNew = treemap(Gplt.Nodes.meanval(nodeloc==itp),r(3,itp),r(4,itp));
        rNew(1,:) = rNew(1,:) + r(1,itp);
        rNew(2,:) = rNew(2,:) + r(2,itp);
        Gplt.Nodes.rsub(nodeloc==itp,:) = rNew';
    end

    % Assign colors

    tf = cellfun(@(x) ismember(Gplt.Nodes.Name, x), vtype.vars, 'uni', 0);
    [rr,cc] = find(cat(2,tf{:}));
    Gplt.Nodes.subtype(rr) = vtype.type(cc);
    Gplt.Nodes.color(rr,:) = vtype.rgb(cc,:);

    % Simplify edges to net between types (but keep out-of-system separate)

    [~,enum] = ismember(Gplt.Edges.EndNodes, Gplt.Nodes.Name);
    etype = ndtypetmp(nodeloc(enum));
 
    for ie = 1:size(etype,1)
        if etype(ie,1) == "X"
            etype(ie,1) = etype(ie,2)+"X";
        end
        if etype(ie,2) == "X"
            etype(ie,2) = etype(ie,1)+"X";
        end
    end

    tf = ismember(etype, ecoord{:,["src" "snk"]}, 'rows') | ...
         ismember(etype(:,[2 1]), ecoord{:,["src" "snk"]}, 'rows');
    if ~all(tf)
        error('missing?');
    end

    edgesum = zeros(nedgeplt,1);
    for ie = 1:nedgeplt
        mask1 = ismember(etype, ecoord{ie,["src" "snk"]}, 'rows');
        mask2 = ismember(etype, ecoord{ie,["snk" "src"]}, 'rows') & ~mask1;
        edgesum(ie) = sum(Gplt.Edges.meanval(mask1)) + sum(Gplt.Edges.meanval(mask2)*-1);
    end

    % Plot shared lines

    axes(h.ax(ig));

    for is = 1:10
        h.seg(is) = patch([ecoord.xseg{is} NaN], [ecoord.yseg{is} NaN], 0);
        hold on;
        if edgesum(is) > 0
            set(h.seg(is), 'cdata', [linspace(0,1,length(ecoord.xseg{is})) NaN]);
        else
            set(h.seg(is), 'cdata', [linspace(1,0,length(ecoord.xseg{is})) NaN]);
        end
    end
    set(h.seg, 'edgecolor','interp','facecolor','none', 'linewidth', 2);
    
    % Plot node boxes

    wnd = sqrt(nodesum(1:4));

    % h.node = patch(([0 0 1 1 0].*wnd - wnd/2 + xnode(1:4)')', ... 
    %                ([0 1 1 0 0].*wnd - wnd/2 + ynode(1:4)')', 'r');
    h.out = plot(xnode(5:end), ynode(5:end), 'kx');
    set(h.out, 'markersize', 8);

    % % Add treemaps

    x1 = Gplt.Nodes.rsub(:,1);
    x2 = Gplt.Nodes.rsub(:,1)+Gplt.Nodes.rsub(:,3);
    y1 = Gplt.Nodes.rsub(:,2);
    y2 = Gplt.Nodes.rsub(:,2)+Gplt.Nodes.rsub(:,4);

    h.nodes = patch([x1 x1 x2 x2 x1]', [y1 y2 y2 y1 y1]', 'w'); %, Gplt.Nodes.color);
    set(h.nodes, 'cdata', permute(Gplt.Nodes.color, [1 3 2]), 'facecolor', 'flat', 'edgecolor', 'none');

    % Plot edge flux circles

    r = sqrt(abs(edgesum/3)./pi);
    h.edge = patch((r(1:10).*cos(th) + ecoord.x(1:10))', ...
                   (r(1:10).*sin(th) + ecoord.y(1:10))', 'g');

    h.elbl = text(ecoord.x(1:10), ecoord.y(1:10), compose("%.2f", abs(edgesum(1:10))), 'horiz', 'center', 'fontsize', 8);
    set(h.edge, 'facecolor', 'w', 'edgecolor', rgb('gray'));

    % Some formatting

    axis tight equal;
    set(h.ax, 'visible', 'off'); 
    set(h.fig, 'color', 'w');

    ttlstr = sprintf("%s: %.2f Tg N", pltlbl{iplt}{ig}, sum(Gplt.Nodes.meanval));
    labelaxes(h.ax(ig), ttlstr, 'northoutside', 'fontweight', 'b', 'interpreter', 'none');

    % Invisible patches for labels

    h.ptmp = gobjects(height(Gplt.Nodes)-1,1);
    for ic = 1:height(Gplt.Nodes)-1
        h.ptmp(ic) = patch([0 1 1 0], [0 0 1 1]+100, Gplt.Nodes.color(ic,:));
    end
    set(h.ptmp, 'visible', 'off');

    lblstr = compose("%9s: ", string(Gplt.Nodes.Name)) + ...
             compose("%5.3f Tg N ", Gplt.Nodes.meanval) + ...
             compose("(%5.2f%%)", 100*Gplt.Nodes.meanval./sum(Gplt.Nodes.meanval));

    h.leg = legendflex(h.ptmp, cellstr(lblstr(1:end-1)), 'ref', h.ax(ig), 'anchor', leganchor{2}, 'buffer', [0 0], ...
        'interpreter', 'none', 'fontsize', 10, 'box', 'off', 'xscale', 0.5, 'fontname', 'menlo');


end

set(h.ax, 'colormap', cmocean('-gray'), 'clim', [0 1.2]);

if iplt == 1
    export_fig(h.fig, 'nbudgets_graph', '-png', '-r150');
else
    export_fig(h.fig, 'nbudgets_graph_bestnpz', '-png', '-r150');
end

% TODO: add arrow directions (maybe line gradient?), treemap the nodes






%% EPOC slide map of domain

xlimfull = minmax(P.xgrd, 'expand', 0.1);
ylimfull = minmax(P.ygrd, 'expand', 0.05);
yheight = 850;

% Create figure

h = plotgrid('size', [1 1], 'margin', 0.0);

h.fig.Position(3:4) = [yheight.*diff(xlimfull)./diff(ylimfull) yheight];

% setpos(h.fig, '# # 800 600');
set(h.fig, 'color', 'w');

% B10K grid

% axes(h.ax(2));

h.roms = pcolor(P.xgrd, P.ygrd, Grd.h(2:end,2:end));
shading flat;
% set(h.roms, 'facealpha', 0.6);
cmapbathy = cmocean('-gray'); %adjustcol(colors{'purple','rgb'}, linspace(1,0,50));
set(h.ax, 'colormap', cmapbathy, 'clim', [0 7000]);
hold(h.ax, 'on');

% Borders

plot(P.xbor, P.ybor, 'color', 'k', 'linewidth', 1.5);


% Colorbar

h.cb = colorbar('south');
h.cb.Position(1) = h.ax.Position(1)+0.1*h.ax.Position(3);
h.cb.Position(2) = h.ax.Position(2)+0.05*h.ax.Position(4);
h.cb.Position(3:4) = [0.3 0.05];
% setpos(h.cb, '0.1 # 0.3 #');
set(h.cb, 'ticks', [], 'edgecolor', colors{'dark','rgb'});
tprops = {'edgecolor', 'none', 'fontsize', 16, 'vert', 'middle'};
h.an(1) = annotation('textbox', h.cb.Position, 'string', '0m', 'horiz', 'left', tprops{:});
h.an(2) = annotation('textbox', h.cb.Position, 'string', '7000m', 'horiz', 'right', tprops{:}, 'color', 'w');

% Axes

set(h.ax, 'xlim', minmax(P.xgrd, 'expand', 0.1), 'ylim', minmax(P.ygrd, 'expand', 0.05), 'dataaspectratio', [1 1 1], ...
    'layer', 'top', 'xtick', [], 'ytick', [], 'visible', 'off');


% exportgraphics(h.fig, '~/Documents/Conferences/202409_EPOC/EPOC_Kearney_slides/B10K_simple_map.png', 'Resolution', 150);

%% ... shelf only, with stations used for metrics grid

% Highlight shelf box, then zoom

xboxlim = minmax(P.xstat, 'expand', 0.2);
yboxlim = minmax(P.ystat, 'expand', 0.15);

h.box(1) = plot(xboxlim([1 1 2 2 1]), yboxlim([1 2 2 1 1]), 'color', cmapbathy(1,:), 'linewidth', 2);
h.box(2) = plot(xboxlim([1 1 2 2 1]), yboxlim([1 2 2 1 1]), 'color', cmapbathy(end,:), 'linewidth', 2, 'linestyle', '--');

exportgraphics(h.fig, '~/Documents/Conferences/202409_EPOC/EPOC_Kearney_slides/B10K_simple_map_box.png', 'Resolution', 150);

% Zoomed version

h.fig.Position(3) = h.fig.Position(3) .* diff(xboxlim)./diff(xlimfull);
h.fig.Position(4) = h.fig.Position(4) .* diff(yboxlim)./diff(ylimfull);

h.con = contour(P.xgrd, P.ygrd, Grd.h(2:end,2:end), [50 100 200], 'color', cmapbathy(100,:));

isgrd = ~(Station.letter == "");

h.stat = text(P.xstat(isgrd), P.ystat(isgrd), Station.letter(isgrd), 'fontname', 'menlo', 'fontsize', 16, 'fontweight', 'b');
set(h.ax, 'xlim', minmax(P.xstat, 'expand', 0.2), 'ylim', minmax(P.ystat, 'expand'));
set(h.box, 'visible', 'off');
set(h.cb, 'visible', 'off');
set(h.an, 'visible', 'off');
pause(1);

set(h.stat(ismember(string({h.stat.String}), ["e","h","k"])), 'color', colors{'orange','rgb'});
set(h.stat(ismember(string({h.stat.String}), ["f","i","l"])), 'color', colors{'yellow','rgb'});
set(h.stat(ismember(string({h.stat.String}), ["g","j","m"])), 'color', colors{'red','rgb'});
set(h.stat(ismember(string({h.stat.String}), ["a","b","c","d"])), 'color', colors{'gray','rgb'});
set(h.stat(ismember(string({h.stat.String}), ["n","o","p","q"])), 'color', colors{'brown','rgb'});

exportgraphics(h.fig, '~/Documents/Conferences/202409_EPOC/EPOC_Kearney_slides/B10K_simple_map_zoomed.png', 'Resolution', 150);

%% Light and temp limitation functions

bappfol = '~/Documents/Repos/ESMs/bering-Apps/';

Param{1} = yaml.loadFile(fullfile(bappfol, "subApps", "BIO_BANAS", "banas_bpar.yaml"));
Param{2} = yaml.loadFile(fullfile(bappfol, "subApps", "BEST_NPZ", "bestnpz_bpar.yaml"));
Param{3} = yaml.loadFile(fullfile(bappfol, "subApps", "BIO_COBALT", "cobalt_bpar.yaml"));

Ilim.temp = [0 5 10];
Ilim.irr = linspace(0, 400, 100);

Ilim.muplt = nan(8, length(Ilim.irr), length(Ilim.temp));

for it = 1:length(Ilim.temp)

    T = Ilim.temp(it);

    % Banas
    
    mu = [...
        Param{1}.Q_P.^(T./10) .* Param{1}.mu0 .* (Param{1}.alpha_win.*Ilim.irr./sqrt(Param{1}.alpha_win.^2.*Ilim.irr.^2 + Param{1}.mu0.^2));
        Param{1}.Q_P.^(T./10) .* Param{1}.mu0 .* (Param{1}.alpha_sum.*Ilim.irr./sqrt(Param{1}.alpha_sum.^2.*Ilim.irr.^2 + Param{1}.mu0.^2));
    ];
    
    % BEST_NPZ
    
    PmaxS = log(2).*Param{2}.DiS .* 10.0.^(Param{2}.DpS.*T);
    PmaxL = log(2).*Param{2}.DiL .* 10.0.^(Param{2}.DpL.*T);
    
    watts2photons = 0.394848; % W m^-2 -> E/m^2/d
    
    mu = [...
        mu
        PmaxS .* tanh(Param{2}.alphaPhS.*Ilim.irr.*watts2photons/(PmaxS.*Param{2}.ccr))
        PmaxL .* tanh(Param{2}.alphaPhL.*Ilim.irr.*watts2photons/(PmaxL.*Param{2}.ccrPhL))
        ];
    
    % COBALT
    
    cfac = 2.77e-18;
    sperd = 86400;
    pcm = [exp(Param{3}.kappa_eppley*T).*Param{3}.P_C_max_Sm; ...
           exp(Param{3}.kappa_eppley*T).*Param{3}.P_C_max_Lg]; % s^-1
    mu = [...
        mu
        (pcm(1).*sperd./(1+0.05)).*(1-exp(-Param{3}.alpha_Sm.*Ilim.irr.*Param{3}.thetamax_Sm./pcm(1)))-Param{3}.bresp_Sm.*sperd
        (pcm(2).*sperd./(1+0.05)).*(1-exp(-Param{3}.alpha_Lg.*Ilim.irr.*Param{3}.thetamax_Lg./pcm(2)))-Param{3}.bresp_Lg.*sperd 
        (pcm(1).*sperd./(1+0.05)).*(1-exp(-Param{3}.alpha_Sm.*Ilim.irr.*Param{3}.thetamin./pcm(1)))-Param{3}.bresp_Sm.*sperd
        (pcm(2).*sperd./(1+0.05)).*(1-exp(-Param{3}.alpha_Lg.*Ilim.irr.*Param{3}.thetamin./pcm(2)))-Param{3}.bresp_Lg.*sperd 
        ];

    Ilim.muplt(:,:,it) = mu;

end


%% ... Idealized comparison of max possible primary production, M2

yr = 1993;

% Water column swrad and temp from one year

F = dir(fullfile(simfolder, 'bgcmip_loop_banas/Out', 'bgcmip_loop_banas_avg_*.nc'));
F = fullfile({F.folder}, {F.name});
F = F(1:40); % Just to speed things up...
tlim = ncdatelim(F, 'ocean_time');
isin = ~(tlim(:,2)<=datetime(yr,1,1) | tlim(:,1)>datetime(yr+1,1,1));
F = F(isin); % first year (loop)
[ixi,ieta] = ind2sub(size(Grd.h), Station.gridcell(2));
Scs = struct('xi_rho', [ixi 1 1], 'eta_rho', [ieta 1 1]);
Test = ncstruct(F, 'temp', 'zeta', 'Akt_bak', Scs);
Test = structfun(@squeeze, Test, 'uni', 0);
Test.t = ncdateread(F, 'ocean_time');
stemp = timetable(Test.t, Test.temp(end,:)');

frc = fullfile(mounteddir('klone'), 'GR011377_bgcmip', 'ROMS_Datasets', 'CFS', string(yr), "CFS-atmos-northPacific-swrad-"+yr+".nc");
Frc = ncstruct(frc, 'lat', 'lon');
[ln,lt] = ndgrid(Frc.lon, Frc.lat);
[~,imin] = pdist2(...
    [wrapTo360(ln(:)) lt(:)], ...
    [Grd.lon_rho(Station.gridcell(2)) Grd.lat_rho(Station.gridcell(2))], ...
    'euclidean', 'smallest', 1);
[iln, ilt] = ind2sub(size(ln), imin);
Swrad = ncstruct(frc, struct('lat', [ilt 1 1], 'lon', [iln 1 1]));
Swrad.t = ncdateread(frc, 'srf_time', Swrad.srf_time);

% Apply diurnal cycle a la ROMS DIURNAL_SRFLUX
% (Note: I shifted longitude so we don't get the weird discontinuity at day
% switches)

swrad = timetable(Swrad.t, squeeze(Swrad.swrad));
swrad = retime(swrad, 'hourly', 'previous');
swrad.Properties.VariableNames = {'swrad'};
yday = doy(swrad.Time);

Dangle=deg2rad(23.44.*cosd((172.0-yday)*2.0*pi/365.2425));
Hangle=(12.0-hour(swrad.Time))*pi/12.0;

LatRad=deg2rad(Grd.lat_rho(ixi,ieta));
cff1=sin(LatRad)*sin(Dangle);
cff2=cos(LatRad)*cos(Dangle);

srflx=max(0, swrad.swrad);
for ii = 1:length(srflx)
    if abs(cff1(ii))>abs(cff2(ii))
        if (cff1(ii)*cff2(ii) > 0)
            cff = cff1(ii);
            srflx(ii) = max(0, srflx(ii)/cff*(cff1(ii)+cff2(ii)*cos(Hangle(ii)-deg2rad(Grd.lon_rho(ixi,ieta)))));
        else
            srflx(ii) = 0;
        end
    else
        cff=(cff1(ii)*acos(-cff1(ii)/cff2(ii))+sqrt((cff2(ii)+cff1(ii))*(cff2(ii)-cff1(ii))))/pi;
        if cff < 10e-10
            srflx(ii) = 0;
        else
            srflx(ii) = max(0, srflx(ii)/cff*(cff1(ii)+cff2(ii)*cos(Hangle(ii)-deg2rad(Grd.lon_rho(ixi,ieta)+180))));
        end
    end
end
swrad.srflx = srflx;

stemp = retime(stemp, swrad.Time, 'nearest');

%% ... Production in the surface-most layer

% BGC parameters

sims = ["banas", "bestnpz", "cobalt"];
berapp = fullfile(mounteddir('klone'), 'GR011377_bgcmip', 'bering-Apps');
% berapp = fullfile(moxdir, 'kearney', 'bering-Apps');
pfile = fullfile(berapp, 'Apps', 'Bering_BGC_variants', "bering_bpar_"+sims+".yaml");
Param = cellfun(@(x) yaml.loadFile(x), pfile, 'uni', 0);

% Environmental variables

E0 = (1-0.58) .* swrad.srflx; % W m^-2, surface
PAR = E0;       % surface PAR, proxy for best case light levels
T = stemp.Var1; % surface temp (proxy for ML temp)

% Banas

maxkappa = Test.Akt_bak(1)/86400;  % diffusivity m d^-1

Eeff = E0.*exp(-Param{1}.att_sw .* sqrt(maxkappa/Param{1}.mu0));
alpha = Param{1}.alpha_win + 0.5.*(Param{1}.alpha_sum - Param{1}.alpha_win) .*  ...
        (1 + tanh((Eeff - Param{1}.Ecrit)./Param{1}.deltaE));
mu = ((alpha.*PAR./sqrt(alpha.^2.*PAR.^2+Param{1}.mu0.^2)) .* Param{1}.mu0).*Param{1}.Q_P.^(T./10);

mu = timetable(swrad.Time, mu, 'variablenames', {'banas_P'});

% BEST_NPZ

watts2photons = 0.394848; % W m^-2 -> E/m^2/d

DrateS = Param{2}.DiS .* 10.0.^(Param{2}.DpS .* T); % doublings d^-1 (temp dependent doubling rate)
DrateL = Param{2}.DiL .* 10.0.^(Param{2}.DpL .* T); % doublings d^-1

PmaxS = log(2.0) .* DrateS;
PmaxL = log(2.0) .* DrateL;

PmaxSs = PmaxS*Param{2}.ccr;    % mg C (mg chl)^-1 d^-1
PmaxLs = PmaxL*Param{2}.ccrPhL; % mg C (mg chl)^-1 d^-1

LightLimS = tanh(Param{2}.alphaPhS * PAR.*watts2photons./PmaxSs);
LightLimL = tanh(Param{2}.alphaPhL * PAR.*watts2photons./PmaxLs);

mu.bestnpz_PhS = PmaxS .* LightLimS;
mu.bestnpz_PhL = PmaxL .* LightLimL;

% COBALT

cfac = 2.77e-18;
sperd = 86400;

expkT = exp(Param{3}.kappa_eppley .* T);
irr_mem = movmean(PAR, [23 0]);

P_C_m = Param{3}.P_C_max_Sm .* expkT + eps;
theta = (Param{3}.thetamax_Sm - Param{3}.thetamin) ./ (1.0 + Param{3}.thetamax_Sm .* Param{3}.alpha_Sm .* irr_mem.*0.5 ./ P_C_m) + Param{3}.thetamin;
irrlim = 1.0-exp(-Param{3}.alpha_Sm.*PAR .* theta ./ P_C_m);

mu.cobalt_sm = (P_C_m ./ (1.0 + Param{3}.zpllgr) .* irrlim - expkT.*Param{3}.bresp_Sm)*sperd;

P_C_m = Param{3}.P_C_max_Lg .* expkT + eps;
theta = (Param{3}.thetamax_Lg - Param{3}.thetamin) ./ (1.0 + Param{3}.thetamax_Lg .* Param{3}.alpha_Lg .* irr_mem.*0.5 ./ P_C_m) + Param{3}.thetamin;
irrlim = 1.0-exp(-Param{3}.alpha_Lg.*PAR .* theta ./ P_C_m);

mu.cobalt_lg = (P_C_m ./ (1.0 + Param{3}.zpllgr) .* irrlim - expkT.*Param{3}.bresp_Lg)*sperd;

%% ... plot temp exponentials

Tplt = linspace(0, 10, 100);
muT = [...
    Param{1}.mu0 .* Param{1}.Q_P.^(Tplt./10)
    log(2.0).*Param{2}.DiS .* 10.0.^(Param{2}.DpS .* Tplt)
    log(2.0).*Param{2}.DiL .* 10.0.^(Param{2}.DpL .* Tplt)
    exp(Param{3}.kappa_eppley .* Tplt) .* Param{3}.P_C_max_Sm.*sperd
    exp(Param{3}.kappa_eppley .* Tplt) .* Param{3}.P_C_max_Lg.*sperd
    ];

%% Idealized N limitation

[no3,nh4] = meshgrid(linspace(0,2,100), linspace(0,1,100)); % mmol/m^3 (/1025/1000 -> mol/kg)
cvt = 1/1025/1000;

% Banas

Ntot = no3 + Param{1}.phi_NH4 .* nh4;
nlim.banas_P = (Ntot./(Param{1}.kmin + 2.*sqrt(Param{1}.kmin.*Ntot) + Ntot));

% BEST_NPZ

NOLimS = no3./((Param{2}.k1PhS + no3) .* (1.0 + nh4./Param{2}.k2PhS));
NHLimS = nh4./(Param{2}.k2PhS + nh4);

nlim.bestnpz_PhS = NOLimS + NHLimS;

NOLimL = no3./((Param{2}.k1PhL + no3) .* (1.0 + nh4./Param{2}.k2PhL));
NHLimL = nh4./(Param{2}.k2PhL + nh4);

nlim.bestnpz_PhL = NOLimL + NHLimL;

% COBALT

no3lim = (no3*cvt)./(Param{3}.k_no3_Sm + (no3*cvt) + (Param{3}.k_no3_Sm./Param{3}.k_nh4_Sm).*(nh4*cvt));
nh4lim = (nh4*cvt)./(Param{3}.k_nh4_Sm + (nh4*cvt) + (Param{3}.k_nh4_Sm./Param{3}.k_no3_Sm).*(no3*cvt));

nlim.cobalt_sm = no3lim + nh4lim;

no3lim = (no3*cvt)./(Param{3}.k_no3_Lg + (no3*cvt) + (Param{3}.k_no3_Lg./Param{3}.k_nh4_Lg).*(nh4*cvt));
nh4lim = (nh4*cvt)./(Param{3}.k_nh4_Lg + (nh4*cvt) + (Param{3}.k_nh4_Lg./Param{3}.k_no3_Lg).*(no3*cvt));

nlim.cobalt_lg = no3lim + nh4lim;

% Assume base growth rates based on early spring (choose a noon time so
% BESTNPZ LightLim is not a complicating factor)

[~,imin] = min(abs(mu.Time - datetime(yr,5,1,12,0,0)));


% Plot

fld = mu.Properties.VariableNames;
nfld = length(fld);

h = plotgrid('size', [nfld+1 nfld+1], 'sp', 0.02, 'mar', 0.08); %, 'ml', 0.08, 'mt', 0.08);
setpos(h.fig, '# # 800 800');
arrayfun(@(x) set(x, 'position', x.Position+[0.03 0 0 0]), h.ax(:,2:end));  % shift diffs away from main
arrayfun(@(x) set(x, 'position', x.Position+[0 -0.03 0 0]), h.ax(2:end,:));
arrayfun(@(x) set(x, 'position', x.Position+[0 -0.012 0 0]), h.ax([3 5],:)); % shift model groups closer
arrayfun(@(x) set(x, 'position', x.Position+[ 0.012 0 0 0]), h.ax(:,[3 5]));

for ii = 1:nfld
    [~,h.lim(ii,1)] = contour(h.ax(ii+1,1), no3, nh4, nlim.(fld{ii}).*mu.(fld{ii})(imin), 0:0.02:1.4);
    [~,h.lim(ii,2)] = contour(h.ax(1,ii+1), no3, nh4, nlim.(fld{ii}).*mu.(fld{ii})(imin), 0:0.02:1.4);
    for jj = 1:nfld
        [~,h.diff(ii,jj)] = contour(h.ax(ii+1,jj+1), no3, nh4, nlim.(fld{ii}).*mu.(fld{ii})(imin) - nlim.(fld{jj}).*mu.(fld{jj})(imin), -1:0.02:1);
    end
end
arrayfun(@(ax) shading(ax, 'flat'), h.ax);
set(h.ax, 'clim', [0 1.4], 'colormap', crameri('-batlow'), 'fontsize', 8, 'tickdir', 'out');
set(h.ax(2:end,2:end), 'clim', [-1 1], 'colormap', cmocean('balance'));

h.lab1 = labelaxes(h.ax(1,2:end), {'single', 'small', 'large', 'small', 'large'}, 'northoutside');
h.lab2 = labelaxes(h.ax(2:end,1), {'single', 'small', 'large', 'small', 'large'}, 'westoutside', 'rotation', 90, 'horiz', 'center', 'hbuffer', 0.3);
suplabel('axes', h.ax(1,2  ), 'title', 'Banas', 'buffert', 0.02);
suplabel('axes', h.ax(1,3:4), 'title', 'BEST\_NPZ', 'buffert', 0.02);
suplabel('axes', h.ax(1,5:6), 'title', 'COBALT', 'buffert', 0.02);

suplabel('axes', h.ax(2  ,1), 'ylabel', 'Banas',     'buffery', 0.05);
suplabel('axes', h.ax(3:4,1), 'ylabel', 'BEST\_NPZ', 'buffery', 0.05);
suplabel('axes', h.ax(5:6,1), 'ylabel', 'COBALT',    'buffery', 0.05);

set(h.ax(1,1), 'visible', 'off');
set(h.ax(:,[4 6]), 'yticklabel', '');
set(h.ax([3 5],:), 'xticklabel', '');

xlabel(h.ax(1,2), 'NO_3 (mmol N m^{-3})');
ylabel(h.ax(1,2), 'NH_4 (mmol N m^{-3})')

h.cb(1) = colorbar(h.ax(1,3), 'south');
h.cb(2) = colorbar(h.ax(4,3), 'north');
h.cb(1).Position(2) = h.ax(2,1).Position(2)+h.ax(2,1).Position(4)-0.04-h.cb(1).Position(4);
h.cb(2).Position(2) = h.ax(2,1).Position(2)+0.04;
xlabel(h.cb(1), 'Growth rate (d^{-1})');
xlabel(h.cb(2), 'Difference (row - column, d^{-1})');

mask = ~(triu(ones(size(h.ax))));
mask(1,:) = 1;
mask(:,1) = 1;
mask = ~mask;

% mask = logical(diag(ones(1,6)));
delete(h.ax(mask));

exportgraphics(h.fig, 'nlim.png', 'Resolution', 300);


%% ... fratio (depth-integrated)

nbudfilem2 = fullfile(simfolder, 'analysis/nbudgets_all_M2single.mat');
NbudM2 = load(nbudfilem2);

[~,mask] = ismember(["../bgcmip_nbudget_banas/Out", ...
                     "../bgcmip_nbudget_bestnpz/Out", ...
                     "../bgcmip_nbudget_cobalt/Out"], NbudM2.outpath);
NbudM2.G = NbudM2.G(mask);
NbudM2.tg = NbudM2.tg(mask);

plt = {1 'phyto'
       2 'PhS'
       2 'PhL'
       3 'nsm'
       3 'nlg'};
nplt = size(plt,1);

fratio = cell(nplt,1);
for ii = 1:nplt
    Tmp = NbudM2.G{plt{ii,1}}.Edges;
    ttmp = NbudM2.tg{plt{ii,1}};
    isss = strcmp(Tmp.EndNodes, plt{ii,2});

    Tmp.val(isss(:,1)) = cellfun(@(x) x*-1, Tmp.val(isss(:,1)), 'uni', 0);
    Tmp = Tmp(any(isss,2),:);

    isnh4 = ismember(lower(Tmp.type), {'gpp','npp'}) & strcmpi(Tmp.EndNodes(:,1),'nh4');
    isno3 = ismember(lower(Tmp.type), {'gpp','npp'}) & strcmpi(Tmp.EndNodes(:,1),'no3');

    fratio{ii} = timetable(ttmp, ...
        max(Tmp.val{isno3},0)./(max(Tmp.val{isno3},0) + max(Tmp.val{isnh4},0)));

end

%% ... plot limitation functions and results at M2

% Colormap for limitation functions

sz = [0 0 1 2 1 2 1 2];
sn = [2 1 0 0 1 1 2 2];
coltmp = colmodel([1 1 2 2 3 3 3 3],:);
for ii = 1:size(coltmp,1)
    if sn(ii) == 1
        coltmp(ii,:) = adjustcol(coltmp(ii,:), 0.5);
    elseif sn(ii) == 2
        coltmp(ii,:) = adjustcol(coltmp(ii,:), -0.2);
    end
end

% Limitation functions

h = plotgrid('size', [2 3], 'sv', 0.1, 'mr', 0.15); %, 'mr', 0.05, 'ml', 0.05, 'mb', 0.15);
% setpos(h.fig, '# # 9in 4in');

for ii = 1:3
    h.ln = plot(h.ax(1,ii), Ilim.irr, Ilim.muplt(:,:,ii));
    
    set(h.ln, {'color'}, num2cell(coltmp,2), 'linewidth', 2);
    set(h.ln(sz==2), 'linestyle', ':');

    xlabel(h.ax(1,ii), 'PAR (W/m^2)');
    ylabel(h.ax(1,ii), 'Growth rate (d^{-1})');
end

h.leg = legendflex(h.ln, {'Banas, winter', 'Banas, summer', ...
       'BEST_NPZ, small', 'BEST_NPZ, large', ...
       'COBALT, small, max chl:C', 'COBALT, large, max chl:C', 'COBALT, small, min chl:C', 'COBALT, large, min chl:C'}, ...
       'interpreter', 'none', 'anchor', {'n','s'}, 'buffer', [0 5], 'nrow', 2, 'box', 'off', 'xscale', 0.5, 'ref', h.ax(1,2));
set(h.fig, 'color', 'w');

set(h.ax, 'ylim', [0 3], 'xlim', minmax(Ilim.irr), 'fontsize', 8);

labelaxes(h.ax(1,:), compose("%d\\circC", Ilim.temp), 'northwest');

% Timeseries

tsmooth = 24*5;

h.ax2 = mergeaxes(h.ax(2,:));
h.ts = plot(h.ax2, mu.Time, movmean(mu{:,:}, tsmooth));

set(h.ts, {'color'}, num2cell(colmodel([1 2 2 3 3],:),2));
set(h.ts([3 5]), 'linestyle', ':');

h.ax2(2) = axes('position', h.ax2(1).Position);

h.par = area(h.ax2(2), swrad.Time, movmax(PAR,tsmooth), 'facecolor', rgb('gold'), 'edgecolor', 'none', 'facealpha', 0.1);
h.ax2(3) = axes('position', h.ax2(1).Position);
h.temp = plot(h.ax2(3), stemp.Time, movmean(T,tsmooth), 'color', rgb('gray'));
h.oax = offsetaxis(h.ax2(3), 'y', 0.1, 'yloc', 'r');

set(h.ax2, 'color', 'none', 'box', 'off', 'fontsize', 8, 'tickdir', 'out');
set(h.ax2(2), 'yaxisloc', 'right', 'xaxisloc', 'top', 'xticklabel', '', 'ycolor', rgb('gold'), 'xcolor', 'k', 'layer', 'top', 'ylim', [0 400]);
set(h.ax2(3), 'visible', 'off');
set(h.oax, 'tickdir', 'out', 'fontsize', 8, 'ycolor', rgb('gray'), 'yticklabelmode', 'auto', 'ytick', -2:2:16);

uistack(h.ax2(1), 'top');

ylabel(h.ax2(1), 'Growth rate (d^{-1})');
ylabel(h.ax2(2), 'PAR (W m^{-2})');
ylabel(h.oax, 'Temperature (\circC)')

set([h.ts; h.temp], 'linewidth', 1.5);

% f-ratio

% h.ax3 = mergeaxes(h.ax(3,:));
% 
% for ii = 1:length(fratio)
%     h.fr(ii) = plot(h.ax3, fratio{ii}.ttmp, fratio{ii}.Var1);
%     hold(h.ax3, 'on');
% end
% 
% colnet = cellfun(@(x) colmodel(x,:), plt(:,1), 'uni', 0);
% set(h.fr(1:5), {'color'}, colnet);
% set(h.fr([3 5]), 'linestyle', ':');
% set(h.ax3, 'xlim', datetime([yr yr+1],1,1));

exportgraphics(h.fig, 'lightlim.png', 'Resolution', 300);


%% ... all fluxes at M2 (flux per group figures)

relflag = true;

nbudfilem2 = fullfile(simfolder, 'analysis/nbudgets_all_M2single.mat');
NbudM2 = load(nbudfilem2);

[~,mask] = ismember(["../bgcmip_loop_nbudget_banas/Out", ...
                     "../bgcmip_loop_nbudget_bestnpz/Out", ...
                     "../bgcmip_loop_nbudget_cobalt/Out"], NbudM2.outpath);
lblname = ["Banas", "BEST_NPZ", "COBALT"];
NbudM2.G = NbudM2.G(mask);
NbudM2.tg = NbudM2.tg(mask);


plt = {1 'phyto'   1 'microzoo'  NaN ''
       2 'PhS'     2 'MZL'         2 'Cop'
       2 'PhL'   NaN ''            2 'NCaS'
       NaN ''    NaN ''            2 'EupS'
       NaN ''    NaN ''            2 'Jel'
       3 'nsm'     3 'nsmz'        3 'nmdz'
       3 'nlg'   NaN ''            3 'nlgz'
       };

% plt = {
%        2 'MZL'
%        3 'nsmz'};

nplt = size(plt,1);
ncol = size(plt,2)/2;

h = plotgrid('size', [nplt+1 ncol], 'sp', 0.05, 'mar', 0.05);
setpos(h.fig, '# # 22in 17in');

% Total fluxes in/our of phyto groups

adj = [0.1 -0.1 0.2 -0.2 0.3 -0.3 0.4 -0.4 0.5]*2; % colors

for ii = 1:nplt
    for jj = 1:ncol
        if isnan(plt{ii,2*jj-1})
            plot(h.ax(ii,jj), NaT, NaN);
        else

            Tmp = NbudM2.G{plt{ii,2*jj-1}}.Edges;
            ttmp = NbudM2.tg{plt{ii,2*jj-1}};
            isss = strcmp(Tmp.EndNodes, plt{ii,2*jj});
        
            Tmp.val(isss(:,1)) = cellfun(@(x) x*-1, Tmp.val(isss(:,1)), 'uni', 0);
            Tmp = Tmp(any(isss,2),:);

            isp = strcmp(NbudM2.G{plt{ii,2*jj-1}}.Nodes.Name, plt{ii,2*jj});
            pbio = NbudM2.G{plt{ii,2*jj-1}}.Nodes.val{isp};

            Tmp.cplt = num2cell(ecmap{Tmp.type,'rgb'},2);
            [unq,~,iunq] = uniquecell(Tmp.cplt);
  
            for it = 1:height(Tmp)
                Tmp.npertp(it) = sum(iunq(1:it)==iunq(it));
                Tmp.cplt{it} = adjustcol(Tmp.cplt{it}, adj(Tmp.npertp(it)));
            end
        
            % Fluxes
            
            axes(h.ax(ii,jj));
            if relflag
                h.ba(ii,jj) = barareaneg(ttmp, cell2mat(Tmp.val')./pbio, 'bar');
            else
                h.ba(ii,jj) = barareaneg(ttmp, cell2mat(Tmp.val'), 'bar');
            end
            set(h.ba(ii,jj).pos, {'facecolor'}, Tmp.cplt);
            set(h.ba(ii,jj).neg, {'facecolor'}, Tmp.cplt);
        
            lbl = cellfun(@(a,b,c) strrep(sprintf('%s, %s\\rightarrow%s',a,b,c), '_', '\_'), Tmp.type, Tmp.EndNodes(:,1), Tmp.EndNodes(:,2), 'uni', 0);
            h.leg(ii,jj) = legendflex(h.ba(ii,jj).pos, lbl, 'ref', h.ax(ii,jj), 'anchor', {'nw','sw'}, 'buffer', [0 0], ...
                'fontsize', 8, 'xscale', 0.5, 'nrow', 3, 'box','off');
        
            labelaxes(h.ax(ii,jj), lblname{plt{ii,2*jj-1}}+": "+string(plt{ii,2*jj}), 'northwest', 'interpreter', 'none', 'fontweight', 'b');

            % Biomass
        
            h.bio(ii,jj) = plot(h.ax(nplt+1,jj), ttmp, pbio, 'color', colmodel(plt{ii,jj*2-1},:));
            hold(h.ax(nplt+1,jj), 'on')

    % % Net production
    % 
    % isnet = ismember(lower(Tmp.type), {'npp', 'gpp', 'exu', 'ege', 'exc', 'res'});
    % h.net(ii) = plot(h2.ax(1), ttmp, sum(cell2mat(Tmp.val(isnet)'),2));
    % hold(h2.ax(1), 'on');
    % 
    % % f-ratio
    % 
    % isnh4 = ismember(lower(Tmp.type), {'gpp','npp'}) & strcmpi(Tmp.EndNodes(:,1),'nh4');
    % isno3 = ismember(lower(Tmp.type), {'gpp','npp'}) & strcmpi(Tmp.EndNodes(:,1),'no3');
    % fratio = Tmp.val{isno3}./(Tmp.val{isno3} + Tmp.val{isnh4});
    % h.frat(ii) = plot(h2.ax(2), ttmp, fratio);
    % hold(h2.ax(2), 'on');

        end
    end
end

% yrplt = 1990;
% set(h.ax, 'xlim', datetime(yrplt, [3 15], 1)); %, 'ylim', [-50 50]);
set(h.ax, 'xlim', [ttmp(1) datetime(1992,7,1)], ...
    'xtick', datetime(1990,4:12:36,1), ...
    'xgrid', 'on', ...
    'xminortick', 'on', ...
    'xminorgrid', 'on');
arrayfun(@(x) set(x.XAxis, 'MinorTickValues', datetime(1990,1:36,1)), h.ax);

set([h.ba.sum], 'barwidth', 1, 'edgecolor', 'none', 'facecolor', 'k', 'facealpha', 0.5);
for ii = 1:length(h.ba)
    for jj = 1:ncol
        if ~isnan(plt{ii,2*jj-1})
            plot(h.ax(ii,jj), h.ba(ii,jj).sum.XData, h.ba(ii,jj).sum.YData, 'k');
        end
    end
end
set([h.ba.sum], 'visible', 'off');

set([h.ba.pos h.ba.neg], 'barwidth', 1);

set(h.ax(1:nplt,1), 'ylim', [-1 1]);
set(h.ax(1:nplt,2), 'ylim', [-2 2]);

set(h.ax(nplt+1,:), 'ylim', [1e-2 400], 'yscale', 'log');

set(h.bio(ismember(plt(:,2:2:end),{'PhL','nlg','Jel','nlgz'})), 'linestyle', '-.');
set(h.bio(ismember(plt(:,2:2:end),{'NCaS'})), 'linestyle', ':');
set(h.fig, 'color', 'w');

set(h.ax(1:nplt-1,:), 'xticklabel', '');
set(h.ax([cellfun(@isnan, plt(:,1:2:end)); false(1,ncol)]), 'visible', 'off');

for jj = 1:ncol
    mask = ismissing(plt(:,jj*2));
    h.leg(jj) = legendflex(h.bio(~mask,jj), plt(~mask,jj*2), ...
        'ref', h.ax(nplt+1,jj), 'anchor', {'nw','sw'}, 'buffer', [0 0], ...
        'fontsize', 8, 'xscale', 1.0, 'nrow', 1, 'box','off');
end

% export_fig(h.fig, "flx_per_biomass_spinup", '-png', '-r300', '-nocrop', '-painters');
export_fig(h.fig, "flx_per_biomass", '-png', '-r300', '-nocrop', '-painters');


% set(h.ax(1:nplt-1), 'xticklabel', '');
% 
% pflag = any(strcmp(plt(:,2), 'phyto'));
% mzflag = any(strcmp(plt(:,2), 'microzoo'));
% if relflag
%     if pflag
%         set(h.ax, 'ylim', [-1 1]);
%     elseif mzflag
%         set(h.ax, 'ylim', [-2 2]);
%     end
% else
%     set(h.ax, 'ylim', [-50 50]);
% end
% 
% if pflag
%     set(h.ax(nplt+1), 'yscale', 'log', 'ylim', [1 400]);
% end
% if mzflag
%     set(h.ax(nplt+1), 'ylim', [1e-2 30], 'yscale', 'log');
% end
% 
% set(h.bio(ismember(plt(:,2:2:end),{'PhS','nsm','Jel'})), 'linestyle', '-.');
% set(h.bio(ismember(plt(:,2:2:end),{'NCaS'})), 'linestyle', ':');
% 
% set(h.fig, 'color', 'w');

% if pflag
%     export_fig(h.fig, "flx_phyto"+yrplt, '-png', '-r150', '-nocrop');
% elseif mzflag
%     export_fig(h.fig, "flx_microzoo"+yrplt, '-png', '-r150', '-nocrop');
% end


% % Plot net
% 
% colnet = cellfun(@(x) colmodel(x,:), plt(:,1), 'uni', 0);
% % colnet([3 5]) = cellfun(@(x) adjustcol(x, 0.5), colnet([3 5]), 'uni', 0);
% % colnet([2 4]) = cellfun(@(x) adjustcol(x, 0.2), colnet([3 5]), 'uni', 0);
% 
% set(h.net(1:5), {'color'}, colnet);
% set(h.net([3 5]), 'linestyle', ':');
% % set(h.net([2 4]), 'linestyle', '-.');
% h.net(6) = plot(h2.ax(1), h.net(2).XData, h.net(2).YData+h.net(3).YData, 'color', adjustcol(colmodel(2,:), 0.8));
% h.net(7) = plot(h2.ax(1), h.net(4).XData, h.net(4).YData+h.net(5).YData, 'color', adjustcol(colmodel(3,:), 0.8));
% set(h.net, 'linewidth', 1.5);
% 
% set(h.frat(1:5), {'color'}, colnet);
% set(h.frat([3 5]), 'linestyle', ':');
% 
% 
% % plot(h2.ax(2), NaT, NaN)
% % h.net(2,:) = copyobj(h.net(1,:), h2.ax(2));
% set(h2.ax, 'xlim', datetime(1993, [4 8], 1));
% set(h2.ax(2), 'ylim', [-0.1 1.1]);
% % set(h2.ax(2), 'xlim', datetime(1993, [3 7], 1));

%% ... plot log-scaled biomass across models

vplt = {...
    ["phyto", "microzoo", "det_small", "det_large"]
    ["PhS", "PhL", "MZL", "Cop", "NCaS", "NCaO", "EupS", "EupO", "Jel", "Det", "DetF"]
    ["ndi", "nsm", "nlg", "nsmz", "nmdz", "nlgz", "ndet", "ldon", "sldon", "srdon"]
    };

h = plotgrid('size', [3 1]);

for ii = 1:3
    tmp = timetable(NbudM2.tg{ii}, NbudM2.G{ii}.Nodes.val{:}, 'variablenames', NbudM2.G{ii}.Nodes.Name);

    [~,loc] = ismember(vplt{ii}, NbudM2.G{ii}.Nodes.Name);
    tmpcmap = ncmap{NbudM2.G{ii}.Nodes.type(loc),'rgb'};

    hold(h.ax(ii), 'on');
    for iv = 1:length(vplt{ii})
        plot(h.ax(ii), tmp, vplt{ii}{iv}, 'color', tmpcmap(iv,:));
    end

end

set(h.ax, 'ylim', [1e-5 300], 'yscale', 'log');
set(h.ax, 'xlim', [ttmp(1) datetime(1994,1,1)], ...
    'xtick', datetime(1990,4:12:48,1), ...
    'xgrid', 'on', ...
    'xminortick', 'on', ...
    'xminorgrid', 'on');
arrayfun(@(x) set(x.XAxis, 'MinorTickValues', datetime(1990,1:48,1)), h.ax);

%% ... same, but with consolidated groups (model x loop/primary panels)

nbudfilem2 = fullfile(simfolder, 'analysis/nbudgets_all_M2single.mat');
NbudM2 = load(nbudfilem2);

[~,mask] = ismember(["../bgcmip_nbudget_banas/Out", ...
                     "../bgcmip_nbudget_bestnpz/Out", ...
                     "../bgcmip_nbudget_cobalt/Out" ...
                     "../bgcmip_loop_nbudget_banas/Out", ...
                     "../bgcmip_loop_nbudget_bestnpz/Out", ...
                     "../bgcmip_loop_nbudget_cobalt/Out"], NbudM2.outpath);
NbudM2.G = NbudM2.G(mask);
NbudM2.tg = NbudM2.tg(mask);

vplt = {...
    []    "phyto"       "microzoo" []
    "PhS" "PhL"         "MZL"      ["Cop","NCaS","NCaO","EupS","EupO","Jel"]
    "nsm" ["nlg" "ndi"] "nsmz"     ["nmdz", "nlgz"]
    };

h = plotgrid('size', [3 2], 'mar', 0.05, 'sv', 0.02, 'mb', 0.1, 'mr', 0.01, 'ml', 0.1, ...
    'collabel', {'Primary simulation', 'Spinup simulation'}, 'rowlabel', {'Banas', 'BEST\_NPZ', 'COBALT'}, ...
    'rowlabeloffset', 0.08, 'collabeloffset', 0.03);
h.ax = fliplr(h.ax);
setpos(h.fig, '# # 8.5in 6in');

for ii = 1:3 % models
    for jj = 1:2 % sims (spinup vs loop)

        sidx = (jj-1)*3+ii;
        tmp = timetable(NbudM2.tg{sidx}, NbudM2.G{sidx}.Nodes.val{:}, 'variablenames', NbudM2.G{sidx}.Nodes.Name);

        vval = nan(height(tmp),4);
        for iv = 1:4
            if ~isempty(vplt{ii,iv})
                vval(:,iv) = sum(tmp{:,vplt{ii,iv}},2);
            end
        end
        
        h.ln(:,ii,jj) = plot(h.ax(ii,jj), tmp.Time, vval);

    end
end

% Prettify

set(h.ax, 'ylim', [1e-3 300], 'yscale', 'log');
set(h.ax, 'xlim', [ttmp(1) datetime(1994,1,1)], ...
    'xtick', datetime(1990,1:3:48,1), ...
    'xgrid', 'on', ...
    'xminortick', 'on', ...
    'xminorgrid', 'on');
arrayfun(@(x) set(x.XAxis, 'MinorTickValues', datetime(1990,1:48,1)), h.ax);
set(h.ax, 'fontsize', 8, 'tickdir', 'out');
set(h.ax(1:end-1,:), 'xticklabel', '');

cmapbio = cptcmap('Paired_12');
set(h.ln(1,:,:), 'color', cmapbio(3,:));
set(h.ln(2,:,:), 'color', cmapbio(4,:));
set(h.ln(3,:,:), 'color', cmapbio(1,:));
set(h.ln(4,:,:), 'color', cmapbio(2,:));
set(h.ln, 'linewidth', 1.5);

% title(h.ax(1,1), 'Spinup simulation', 'fontsize', 10);
% title(h.ax(1,2), 'Primary simulation', 'fontsize', 10);

ylabel(h.ax(1,2), "Biomass (mmol N m^{-2})")
h.leg = legend(h.ln(:,3,2), ["Small phyto.", "Large phyto.", "Microzoo.", "Mesozoo."], 'location', 'best');


exportgraphics(h.fig, 'biomass_pz_logscaled.png', 'resolution', 150);

%% ... same, but with consolidated groups (P/Z x loop/primary panels)

% nbudfilem2 = fullfile(simfolder, 'analysis/nbudgets_all_M2single.mat');
% NbudM2 = load(nbudfilem2);
% 
% [~,mask] = ismember(["../bgcmip_nbudget_banas/Out", ...
%                      "../bgcmip_nbudget_bestnpz/Out", ...
%                      "../bgcmip_nbudget_cobalt/Out" ...
%                      "../bgcmip_loop_nbudget_banas/Out", ...
%                      "../bgcmip_loop_nbudget_bestnpz/Out", ...
%                      "../bgcmip_loop_nbudget_cobalt/Out"], NbudM2.outpath);
% NbudM2.G = NbudM2.G(mask);
% NbudM2.tg = NbudM2.tg(mask);

vplt = {...
    []    "phyto"       "microzoo" []
    "PhS" "PhL"         "MZL"      ["Cop","NCaS","NCaO","EupS","EupO","Jel"]
    "nsm" ["nlg" "ndi"] "nsmz"     ["nmdz", "nlgz"]
    };

h = plotgrid('size', [3 2], 'mar', 0.05, 'sv', 0.02, 'mb', 0.1, 'mr', 0.01, 'ml', 0.1, ...
    'collabel', {'Primary simulation', 'Spinup simulation'}, 'rowlabel', {'Phytoplankton', 'Microzooplankton', 'Mesozooplankton'}, ...
    'rowlabeloffset', 0.08, 'collabeloffset', 0.03);
h.ax = fliplr(h.ax);
setpos(h.fig, '# # 8.5in 6in');

arrayfun(@(x) hold(x,'on'), h.ax);


for ii = 1:3
    for jj = 1:2 % sims (spinup vs loop)

        sidx = (jj-1)*3+ii;
        tmp = timetable(NbudM2.tg{sidx}, NbudM2.G{sidx}.Nodes.val{:}, 'variablenames', NbudM2.G{sidx}.Nodes.Name);

        vval = nan(height(tmp),4);
        for iv = 1:4
            if ~isempty(vplt{ii,iv})
                vval(:,iv) = sum(tmp{:,vplt{ii,iv}},2);
            end
        end

        % h.ln(1:2,ii,jj) = plot(h.ax(1,jj), tmp.Time, vval(:,1:2), 'color', adjustcol(colmodel(ii,:),0.5));
        h.ln(1,ii,jj) = plot(h.ax(1,jj), tmp.Time, sum(vval(:,1:2),2,'omitnan'), 'color', colmodel(ii,:));
        h.ln(2,ii,jj) = plot(h.ax(2,jj), tmp.Time, vval(:,3), 'color', colmodel(ii,:));
        h.ln(3,ii,jj) = plot(h.ax(3,jj), tmp.Time, vval(:,4), 'color', colmodel(ii,:));

        % set(h.ln(1,ii,jj), 'linestyle', '-.');
        % set(h.ln(2,ii,jj), 'linestyle', '--');

        % 
        % return
        % 
        % 
        % h.ln(:,ii,jj) = plot(h.ax(ii,jj), tmp.Time, vval);
    end
end

% Prettify

set(h.ax, 'ylim', [1e-3 300], 'yscale', 'log');
set(h.ax, 'xlim', [ttmp(1) datetime(1994,1,1)], ...
    'xtick', datetime(1990,1:3:48,1), ...
    'xgrid', 'on', ...
    'xminortick', 'on', ...
    'xminorgrid', 'on');
arrayfun(@(x) set(x.XAxis, 'MinorTickValues', datetime(1990,1:48,1)), h.ax);
set(h.ax, 'fontsize', 8, 'tickdir', 'out');
set(h.ax(1:end-1,:), 'xticklabel', '');

ylabel(h.ax(1,2), "Biomass (mmol N m^{-2})");
h.leg = legend(h.ln(3,1:3,2), ["Banas", "BEST\_NPZ", "COBALT"], 'location', 'best');

exportgraphics(h.fig, 'biomass_pz_logscaled_2.png', 'resolution', 150);


%% Are ammonium levels different in the basin?


fname = {...
fullfile(simfolder, 'bgcmip_loop_banas', 'Out', 'bgcmip_loop_banas_avg_00001.nc')
fullfile(simfolder, 'bgcmip_loop_bestnpz', 'Out', 'bgcmip_loop_bestnpz_avg_00001.nc')
fullfile(simfolder, 'bgcmip_loop_cobalt', 'Out', 'bgcmip_loop_cobalt_avg_00001.nc')
};

amm = cell(3,1);
amm{1} = ncread(fname{1}, 'NH4');
amm{2} = ncread(fname{2}, 'NH4');
amm{3} = ncread(fname{3}, 'nh4');

ammavg = cellfun(@(x) mean(x,4), amm, 'uni', 0);

nz = 30;

[zr,zw] = calcromsz(Grd.h, 0, nz);

zplt = [-500 -100 -50 -5];
zidx = ones(size(Grd.h)).*permute(1:nz,[3 1 2]);

[nxi,neta] = size(Grd.h);
[xx,ee] = ndgrid(1:nxi, 1:neta);

for iz = 1:length(zplt)

    tmp = zidx;
    tmp(zw(:,:,2:end)<zplt(iz)) = NaN;
    tmp = bottom(tmp(:,:,end:-1:1));

    ind = sub2ind([nxi neta nz], xx,ee,tmp);

    h = plotgrid('size', [1 3], 'sp', 0, 'mar', 0.01);
    h.fig.Position(3:4) = [1220 400];
    % setpos(h.fig, '# # ')
    for iax = 1:3
        val = ammavg{iax}(ind);
        val(zw(:,:,1)>zplt(iz)) = NaN;

        axes(h.ax(iax));
        plotromsrho(Grd, val);
    end
    set(h.ax, 'clim', [0 10])

end

%% Timing of spring bloom

%% ... ice cover (for organizing by early/late)

% M2, M8 ice cover

% outfile = dir(fullfile(simfolder, 'bgcmip_loop_banas', 'Out', '*avg*'));
% outfile = fullfile({outfile.folder}, {outfile.name});
% 
% taice = ncdateread(outfile, 'ocean_time');
% Scs = struct('xi_rho', [Station.xi(2) 1 1], 'eta_rho', [Station.eta(2) 1 1]);
% aice = ncstruct(outfile, 'aice', Scs);
% Scs = struct('xi_rho', [Station.xi(6) 1 1], 'eta_rho', [Station.eta(6) 1 1]);
% aice2 = ncstruct(outfile, 'aice', Scs);
% 
% aice = timetable(taice, squeeze(aice.aice), squeeze(aice2.aice));
% 
% save aice_m2m8 aice;


% Regrid to year x doy

[Ice.grd, Ice.yr, Ice.tmid] = reshapetimeseries(aice.taice, aice{:,1}, 'bin', 52);
Ice.grd(:,:,2) = reshapetimeseries(aice.taice, aice{:,2}, 'bin', 52);
% Ice.grd(Ice.grd <= 0) = NaN;

% Date of retreat

yy = unique(year(aice.taice));

Ice.retreat = nan(length(yy),2);
for iy = 1:length(yy)
    for ii = 1:2

        [b,n,bi] = RunLength(Ice.grd(:,iy,ii)>0.1);
        n = n(b==1);
        bi = bi(b==1);
    
        if ~isempty(n)
            if isscalar(n)
                Ice.retreat(iy,ii) = doy(Ice.tmid(bi+n));
            else
                idx = find(Ice.tmid(bi) < datetime(1990,7,1),1,'last');
                Ice.retreat(iy,ii) = doy(Ice.tmid(bi(idx)+n(idx)));
            end
        end
    end
end

Ice.grd(Ice.grd <= 0) = NaN;

for ii = 1:2
    [~,Ice.srt(:,ii)] = sort(Ice.retreat(:,ii), 'missingplacement', 'first');
end

[Ice.yrg, Ice.doy] = ndgrid(Ice.yr, doy(Ice.tmid));


% [aicegrd, yyice, tmidice] = reshapetimeseries(aice.taice, aice{:,1}, 'bin', 52);
% % aicegrd = aicegrd';
% [aicegrd(:,:,2), yyice, tmidice] = reshapetimeseries(aice.taice, aice{:,1}, 'bin', 52);
% aicegrd(aicegrd <= 0) = NaN;
% [yyice, tmidice] = ndgrid(yyice, doy(tmidice));




% return
% 
% [aice, taice] = deal(cell(length(outfile),1));
% 
% 
% w = Grd.area_feast./sum(Grd.area_feast);
% w = w(:);
% 
% for ii = 1:length(outfile)
%     aice{ii} = nansum(reshape(ncread(outfile{ii}, 'aice'), nxi*neta,[]).*w, 1)';
%     taice{ii} = ncdateread(outfile{ii}, 'ocean_time');
% end
% 
% aice = timetable(cat(1, taice{:}), cat(1, aice{:}));
% 
% aicemax = retime(aice, 'yearly', 'mean');

%% ... Plot NPP vs year, doy

% M2 = 2, M8 = 6

midx = [2 6];

h = plotgrid('size', [3 2], 'mb', 0.05, 'mt', 0.08, ...
    'collabel', {'Southeastern Bering Sea (M2)', 'Northern Bering Sea (M8)'}, ...
    'rowlabel', {'Banas', 'BEST\_NPZ', 'COBALT'}, 'rowlabeloffset', 0.08, 'collabeloffset', 0.03);
setpos(h.fig, '# # 8.5in 8in');

for im = 1:length(midx)

    for ii = 1:nbgc
        tmp = Metrics(midx(im),2).(bgcname{ii});

        tmp = retime(tmp, 'daily', 'nearest');
        [xx,yy,tmid] = reshapetimeseries(tmp.Time, tmp.npp, 'bin', 365);

        pcolorpad(h.ax(ii,im), doy(tmid), 1:length(yy)+1, xx(:,Ice.srt(:,im))');
        set(h.ax(ii,im), 'ytick', (1:length(yy))+0.5, 'yticklabel', string(yy(isrt)));
        hold(h.ax(ii,im), 'on');
        % plot(h.ax(ii,im), doysrt, (1:length(yy))+0.5, 'marker', '|', 'color', rgb('medium blue'), 'linestyle', 'none');
        
        [yplt, xplt] = ndgrid((1:length(yy))+0.5, doy(Ice.tmid));
        isz = discretize(reshape(Ice.grd(:,Ice.srt(:,im),im)',[],1), max(0.01,0:.1:0.5));
        scatter(h.ax(ii,im), xplt(:), yplt(:), isz*2, rgb('blue'));

    end
end

arrayfun(@(x) shading(x, 'flat'), h.ax);
tk = datetime(1990,1:12,1);
set(h.ax, 'colormap', cmocean('speed'), 'clim', [0 50], 'xtick', doy(tk), ...
    'xticklabel', datestr(tk, 'm'), 'fontsize', 6, 'layer', 'top', 'tickdir', 'out');

h.cb = colorbar(h.ax(end), 'east');
h.cb.Position(1) = 0.95;
ylabel(h.cb, 'NPP (mmol N m^-3)')


exportgraphics(h.fig, 'bloomvariability.png', 'resolution', 150);

%% Double-check compilation versions

LogFile = dir(fullfile(simfolder, '*', 'Log', '*_log.txt'));
logname = fullfile({LogFile.folder}, {LogFile.name})';
gitrevis = cell(size(logname));
for ii = 1:length(LogFile)
    fid = fopen(logname{ii});
    flag = 0;
    while ~flag
        txt = fgetl(fid);
        flag = contains(txt, 'GIT Revision');
    end
    fclose(fid);
    gitrevis{ii} = txt;
end








