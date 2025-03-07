% Read in data for comparable metrics

bgcmip_setup;

% sims = ["banas", "bestnpz", "cobalt"];
% nsim = length(sims);
% 
% simdir = fullfile('..', "bgcmip_loop"+sims);  
% 
% grdfile = fullfile(moxdir, 'ROMS_Datasets', 'grids', 'AlaskaGrids_Bering10K.nc');
% Grd = ncstruct(grdfile);
% 
% nlayer = 30;
% 
% berapp = fullfile(moxdir, 'kearney', 'bering-Apps');
% pfile = fullfile(berapp, 'Apps', 'Bering_BGC_variants', "bering_bpar_"+sims+".yaml");
% Param = cellfun(@(x) yaml.loadFile(x), pfile, 'uni', 0);
% % save bgcmip_params Param;
% 
% % Param = load('bgcmip_params.mat');
% % Param = Param.Param;
% 
% %% Locations to extract stuff
% 
% %% ... Mask: shelf-only
% 
% warning('off', 'MATLAB:polyshape:repairedBySimplify');
% 
% M = load(fullfile(moxdir, 'bering10k', 'output', 'forecasts', 'cmip6', 'analysis_kak', 'nbudget_masks.mat'));
% mask = nan(size(Grd.h));
% for ii = 1:size(M.masks,1)
%     mask(M.masks{ii,2}) = ii;
%     [ln, lt] = mask2poly(Grd.lon_psi, Grd.lat_psi, M.masks{ii,2}(2:end-1,2:end-1));
%     pmask{ii} = polyshape(ln, lt);
% end
% pshelf = union([pmask{1:8}]);
% 
% shelfmask = reshape(pshelf.isinterior(Grd.lon_rho(:), Grd.lat_rho(:)), size(Grd.h)) & Grd.mask_rho==1;
% [nxi, neta] = size(Grd.h);
% nshelf = nnz(shelfmask);
% 
% %% ... Mask: transect plots (eliminate anything south of Aleutians)
% 
% pland = struct('lat', {cell(2,1)}, 'lon', {cell(2,1)});
% [pland.lat{1}, pland.lon{1}] = borders('alaska');
% [pland.lat{2}, pland.lon{2}] = borders('russia');
% pland = polyshape(wrapTo360([pland.lon{:}]), [pland.lat{:}]);
% ltbox = [50.5 56.76];
% lnbox = [163 196];
% aleut = intersect(polyshape(lnbox([1 1 2 2 1]), ltbox([1 2 2 1 1])), pland);
% % aleut = sortregions(aleut, 'centroid', 'ascend');
% tmp = regions(aleut);
% [clon,clat] = centroid(tmp);
% [~,imin] = min((clat - 56.594).^2 + (clon - 190.52).^2); % remove Pribs (in bounding box)
% aleut = rmboundary(aleut, imin);
% 
% avert = cat(1, aleut.Vertices);
% avert = unique(avert, 'rows');
% avert = avert(~any(isnan(avert),2),:);
% a = alphaShape(avert(:,1), avert(:,2));
% a.Alpha = 10;
% [~, apt] = boundaryFacets(a);
% paleutstrip = polybuffer(polyshape(apt(:,1), apt(:,2)), 0.1);
% 
% pbasin = subtract(pmask{10}, polybuffer(polyshape(apt(:,1), apt(:,2)), 0.1));
% pbasin = regions(sortregions(pbasin, 'area', 'descend'));
% pbasin = pbasin(1);
% 
% pmain = union(pbasin, pshelf);
% 
% mask = double(~reshape(isinterior(paleutstrip, Grd.lon_rho(:), Grd.lat_rho(:)), size(Grd.h)));
% % mask = double(reshape(isinterior(pmain, Grd.lon_rho(:), Grd.lat_rho(:)), size(Grd.h)));
% mask(mask == 0) = NaN;
% 
% %% ... Stations
% 
% statfile = fullfile(moxdir, 'kearney/runscripts/bering_stations_updated.csv');
% Station = readtable(statfile);
% 
% % Just EBS ones for now
% 
% onebs = pshelf.isinterior(wrapTo360(Station.lon), Station.lat);
% Station = Station(onebs,:);
% 
% % Just a few select ones (ECCWO poster)
% 
% % samplestat = ["BS-2", "BS-8", "14", "2 low"];
% % [~,loc] = ismember(samplestat, Station.name);
% % Station = Station(loc,:);
% 
% [Station.xi, Station.eta] = ind2sub(size(Grd.h), Station.gridcell);

%% Extract comparable metrics (new version using post-processed data)

% All versions: initial 30-year, looped 30-year, plus burial variants
 
suffix = ["_CFS", "_CFS_loop", "bury_CFS", "bury_CFS_loop", "burynoinfauna_CFS", "burynoinfauna_CFS_loop"];

% Make robust to restart

nsuffix = length(suffix);
nstation = height(Station);

outfile = 'bgcmip_metrics_at_stations.mat';

canread = true(3, nstation, nsuffix);
canread(1,:,contains(suffix,"bury")) = false; % No bury variants for banas
canread(3,:,contains(suffix,"noinfauna")) = false; % No noinfuauna variant for cobalt

if exist(outfile, 'file')
    M = load(outfile);
    Metrics = M.Metrics;
    Metrics = orderfields(Metrics, {'banas', 'bestnpz', 'cobalt'});
    
    hasdata = ~cellfun(@isempty, struct2cell(Metrics));
else
    hasdata = false(3, nstation, nsuffix);
end

readflag = canread & ~hasdata;

for isuf = nsuffix:-1:1
    for ig = nstation:-1:1

        fprintf('%s, Station %d/%d\n', suffix{isuf}, ig, height(Station));
        % if ig == nstation && isuf == nsuffix
        %     sample = dir(fullfile(moxdir, 'roms_for_public', 'B10K-BGCMIPbanas_CFS_loop', 'Level2', '*_average_phyto_integrated.nc'));
        %     t = ncdateread(fullfile({sample.folder}, {sample.name}), 'ocean_time');
        %     nt = length(t);
        % end
        
        Scs = struct('xi_rho', [Station.xi(ig) 1 1], 'eta_rho', [Station.eta(ig) 1 1]);
    
        % Banas model (uses mmol N m^-3 d^-1)

        if readflag(1,ig,isuf)

            disp("  "+sims(1));

            filebase = fullfile(moxdir, 'roms_for_public', "B10K-BGCMIP"+sims(1)+suffix(isuf), 'Level2');
    
            sample = dir(fullfile(fileparts(filebase), 'Level1', '*temp.nc'));
            t = ncdateread(fullfile({sample.folder}, {sample.name}), 'ocean_time');

            Tmp = struct;
            Tmp.phyto       = readprocessed(fullfile(filebase, '*_average_phyto_integrated.nc'), 'phyto', Scs);
            Tmp.plfrac      = nan(1,1,nt);
            Tmp.zmicro      = readprocessed(fullfile(filebase, '*_average_microzoo_integrated.nc'), 'microzoo', Scs);
            Tmp.zmeso       = nan(1,1,nt);
            Tmp.npp         = readprocessed(fullfile(filebase, '*_diagnos_npp_integrated.nc'),  'npp', Scs);
            Tmp.zprod_micro = readprocessed(fullfile(filebase, '*_diagnos_gra_integrated.nc'),  'gra', Scs);
            Tmp.zprod_meso  = readprocessed(fullfile(filebase, '*_diagnos_pmor_integrated.nc'), 'pmor', Scs); % max bound: pmor
            Tmp.NH4s        = readprocessed(fullfile(filebase, '*_average_NH4_surface5m.nc'), 'NH4', Scs);
            Tmp.NH4b        = readprocessed(fullfile(filebase, '*_average_NH4_bottom5m.nc'),  'NH4', Scs);
            Tmp.NO3s        = readprocessed(fullfile(filebase, '*_average_NO3_surface5m.nc'), 'NO3', Scs);
            Tmp.NO3b        = readprocessed(fullfile(filebase, '*_average_NO3_bottom5m.nc'),  'NO3', Scs);
            
            Tmp = structfun(@(x) permute(x, [3 1 2]), Tmp, 'uni', 0);
            Metrics(ig,isuf).banas = table2timetable(struct2table(Tmp), 'rowtimes', t);

            save(outfile, 'Metrics');
        end
            
        % BESTNPZ (uses mg C m^-3 d^-1)
            
        if readflag(2,ig,isuf)
            disp("  "+sims(2));
            
            filebase = fullfile(moxdir, 'roms_for_public', "B10K-BGCMIP"+sims(2)+suffix(isuf), 'Level2');
            
            sample = dir(fullfile(fileparts(filebase), 'Level1', '*temp.nc'));
            t = ncdateread(fullfile({sample.folder}, {sample.name}), 'ocean_time');

            Tmp = struct;
            Tmp.phyto  = (readprocessed(fullfile(filebase, '*_average_PhS_integrated.nc'),  'PhS', Scs) + ...
                          readprocessed(fullfile(filebase, '*_average_PhL_integrated.nc'),  'PhL', Scs)) .* Param{2}.xi;
            Tmp.plfrac = (readprocessed(fullfile(filebase, '*_average_PhL_integrated.nc'),  'PhL', Scs).*Param{2}.xi)./Tmp.phyto;
            Tmp.zmicro =  readprocessed(fullfile(filebase, '*_average_MZL_integrated.nc'),  'MZL', Scs).*Param{2}.xi;
            Tmp.zmeso  = (readprocessed(fullfile(filebase, '*_average_Cop_integrated.nc'),  'Cop', Scs) + ...
                          readprocessed(fullfile(filebase, '*_average_NCaS_integrated.nc'), 'NCaS', Scs) + ...
                          readprocessed(fullfile(filebase, '*_average_NCaO_integrated.nc'), 'NCaO', Scs) + ...
                          readprocessed(fullfile(filebase, '*_average_EupS_integrated.nc'), 'EupS', Scs) + ...
                          readprocessed(fullfile(filebase, '*_average_EupO_integrated.nc'), 'EupO', Scs)).*Param{2}.xi;
            Tmp.npp =         (readprocessed(fullfile(filebase, '*_diagnos_prod_PhS_integrated.nc'),  'prod_PhS',  Scs) + ...
                               readprocessed(fullfile(filebase, '*_diagnos_prod_PhL_integrated.nc'),  'prod_PhL',  Scs)).*Param{2}.xi;
            Tmp.zprod_micro =  readprocessed(fullfile(filebase, '*_diagnos_prod_MZL_integrated.nc'),  'prod_MZL',  Scs).*Param{2}.xi;
            Tmp.zprod_meso  = (readprocessed(fullfile(filebase, '*_diagnos_prod_Cop_integrated.nc'),  'prod_Cop',  Scs) + ...
                               readprocessed(fullfile(filebase, '*_diagnos_prod_NCaS_integrated.nc'), 'prod_NCaS', Scs) + ...
                               readprocessed(fullfile(filebase, '*_diagnos_prod_NCaO_integrated.nc'), 'prod_NCaO', Scs) + ...
                               readprocessed(fullfile(filebase, '*_diagnos_prod_EupS_integrated.nc'), 'prod_EupS', Scs) + ...
                               readprocessed(fullfile(filebase, '*_diagnos_prod_EupO_integrated.nc'), 'prod_EupO', Scs)).*Param{2}.xi;
            Tmp.NH4s        = readprocessed(fullfile(filebase, '*_average_NH4_surface5m.nc'), 'NH4', Scs);
            Tmp.NH4b        = readprocessed(fullfile(filebase, '*_average_NH4_bottom5m.nc'),  'NH4', Scs);
            Tmp.NO3s        = readprocessed(fullfile(filebase, '*_average_NO3_surface5m.nc'), 'NO3', Scs);
            Tmp.NO3b        = readprocessed(fullfile(filebase, '*_average_NO3_bottom5m.nc'),  'NO3', Scs);
            
            Tmp = structfun(@(x) permute(x, [3 1 2]), Tmp, 'uni', 0);
            Metrics(ig,isuf).bestnpz = table2timetable(struct2table(Tmp), 'rowtimes', t);

            save(outfile, 'Metrics');
        end
            
        % COBALT debugging test (uses mol kg^-1 s^-1)
            
        if readflag(3,ig,isuf)

            disp("  "+sims(3));
            
            filebase = fullfile(moxdir, 'roms_for_public', "B10K-BGCMIP"+sims(3)+suffix(isuf), 'Level2');
            
            sample = dir(fullfile(fileparts(filebase), 'Level1', '*temp.nc'));
            t = ncdateread(fullfile({sample.folder}, {sample.name}), 'ocean_time');

            Tmp = struct;
            
            Tmp.phyto  = (readprocessed(fullfile(filebase, '*_average_nsm_integrated.nc'),  'nsm', Scs) + ...
                          readprocessed(fullfile(filebase, '*_average_nlg_integrated.nc'),  'nlg', Scs) + ...
                          readprocessed(fullfile(filebase, '*_average_ndi_integrated.nc'),  'ndi', Scs)).*1025.*1000;
            Tmp.plfrac =  readprocessed(fullfile(filebase, '*_average_nlg_integrated.nc'),  'nlg', Scs).*1025.*1000./Tmp.phyto;
            Tmp.zmicro =  readprocessed(fullfile(filebase, '*_average_nsmz_integrated.nc'), 'nsmz',Scs).*1025.*1000;
            Tmp.zmeso  = (readprocessed(fullfile(filebase, '*_average_nmdz_integrated.nc'), 'nmdz',Scs) + ...
                          readprocessed(fullfile(filebase, '*_average_nlgz_integrated.nc'), 'nlgz',Scs))*1025.*1000;
            Tmp.npp    = (readprocessed(fullfile(filebase, '*_diagnos_npp_sm_integrated.nc'), 'npp_sm',Scs) + ...
                          readprocessed(fullfile(filebase, '*_diagnos_npp_lg_integrated.nc'), 'npp_lg',Scs) + ...
                          readprocessed(fullfile(filebase, '*_diagnos_npp_di_integrated.nc'), 'npp_di',Scs)).*1025.*1000.*seconds(days(1));
            Tmp.zprod_micro = readprocessed(fullfile(filebase, '*_diagnos_zprod_sm_integrated.nc'), 'zprod_sm',Scs).*1025.*1000.*seconds(days(1));
            Tmp.zprod_meso = (readprocessed(fullfile(filebase, '*_diagnos_zprod_md_integrated.nc'), 'zprod_md',Scs) + ...
                              readprocessed(fullfile(filebase, '*_diagnos_zprod_lg_integrated.nc'), 'zprod_lg',Scs)).*1025.*1000.*seconds(days(1));
            Tmp.NH4s        = readprocessed(fullfile(filebase, '*_average_nh4_surface5m.nc'), 'nh4', Scs).*1025.*1000;
            Tmp.NH4b        = readprocessed(fullfile(filebase, '*_average_nh4_bottom5m.nc'),  'nh4', Scs).*1025.*1000;
            Tmp.NO3s        = readprocessed(fullfile(filebase, '*_average_no3_surface5m.nc'), 'no3', Scs).*1025.*1000;
            Tmp.NO3b        = readprocessed(fullfile(filebase, '*_average_no3_bottom5m.nc'),  'no3', Scs).*1025.*1000;
             
            Tmp = structfun(@(x) permute(x, [3 1 2]), Tmp, 'uni', 0);
            Metrics(ig,isuf).cobalt = table2timetable(struct2table(Tmp), 'rowtimes', t);

            save(outfile, 'Metrics');

        end

    end
end

%% Quick fix... eliminate 

% for isuf = [4 6]
%     for ig = nstation:-1:1
%         Metrics(ig,isuf).bestnpz = [];
%     end
% end

%% Subfunctions
%% ... readprocessed

function data = readprocessed(pth, var, Scs)

F = dir(pth);
Tmp = ncstruct(fullfile({F.folder},{F.name}), var, Scs);
data = Tmp.(var);

end

