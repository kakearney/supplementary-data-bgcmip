%bgcmip_setup

sims = ["banas", "bestnpz", "cobalt"];
nsim = length(sims);

simfolder = fullfile(mounteddir('klone'), 'GR011377_bgcmip', 'BGC_hindcasts_workdir');
% simfolder = fullfile(moxdir, 'kearney', 'BGC_hindcasts_workdir');

grdfile = fullfile(mounteddir('klone'), 'GR011377_bgcmip', 'ROMS_Datasets', 'grids', 'AlaskaGrids_Bering10K.nc');
% grdfile = fullfile(moxdir, 'ROMS_Datasets', 'grids', 'AlaskaGrids_Bering10K.nc');
Grd = ncstruct(grdfile);

nlayer = 30;

berapp = fullfile(mounteddir('klone'), 'GR011377_bgcmip', 'bering-Apps');
% berapp = fullfile(moxdir, 'kearney', 'bering-Apps');

pfile = fullfile(berapp, 'Apps', 'Bering_BGC_variants', "bering_bpar_"+sims+".yaml");
Param = cellfun(@(x) yaml.loadFile(x), pfile, 'uni', 0);

%% Locations to extract stuff

%% ... Mask: shelf-only

warning('off', 'MATLAB:polyshape:repairedBySimplify');

M = load(fullfile(mounteddir('klone'), 'GR011836_aclim','cmip6','analysis_kak', 'nbudget_masks.mat'));
% M = load(fullfile(moxdir, 'bering10k', 'output', 'forecasts', 'cmip6', 'analysis_kak', 'nbudget_masks.mat'));

mask = nan(size(Grd.h));
for ii = 1:size(M.masks,1)
    mask(M.masks{ii,2}) = ii;
    [ln, lt] = mask2poly(Grd.lon_psi, Grd.lat_psi, M.masks{ii,2}(2:end-1,2:end-1));
    pmask{ii} = polyshape(ln, lt);
end
pshelf = union([pmask{1:8}]);

shelfmask = reshape(pshelf.isinterior(Grd.lon_rho(:), Grd.lat_rho(:)), size(Grd.h)) & Grd.mask_rho==1;
[nxi, neta] = size(Grd.h);
nshelf = nnz(shelfmask);

%% ... Mask: transect plots (eliminate anything south of Aleutians)

pland = struct('lat', {cell(2,1)}, 'lon', {cell(2,1)});
[pland.lat{1}, pland.lon{1}] = borders('alaska');
[pland.lat{2}, pland.lon{2}] = borders('russia');
pland = polyshape(wrapTo360([pland.lon{:}]), [pland.lat{:}]);
ltbox = [50.5 56.76];
lnbox = [163 196];
aleut = intersect(polyshape(lnbox([1 1 2 2 1]), ltbox([1 2 2 1 1])), pland);
% aleut = sortregions(aleut, 'centroid', 'ascend');
tmp = regions(aleut);
[clon,clat] = centroid(tmp);
[~,imin] = min((clat - 56.594).^2 + (clon - 190.52).^2); % remove Pribs (in bounding box)
aleut = rmboundary(aleut, imin);

avert = cat(1, aleut.Vertices);
avert = unique(avert, 'rows');
avert = avert(~any(isnan(avert),2),:);
a = alphaShape(avert(:,1), avert(:,2));
a.Alpha = 10;
[~, apt] = boundaryFacets(a);
paleutstrip = polybuffer(polyshape(apt(:,1), apt(:,2)), 0.1);

pbasin = subtract(pmask{10}, polybuffer(polyshape(apt(:,1), apt(:,2)), 0.1));
pbasin = regions(sortregions(pbasin, 'area', 'descend'));
pbasin = pbasin(1);

pmain = union(pbasin, pshelf);

mask = double(~reshape(isinterior(paleutstrip, Grd.lon_rho(:), Grd.lat_rho(:)), size(Grd.h)));
% mask = double(reshape(isinterior(pmain, Grd.lon_rho(:), Grd.lat_rho(:)), size(Grd.h)));
mask(mask == 0) = NaN;

%% ... Stations

statfile = fullfile(mounteddir('klone'), 'GR011846_reem/kearney/runscripts/bering_stations_updated.csv');
% statfile = fullfile(moxdir, 'kearney/runscripts/bering_stations_updated.csv');
Station = readtable(statfile);

% Just EBS ones for now

onebs = pshelf.isinterior(wrapTo360(Station.lon), Station.lat);
Station = Station(onebs,:);

% Just a few select ones (ECCWO poster)

% samplestat = ["BS-2", "BS-8", "14", "2 low"];
% [~,loc] = ismember(samplestat, Station.name);
% Station = Station(loc,:);

[Station.xi, Station.eta] = ind2sub(size(Grd.h), Station.gridcell);