% Create initialization file for second loop of simulations

simroms  = ["BIO_BANAS", "BEST_NPZ", "BIO_COBALT"];
sims     = ["banas",     "bestnpz",  "cobalt"];
nsim = length(sims);

datafol = fullfile(moxdir, "ROMS_Datasets");

%%

dryrun = 'run';

for ii = 1:nsim
    rst = dir(fullfile('..', "bgcmip_"+sims{ii}, 'Out', '*rst*.nc'));

    oldrst = fullfile(rst(end).folder, rst(end).name);
    newini = fullfile('..', "ini_hindcastloop2_"+simroms{ii}+".nc");

    t = ncdateread(oldrst, 'ocean_time');
    [~,imax] = max(t);

    cmd = sprintf('ncks -F -d ocean_time,%d,%d %s %s', imax, imax, oldrst, newini);
    safesystem(cmd, dryrun);

    tunit = ncreadatt(oldrst, 'ocean_time', 'units');
    tnew = cftime(datetime(1990,1,15), tunit, [], 'reverse');
    if strcmp(dryrun, 'run')
        ncwrite(newini, 'ocean_time', tnew);
    end

end

%% Compare old to new

inifile = fullfile(datafol, 'initial', "ini_hindcast_unnested_Bering10K_" + simroms + ".nc");
newfile = fullfile('..', "ini_hindcastloop2_" + simroms + ".nc");

I1 = ncinfo(inifile{1});
I2 = ncinfo(newfile{1});


%% Check bry

% How did my bry files lose their fill values?  Ugh.

B = dir(fullfile(datafol, 'CFS', '*', '*bryocn*.nc'));
B = B(~contains({B.folder}, 'old') & ~contains({B.name}, 'old'));

cmd = {};
for ib = 1:length(B)
    bname =  fullfile(B(ib).folder, B(ib).name);
    cmd = [...
        cmd
        sprintf('echo "%s"', B(ib).name)
        sprintf('cp %s %s', bname, strrep(bname, '.nc', '.old2.nc'));
        sprintf('cdo setmissval,-999999 %s test.nc', bname)
        sprintf('mv test.nc %s', bname)
        ''
        ];
end

cmd = strrep(cmd, GetFullPath(datafol)+"/", '');

fid = fopen(fullfile(datafol, 'fixbrynans.sh'), 'wt');
fprintf(fid, '%s\n', cmd{:});
fclose(fid);




% testfile = fullfile(moxdir, 'ROMS_Datasets/CFS/1990/CFS-ocean-Bering10K-N30-bryocn-1990.nc');