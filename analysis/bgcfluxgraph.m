function G = bgcfluxgraph(bgcname, varargin)

p = inputParser;
p.addParameter('bappfol', fullfile(moxdir, 'kearney', 'bering-Apps'));
p.addParameter('param', [])
p.parse(varargin{:});

Opt = p.Results;

switch bgcname
    case 'banas'
        pfile = fullfile(Opt.bappfol, "subApps", "BIO_BANAS", "banas_bpar.yaml");
        iofile = fullfile(Opt.bappfol, "subApps", "BIO_BANAS", "banas_io.yaml");

    case 'bestnpz'
        pfile = fullfile(Opt.bappfol, "subApps", "BEST_NPZ", "bestnpz_bpar.yaml");
        iofile = fullfile(Opt.bappfol, "subApps", "BEST_NPZ", "bestnpz_io.yaml");

    case 'cobalt'
        pfile = fullfile(Opt.bappfol, "subApps", "BIO_COBALT", "cobalt_bpar.yaml");
        iofile = fullfile(Opt.bappfol, "subApps", "BIO_COBALT", "cobalt_io.yaml");
end

% bappfol = "../../bering-Apps/";
% sims = ["banas", "bestnpz", "cobalt"];
% nsim = length(sims);

% P = load('bgcmip_params.mat');
% P = P.Param;

P = yaml.loadFile(pfile);
Io = yaml.loadFile(iofile);
Io = struct2table(cat(1, Io.metadata{:}));

% 
% iofiles = [...
%     fullfile(bappfol, "subApps", "BIO_BANAS", "banas_io.yaml")
%     fullfile(bappfol, "subApps", "BEST_NPZ", "bestnpz_io.yaml")
%     fullfile(bappfol, "subApps", "BIO_COBALT", "cobalt_io.yaml")];
% 
% Io = cell(nsim,1);
% for ii = 1:nsim
%     Tmp = yaml.loadFile(iofiles{ii});
%     Io{ii} = struct2table(cat(1, Tmp.metadata{:}));
% end

sperday = seconds(days(1));

switch bgcname

    case 'banas'

        %----------
        % Banas
        %----------
        
        ntype = {...
            "nh4"       "N"     @(A,dz) dz.*A.NH4
            "no3"       "N"     @(A,dz) dz.*A.NO3
            "phyto"     "P"     @(A,dz) dz.*A.phyto
            "microzoo"  "Z"     @(A,dz) dz.*A.microzoo
            "det_small" "D"     @(A,dz) dz.*A.det_small
            "det_large" "D"     @(A,dz) dz.*A.det_large
            "out"       "X"     @(A,dz) zeros(size(A.NH4))
            };   
        ntbl = table(string(ntype(:,1)), string(ntype(:,2)), ntype(:,3), 'variablenames', {'Name', 'type', 'fun'});
        
        tmp = {...
            "npp"   "nh4"       "phyto"         @(D,dz)          dz.*D.npp.*(1-D.fratio)
            "npp"   "no3"       "phyto"         @(D,dz)          dz.*D.fratio.*D.npp
            "gra"   "phyto"     "microzoo"      @(D,dz)          dz.*D.gra
            "ege"   "microzoo"  "det_small"     @(D,dz)          dz.*D.gra.*(1-P.epsil-P.fex)
            "ege"   "microzoo"  "nh4"           @(D,dz)          dz.*D.gra.*P.fex
            "mor"   "phyto"     "det_small"     @(D,dz)          dz.*D.pmor
            "agg"   "phyto"     "det_large"     @(D,dz)          dz.*D.agg
            "mor"   "microzoo"  "out"           @(D,dz)          dz.*D.zmor
            "rem"   "det_small" "nh4"           @(D,dz)          dz.*D.srem
            "rem"   "det_large" "nh4"           @(D,dz)          dz.*D.lrem
            "nit"   "nh4"       "no3"           @(D,dz)          dz.*D.nit
            "hadv"  "out"       "nh4"           @(D,dz) sperday.*dz.*D.NH4_hadv  
            "hadv"  "out"       "no3"           @(D,dz) sperday.*dz.*D.NO3_hadv 
            "hadv"  "out"       "phyto"         @(D,dz) sperday.*dz.*D.phyto_hadv 
            "hadv"  "out"       "microzoo"      @(D,dz) sperday.*dz.*D.microzoo_hadv
            "hadv"  "out"       "det_small"     @(D,dz) sperday.*dz.*D.det_small_hadv
            "hadv"  "out"       "det_large"     @(D,dz) sperday.*dz.*D.det_large_hadv
            "hdif"  "out"       "nh4"           @(D,dz) sperday.*dz.*D.NH4_hdiff  
            "hdif"  "out"       "no3"           @(D,dz) sperday.*dz.*D.NO3_hdiff 
            "hdif"  "out"       "phyto"         @(D,dz) sperday.*dz.*D.phyto_hdiff 
            "hdif"  "out"       "microzoo"      @(D,dz) sperday.*dz.*D.microzoo_hdiff
            "hdif"  "out"       "det_small"     @(D,dz) sperday.*dz.*D.det_small_hdiff
            "hdif"  "out"       "det_large"     @(D,dz) sperday.*dz.*D.det_large_hdiff
            "vadv"  "out"       "nh4"           @(D,dz) sperday.*dz.*D.NH4_vadv  
            "vadv"  "out"       "no3"           @(D,dz) sperday.*dz.*D.NO3_vadv 
            "vadv"  "out"       "phyto"         @(D,dz) sperday.*dz.*D.phyto_vadv 
            "vadv"  "out"       "microzoo"      @(D,dz) sperday.*dz.*D.microzoo_vadv
            "vadv"  "out"       "det_small"     @(D,dz) sperday.*dz.*D.det_small_vadv
            "vadv"  "out"       "det_large"     @(D,dz) sperday.*dz.*D.det_large_vadv
            "vdif"  "out"       "nh4"           @(D,dz) sperday.*dz.*D.NH4_vdiff  
            "vdif"  "out"       "no3"           @(D,dz) sperday.*dz.*D.NO3_vdiff 
            "vdif"  "out"       "phyto"         @(D,dz) sperday.*dz.*D.phyto_vdiff 
            "vdif"  "out"       "microzoo"      @(D,dz) sperday.*dz.*D.microzoo_vdiff
            "vdif"  "out"       "det_small"     @(D,dz) sperday.*dz.*D.det_small_vdiff
            "vdif"  "out"       "det_large"     @(D,dz) sperday.*dz.*D.det_large_vdiff
            "snk"   "out"       "det_large"     @(D,dz)          dz.*D.sinkl
            "snk"   "out"       "det_small"     @(D,dz)          dz.*D.sinks
            };
            
        etbl = table(string(tmp(:,2:3)), string(tmp(:,1)), tmp(:,4), 'variablenames', {'EndNodes', 'type', 'fun'});
        
        G = digraph(etbl, ntbl);

    case 'bestnpz'
        %----------
        % BEST_NPZ
        %----------
        
        % Parse edges from flux diagnostics
        
        isflx = startsWith(Io.index_code, "iDbio3(iflx_");
        flx = Io.variable(isflx);
        flx = setdiff(flx, ["ice_Alk_flux", ...
                            "ice_TIC_flux", ...
                            "advdiff_TIC", ...
                            "advdiff_alkalinity", ...
                            "advdiff_oxygen"]);
        flx = flx(~startsWith(flx, 'advdiff'));
        
        tmp = regexp(strrep(strrep(strrep(flx, "INO3", "IceNO3"), "INH4", "IceNH4"), "IPhL", "IcePhL"), '_', 'split');
        tmp = tmp(cellfun(@length, tmp) == 3);
        tmp = cat(1, tmp{:});
        
        % State variables
        
        ntype = {...
            "NH4"    "N"    @(A,dz) dz    .*      A.NH4    
            "NO3"    "N"    @(A,dz) dz    .*      A.NO3
            "IceNH4" "N"    @(A,dz) P.aidz.*      A.IceNH4
            "IceNO3" "N"    @(A,dz) P.aidz.*      A.IceNO3
            "PhS"    "P"    @(A,dz) dz    .*P.xi.*A.PhS
            "PhL"    "P"    @(A,dz) dz    .*P.xi.*A.PhL
            "IcePhL" "P"    @(A,dz) P.aidz.*P.xi.*A.IcePhL
            "MZL"    "Z"    @(A,dz) dz    .*P.xi.*A.MZL
            "Cop"    "Z"    @(A,dz) dz    .*P.xi.*A.Cop
            "NCaS"   "Z"    @(A,dz) dz    .*P.xi.*A.NCaS
            "NCaO"   "Z"    @(A,dz) dz    .*P.xi.*A.NCaO
            "EupS"   "Z"    @(A,dz) dz    .*P.xi.*A.EupS
            "EupO"   "Z"    @(A,dz) dz    .*P.xi.*A.EupO
            "Jel"    "Z"    @(A,dz) dz    .*P.xi.*A.Jel
            "Ben"    "Z"    @(A,dz)         P.xi.*A.Ben
            "Det"    "D"    @(A,dz) dz    .*P.xi.*A.Det
            "DetF"   "D"    @(A,dz) dz    .*P.xi.*A.DetF
            "DetBen" "D"    @(A,dz)         P.xi.*A.DetBen
            "Out"    "X"    @(A,dz) zeros(size(A.NH4))};
        ntbl = table(string(ntype(:,1)), string(ntype(:,2)), ntype(:,3), 'variablenames', {'Name', 'type', 'fun'});
        
        % Add adv/dif fluxes for each pelagic variables
        
        tvar = setdiff(string(ntype(:,1)), ["Ben", "DetBen", "IceNH4", "IceNO3", "IcePhL", "Out"]);
        ntv = length(tvar);
        tmp = [tmp; ...
            [[repmat("hadv",ntv,1); repmat("vadv",ntv,1); repmat("hdif",ntv,1); repmat("vdif",ntv,1)], ...
            repmat("Out",ntv*4,1) repmat(tvar,4,1)]];
        flx = [flx; tvar+"_hadv"; tvar+"_vadv"; tvar+"_hdiff"; tvar+"_vdiff"];
        
        etbl = table(tmp(:,2:3), tmp(:,1), 'variablenames', {'EndNodes', 'type'});
        
        for ie = 1:height(etbl)
            if contains(etbl.type{ie}, {'adv', 'dif'})
                % Adv/dif are in mmolN/m3/s or mgC/m3/s
                if ismember(etbl.EndNodes{ie,2}, {'NO3', 'NH4'})
                    etbl.fun{ie} = @(D,dz) sperday.*dz.*D.(flx{ie});
                else
                    etbl.fun{ie} = @(D,dz) sperday.*dz.*D.(flx{ie}).*P.xi;
                end
            else
                % Others are in mgC/m2/d (even NO3 and NH4)
                etbl.fun{ie} = @(D,dz) D.(flx{ie}).*P.xi;
            end
        end
        
        % Build graph
        
        G = digraph(etbl, ntbl);

    case 'cobalt'

        %----------
        % COBALT
        %----------
        
        % Parse edges from flux diagnostics
        
        isflx = startsWith(Io.index_code, "iDbio3(iflx_");
        flx = Io.variable(isflx);
        
        tmp = regexp(strrep(flx, "bac", "nbact"), '_', 'split');
        tmp = tmp(cellfun(@length, tmp) == 3);
        tmp = cat(1, tmp{:});
        
        % State variables
        
        mll2mlr = 1000 * 1025; % molal to molar, i.e. mol/kg -> mmol/m3
        
        ntype = {...
            "nh4"   "N" @(A,dz) mll2mlr.*dz.*A.nh4
            "no3"   "N" @(A,dz) mll2mlr.*dz.*A.no3
            "n2"    "N" @(A,dz) zeros(size(A.nh4))
            "ndi"   "P" @(A,dz) mll2mlr.*dz.*A.ndi
            "nsm"   "P" @(A,dz) mll2mlr.*dz.*A.nsm
            "nlg"   "P" @(A,dz) mll2mlr.*dz.*A.nlg
            "nsmz"  "Z" @(A,dz) mll2mlr.*dz.*A.nsmz
            "nmdz"  "Z" @(A,dz) mll2mlr.*dz.*A.nmdz
            "nlgz"  "Z" @(A,dz) mll2mlr.*dz.*A.nlgz
            "hip"   "Z" @(A,dz) zeros(size(A.nh4))
            "ndet"  "D" @(A,dz) mll2mlr.*dz.*A.ndet
            "ldon"  "D" @(A,dz) mll2mlr.*dz.*A.ldon
            "sldon" "D" @(A,dz) mll2mlr.*dz.*A.sldon
            "srdon" "D" @(A,dz) mll2mlr.*dz.*A.srdon
            "nbact" "B" @(A,dz) mll2mlr.*dz.*A.nbact
            "out"   "X" @(A,dz) zeros(size(A.nh4))
            };
        
        ntbl = table(string(ntype(:,1)), string(ntype(:,2)), ntype(:,3), 'variablenames', {'Name', 'type', 'fun'});
        
        % Add adv/dif fluxes for each pelagic variable
        
        tvar = setdiff(string(ntype(:,1)), ["hip", "out", "n2"]);
        ntv = length(tvar);
        tmp = [tmp; ...
            [[repmat("hadv",ntv,1); repmat("vadv",ntv,1); repmat("hdif",ntv,1); repmat("vdif",ntv,1)], ...
            repmat("out",ntv*4,1) repmat(tvar,4,1)]];
        flx = [flx; tvar+"_hadv"; tvar+"_vadv"; tvar+"_hdiff"; tvar+"_vdiff"];
        
        
        etbl = table(tmp(:,2:3), tmp(:,1), 'variablenames', {'EndNodes', 'type'});
        
        for ie = 1:height(etbl)
            etbl.fun{ie} = @(D,dz) D.(flx{ie}).*sperday.*mll2mlr.*dz;
        end
        
        G = digraph(etbl, ntbl);
end
