%% readnbudget

bgcmip_setup;

parforArg = Inf;

w = warning('off'); % suppress local.m warning about 1 grid cell

%% Budget masks

Budget(1).name = "M2single";
Budget(1).mask = false(size(Grd.h));

[ixi,ieta] = ind2sub([nxi neta], Station.gridcell(strcmp(Station.name, 'BS-2')));
Budget(1).mask(ixi,ieta) = true;

% SEBS: by shelf region

reg = ["SEBSinner", "SEBSmiddle", "SEBSouter"];
[~,rloc] = ismember(reg, M.masks(:,1));
for ir = 1:length(reg)
    Budget(ir+1).name = reg(ir);
    Budget(ir+1).mask = M.masks{rloc(ir),2};
end

% SEBS: total

Budget(5).name = "SEBS";
Budget(5).mask = any(cat(3, M.masks{rloc,2}),3);

%% Simulations to read

outpath = [fullfile("..", "bgcmip_nbudget_"+sims, "Out"); ...
           fullfile("..", "bgcmip_nbudget_"+sims+"bury", "Out"); ...
           fullfile("..", "bgcmip_loop_nbudget_"+sims, "Out"); ...
           fullfile("..", "bgcmip_loop_nbudget_"+sims+"bury", "Out")
           fullfile("..", "bgcmip_nbudget_"+sims+"burynoinfauna", "Out"); ...
           fullfile("..", "bgcmip_loop_nbudget_"+sims+"burynoinfauna", "Out")];
outtype = repmat(sims, length(outpath), 1);

keepmask = true(size(outtype));
keepmask(2,1) = false; % No bury variant for banas
keepmask(4,1) = false;
keepmask(2,3) = false; % Cobalt bury variant not yet useful
keepmask(4,3) = false;
keepmask(5:6,1) = false; % burynoinfauna only for bestnpz
keepmask(5:6,3) = false; 

outpath = outpath(keepmask);
outtype = outtype(keepmask);

nout = length(outtype);

%% ... Read data, store values in graph tables

for ib = 5:5 %1:length(Budget)

    G = cell(nout,1);
    tg = cell(nout,1);

    budgetfile = "nbudgets_all_"+Budget(ib).name+".mat";
    if exist(budgetfile, 'file')
        Tmp = load(budgetfile);
        [tf, loc] = ismember(Tmp.outpath, outpath);
        G(loc) = Tmp.G;
        tg(loc) = Tmp.tg;
    end
        
    % Subgrid stuff (didn't have anything other than rho-variables, so we 
    % won't worry about other coordinates)

    xrlim = [find(any(Budget(ib).mask,2),1,'first') find(any(Budget(ib).mask,2),1,'last')]; 
    erlim = [find(any(Budget(ib).mask,1),1,'first') find(any(Budget(ib).mask,1),1,'last')]; 

    Scs = struct('xi_rho',  [xrlim(1) diff(xrlim)+1 1], ...
                 'eta_rho', [erlim(1) diff(erlim)+1 1]);

    GrdSub = ncstruct(grdfile, Scs, 'h', 'area_feast');
    [nxisub, netasub] = size(GrdSub.h);

    masksub = Budget(ib).mask(xrlim(1):xrlim(2), erlim(1):erlim(2));

    for ii = 1:nout
        if isempty(G{ii})
            G{ii} = bgcfluxgraph(outtype{ii});

            nnode = numnodes(G{ii});
            nedge = numedges(G{ii});
    
            % Read state variable depth-integrated concentration and flux variable
            % depth-integrated flux concentration
            % This step assumes that averages and diagnostics files are archived
            % with the same time step, which is necessary to match up the zeta
            % values
    
            afiles = dir(fullfile(outpath{ii}, '*avg*'));
            afiles = fullfile({afiles.folder}, {afiles.name});
    
            dfiles = dir(fullfile(outpath{ii}, '*dia*'));
            dfiles = fullfile({dfiles.folder}, {dfiles.name});
    
            ta = cellfun(@(x) ncdateread(x, 'ocean_time'), afiles, 'uni', 0);
            td = cellfun(@(x) ncdateread(x, 'ocean_time'), dfiles, 'uni', 0);
    
            if ~isequal(ta, td)
                warning('Mismatched avg vs dia times: %s\n', sims{ii});
            end
    
            nfile = length(dfiles);
    
            t = cat(1, td{1:nfile});
    
            nodeval = cell(nnode, nfile);
            edgeval = cell(nedge, nfile);
    
            parfor (ifl = 1:nfile, parforArg)
                fprintf('Reading data, sim %d/%d, file %d/%d\n', ii, nsim, ifl, nfile);
    
                % Read all variables
                A = ncstruct(afiles{ifl}, Scs);
                D = ncstruct(dfiles{ifl}, Scs);
    
                % Depth of layers
                [~,zw] = calcromsz(GrdSub.h, A.zeta, nlayer); 
                dz = diff(zw, 1, 3);
                nt = size(A.zeta, 3);
    
                % Permute all rho-variables to xi x eta x depth x time
                fld = fieldnames(A);
                for id = 1:length(fld)
                    if isequal(size(A.(fld{id})), [nxisub netasub nt])
                        A.(fld{id}) = permute(A.(fld{id}), [1 2 4 3]);
                    end
                end
    
                % Values per grid cell
                for in = 1:nnode
                    nodeval{in,ifl} = G{ii}.Nodes.fun{in}(A, dz);
                    % Depth sum, spatial mean
                    nodeval{in,ifl} = local(permute(sum(nodeval{in,ifl},3), [1 2 4 3]), masksub, 'weight', GrdSub.area_feast, 'omitnan');
                end
                for ie = 1:nedge
                    edgeval{ie,ifl} = G{ii}.Edges.fun{ie}(D, dz);
                    % Depth sum, spatial mean
                    edgeval{ie,ifl} = local(permute(sum(edgeval{ie,ifl},3), [1 2 4 3]), masksub, 'weight', GrdSub.area_feast, 'omitnan');
                end
    
            end
    
            for in = 1:nnode
                G{ii}.Nodes.val{in} = cat(1, nodeval{in,:});
            end
            for ie = 1:nedge
                G{ii}.Edges.val{ie} = cat(1, edgeval{ie,:});
            end
    
            % Fix any repeat times (keep first)
    
            [tunq, iunq] = unique(t, 'first');
    
            if length(tunq) < length(t)
                G{ii}.Nodes.val = cellfun(@(x) x(iunq), G{ii}.Nodes.val, 'uni', 0);
                G{ii}.Edges.val = cellfun(@(x) x(iunq), G{ii}.Edges.val, 'uni', 0);
            end
    
            tg{ii} = tunq;
        end
    end

    save("nbudgets_all_"+Budget(ib).name, 'G', 'tg', 'outpath');
end
 
%%

warning(w);

% end