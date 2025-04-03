function error = func_chemo_error(WT_y_steady_clim,   WT_met_reac_steady_clim, ...
                                  WT_y_steady_cnlim1, WT_met_reac_steady_cnlim1, ...
                                  WT_y_steady_cnlim2, WT_met_reac_steady_cnlim2, ...
                                  WT_y_steady_nlim,   WT_met_reac_steady_nlim, ...
                                  num_y, num_flux, num_prot, ...
                                  pt_uM_to_gPerL, lp_uM_to_gPerL, lipid_sc_c,...
                                  D_chem, par)
%
% This function calculate the error between experiment and simulation for chemostat
% Relative Error = (abs(exp data/sim data - 1))^2
% Absolute Error = (abs(exp data - sim data))^2
%
%

%% ----------- load exp data ---------------------------------------------
[D_clim, cell_clim, glucose_clim, Jgy_clim, Jeh_clim,... 
precursor_clim, atp_clim, aa_clim,  prot_sec_clim, ...
lipid_clim, protein_clim,carbo_clim, RNA_clim, paperinfo_clim,...
color_clim, shape_clim, gas_clim] = SC_C_lim_chem();   

[D_nlim, cell_nlim, glucose_nlim, Jgy_nlim, Jeh_nlim, nh4_nlim, gas_nlim, ... 
precursor_nlim, atp_nlim, aa_nlim,  prot_sec_nlim, ...
lipid_nlim, protein_nlim,carbo_nlim, RNA_nlim, paperinfo_nlim,...                    
color_nlim, shape_nlim, Jnh4_nlim]  = SC_N_lim_chem(); 

lipid_constant_part = lipid_sc_c; 

%% -------- load simu data and processing --------------------------------
%%% proteins frac (sim)
total_protein_con_clim = sum(WT_y_steady_clim(:,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);
total_protein_con_nlim = sum(WT_y_steady_nlim(:,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);
total_protein_con_cnlim1 = sum(WT_y_steady_cnlim1(:,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);
total_protein_con_cnlim2 = sum(WT_y_steady_cnlim2(:,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);

WT_y_steady_clim(:,[num_y.r: num_y.lo]) = (WT_y_steady_clim(:,[num_y.r: num_y.lo]).*par.l')./ total_protein_con_clim;
WT_y_steady_nlim(:,[num_y.r: num_y.lo]) = (WT_y_steady_nlim(:,[num_y.r: num_y.lo]).*par.l')./ total_protein_con_nlim;
WT_y_steady_cnlim1(:,[num_y.r: num_y.lo]) = (WT_y_steady_cnlim1(:,[num_y.r: num_y.lo]).*par.l')./ total_protein_con_cnlim1;
WT_y_steady_cnlim2(:,[num_y.r: num_y.lo]) = (WT_y_steady_cnlim2(:,[num_y.r: num_y.lo]).*par.l')./ total_protein_con_cnlim2;

WT_y_steady_clim_frac = WT_y_steady_clim;
WT_y_steady_nlim_frac = WT_y_steady_nlim;
WT_y_steady_cnlim1_frac = WT_y_steady_cnlim1;
WT_y_steady_cnlim2_frac = WT_y_steady_cnlim2;


%%% REZ fraction calculation before normalizing
% fraction
R_chemo_sim_clim = WT_y_steady_clim(:,num_y.r).*par.l(num_prot.r)' ./ total_protein_con_clim;
Z_chemo_sim_clim = WT_y_steady_clim(:,num_y.z).*par.l(num_prot.z)' ./ total_protein_con_clim;
E_chemo_sim_clim = sum(WT_y_steady_clim(:,num_y.gy:num_y.lo).*par.l(num_prot.gy:num_prot.lo)',2) ./ total_protein_con_clim;

R_chemo_sim_nlim = WT_y_steady_nlim(:,num_y.r).*par.l(num_prot.r)' ;
Z_chemo_sim_nlim = WT_y_steady_nlim(:,num_y.z).*par.l(num_prot.z)' ;
E_chemo_sim_nlim = sum(WT_y_steady_nlim(:,num_y.gy:num_y.lo), 2);

% relative fold change
R_chemo_sim_fc = R_chemo_sim_clim./R_chemo_sim_clim(1); 
Z_chemo_sim_fc = Z_chemo_sim_clim./Z_chemo_sim_clim(1);
E_chemo_sim_fc = E_chemo_sim_clim./E_chemo_sim_clim(1);

R_chemo_sim_fc_nlim = R_chemo_sim_nlim./R_chemo_sim_nlim(1); 
Z_chemo_sim_fc_nlim = Z_chemo_sim_nlim./Z_chemo_sim_nlim(1);
E_chemo_sim_fc_nlim = E_chemo_sim_nlim./E_chemo_sim_nlim(1);


%%% normalize precursor, amino acid, atp, protein with respect to first data point 
D_exp = 0.05; % first D in exp
[norm_D, norm_indx] = min(abs(D_chem - D_exp)); 
% aa atp (sim)
WT_y_steady_clim(:,[num_y.aa_in, num_y.ae]) = WT_y_steady_clim(:,[num_y.aa_in, num_y.ae]) ./ WT_y_steady_clim(norm_indx,[num_y.aa_in, num_y.ae]);
WT_y_steady_nlim(:,[num_y.aa_in, num_y.ae]) = WT_y_steady_nlim(:,[num_y.aa_in, num_y.ae]) ./ WT_y_steady_nlim(norm_indx,[num_y.aa_in, num_y.ae]);
WT_y_steady_cnlim1(:,[num_y.aa_in, num_y.ae]) = WT_y_steady_cnlim1(:,[num_y.aa_in, num_y.ae]) ./ WT_y_steady_cnlim1(norm_indx,[num_y.aa_in, num_y.ae]);
WT_y_steady_cnlim2(:,[num_y.aa_in, num_y.ae]) = WT_y_steady_cnlim2(:,[num_y.aa_in, num_y.ae]) ./ WT_y_steady_cnlim2(norm_indx,[num_y.aa_in, num_y.ae]);

% precursor (sim)
%norm_indx = numel(WT_y_steady_clim(:,[num_y.pc])); 
WT_y_steady_clim(:,[num_y.pc]) = WT_y_steady_clim(:,[num_y.pc]) ./ WT_y_steady_clim(norm_indx,[num_y.pc]);
WT_y_steady_nlim(:,[num_y.pc]) = WT_y_steady_nlim(:,[num_y.pc]) ./ WT_y_steady_nlim(norm_indx,[num_y.pc]);
WT_y_steady_cnlim1(:,[num_y.pc]) = WT_y_steady_cnlim1(:,[num_y.pc]) ./ WT_y_steady_cnlim1(norm_indx,[num_y.pc]);
WT_y_steady_cnlim2(:,[num_y.pc]) = WT_y_steady_cnlim2(:,[num_y.pc]) ./ WT_y_steady_cnlim2(norm_indx,[num_y.pc]);

% protein log2fc 
WT_y_steady_clim(:,[num_y.r: num_y.lo]) = log2(WT_y_steady_clim(:,[num_y.r: num_y.lo])./WT_y_steady_clim(norm_indx,[num_y.r: num_y.lo])); 
WT_y_steady_nlim(:,[num_y.r: num_y.lo]) = log2(WT_y_steady_nlim(:,[num_y.r: num_y.lo])./WT_y_steady_nlim(norm_indx,[num_y.r: num_y.lo])); 
WT_y_steady_cnlim1(:,[num_y.r: num_y.lo]) = log2(WT_y_steady_cnlim1(:,[num_y.r: num_y.lo])./WT_y_steady_cnlim1(norm_indx,[num_y.r: num_y.lo])); 
WT_y_steady_cnlim2(:,[num_y.r: num_y.lo]) = log2(WT_y_steady_cnlim2(:,[num_y.r: num_y.lo])./WT_y_steady_cnlim2(norm_indx,[num_y.r: num_y.lo])); 

%% Organize N-conditions
% N-conditions
cnlim1 = {'Yu30', 'Yu50', 'Yu2021'};
cnlim2 = {'Yu115', 'kumar_ln'};
nlim = {'hackett','boer'};

fields_to_copy = ismember(fieldnames(D_nlim), cnlim1); % Logical array
all_fields = fieldnames(D_nlim); % Get all fields of the original structure
for i = 1:length(all_fields)
    if fields_to_copy(i)
        field_name = all_fields{i};
        if isfield(D_nlim, field_name)
            D_cnlim1.(field_name) = D_nlim.(field_name);
        end
        if isfield(cell_nlim, field_name)
            cell_cnlim1.(field_name) = cell_nlim.(field_name);
        end
        if isfield(Jgy_nlim, field_name)
            Jgy_cnlim1.(field_name) = Jgy_nlim.(field_name);
        end
        if isfield(Jeh_nlim, field_name)
            Jeh_cnlim1.(field_name) = Jeh_nlim.(field_name);
        end
        if isfield(nh4_nlim, field_name)
            nh4_cnlim1.(field_name) = nh4_nlim.(field_name);
        end
        precs = fieldnames(precursor_nlim);    
        for k = 1: length(precs)
            prec = precs{k};
            if isfield(precursor_nlim.(prec), field_name)
                precursor_cnlim1.(prec).(field_name) = precursor_nlim.(prec).(field_name);
            end   
        end
        aas = fieldnames(aa_nlim);    
        for k = 1: length(aas)
            aa = aas{k};
            if isfield(aa_nlim.(aa), field_name)
                aa_cnlim1.(aa).(field_name) = aa_nlim.(aa).(field_name);
            end   
        end
        if isfield(prot_sec_nlim, field_name)
            prot_sec_cnlim1.(field_name) = prot_sec_nlim.(field_name);
        end
        if isfield(lipid_nlim, field_name)
            lipid_cnlim1.(field_name) = lipid_nlim.(field_name);
        end
        if isfield(protein_nlim, field_name)
            protein_cnlim1.(field_name) = protein_nlim.(field_name);
        end
        proteins = fieldnames(prot_sec_nlim);    
        for k = 1: length(proteins)
            pro = proteins{k};
            if isfield(prot_sec_nlim.(pro), field_name)
                prot_sec_cnlim1.(pro).(field_name) = prot_sec_nlim.(pro).(field_name);
            end   
        end
    end
end

fields_to_copy = ismember(fieldnames(D_nlim), cnlim2); % Logical array
all_fields = fieldnames(D_nlim); % Get all fields of the original structure
for i = 1:length(all_fields)
    if fields_to_copy(i)
        field_name = all_fields{i};
        if isfield(D_nlim, field_name)
            D_cnlim2.(field_name) = D_nlim.(field_name);
        end
        if isfield(cell_nlim, field_name)
            cell_cnlim2.(field_name) = cell_nlim.(field_name);
        end
        if isfield(Jgy_nlim, field_name)
            Jgy_cnlim2.(field_name) = Jgy_nlim.(field_name);
        end
        if isfield(Jeh_nlim, field_name)
            Jeh_cnlim2.(field_name) = Jeh_nlim.(field_name);
        end
        if isfield(nh4_nlim, field_name)
            nh4_cnlim2.(field_name) = nh4_nlim.(field_name);
        end
        precs = fieldnames(precursor_nlim);    
        for k = 1: length(precs)
            prec = precs{k};
            if isfield(precursor_nlim.(prec), field_name)
                precursor_cnlim2.(prec).(field_name) = precursor_nlim.(prec).(field_name);
            end   
        end
        aas = fieldnames(aa_nlim);    
        for k = 1: length(aas)
            aa = aas{k};
            if isfield(aa_nlim.(aa), field_name)
                aa_cnlim2.(aa).(field_name) = aa_nlim.(aa).(field_name);
            end   
        end
        if isfield(prot_sec_nlim, field_name)
            prot_sec_cnlim2.(field_name) = prot_sec_nlim.(field_name);
        end
        if isfield(lipid_nlim, field_name)
            lipid_cnlim2.(field_name) = lipid_nlim.(field_name);
        end
        if isfield(protein_nlim, field_name)
            protein_cnlim2.(field_name) = protein_nlim.(field_name);
        end
        proteins = fieldnames(prot_sec_nlim);    
        for k = 1: length(proteins)
            pro = proteins{k};
            if isfield(prot_sec_nlim.(pro), field_name)
                prot_sec_cnlim2.(pro).(field_name) = prot_sec_nlim.(pro).(field_name);
            end   
        end
    end
end

fields_to_copy = ismember(fieldnames(D_nlim), nlim); % Logical array
all_fields = fieldnames(D_nlim); % Get all fields of the original structure
for i = 1:length(all_fields)
    if fields_to_copy(i)
        field_name = all_fields{i};
        if isfield(D_nlim, field_name)
            D_nlimf.(field_name) = D_nlim.(field_name);
        end
        if isfield(cell_nlim, field_name)
            cell_nlimf.(field_name) = cell_nlim.(field_name);
        end
        if isfield(Jgy_nlim, field_name)
            Jgy_nlimf.(field_name) = Jgy_nlim.(field_name);
        end
        if isfield(Jeh_nlim, field_name)
            Jeh_nlimf.(field_name) = Jeh_nlim.(field_name);
        end
        if isfield(nh4_nlim, field_name)
            nh4_nlimf.(field_name) = nh4_nlim.(field_name);
        end
        precs = fieldnames(precursor_nlim);    
        for k = 1: length(precs)
            prec = precs{k};
            if isfield(precursor_nlim.(prec), field_name)
                precursor_nlimf.(prec).(field_name) = precursor_nlim.(prec).(field_name);
            end   
        end
        aas = fieldnames(aa_nlim);    
        for k = 1: length(aas)
            aa = aas{k};
            if isfield(aa_nlim.(aa), field_name)
                aa_nlimf.(aa).(field_name) = aa_nlim.(aa).(field_name);
            end   
        end
        if isfield(prot_sec_nlim, field_name)
            prot_sec_nlimf.(field_name) = prot_sec_nlim.(field_name);
        end
        if isfield(lipid_nlim, field_name)
            lipid_nlimf.(field_name) = lipid_nlim.(field_name);
        end
        if isfield(protein_nlim, field_name)
            protein_nlimf.(field_name) = protein_nlim.(field_name);
        end
        proteins = fieldnames(prot_sec_nlim);    
        for k = 1: length(proteins)
            pro = proteins{k};
            if isfield(prot_sec_nlim.(pro), field_name)
                prot_sec_nlimf.(pro).(field_name) = prot_sec_nlim.(pro).(field_name);
            end   
        end
    end
end

%% ------------------ pair up C-lim simu and exp -----------------------------
err = zeros(4,1);
err(1) = simu_exp_pair(D_clim, 'clim', WT_met_reac_steady_clim, WT_y_steady_clim,...
                       Jgy_clim, Jeh_clim, cell_clim, precursor_clim, aa_clim, [], protein_clim, lipid_clim, prot_sec_clim, total_protein_con_clim,...
                       D_chem, num_flux, num_y);

%% ------------------ pair up N-lim simu and exp -----------------------------
err(2) = simu_exp_pair(D_cnlim1, 'cnlim1', WT_met_reac_steady_cnlim1, WT_y_steady_cnlim1,...
                       Jgy_cnlim1, Jeh_cnlim1, [], [], [], nh4_cnlim1, protein_cnlim1, [], prot_sec_cnlim1, total_protein_con_cnlim1,...
                       D_chem, num_flux, num_y);
err(3) = simu_exp_pair(D_cnlim2, 'cnlim2', WT_met_reac_steady_cnlim2, WT_y_steady_cnlim2,...
                       Jgy_cnlim2, Jeh_cnlim2, cell_cnlim2, precursor_cnlim2, aa_cnlim2, nh4_cnlim2, protein_cnlim2, [], prot_sec_cnlim2, total_protein_con_cnlim2,...
                       D_chem, num_flux, num_y);
err(4) = simu_exp_pair(D_nlimf, 'nlimf', WT_met_reac_steady_cnlim1, WT_y_steady_cnlim1,...
                       Jgy_nlimf, Jeh_nlimf, cell_nlimf, precursor_nlimf, aa_nlimf, [], protein_nlimf, lipid_nlimf, prot_sec_nlimf, total_protein_con_nlim,...
                       D_chem, num_flux, num_y);

error = sum(sum(err));


% function finds  find the corresponding dilution index for simulation 
function d_index = find_d_index(simu_d,exp_d)
    d_index = zeros(length(exp_d),1);
    for i = 1: length(exp_d)
        if exp_d(i) <= max(simu_d)
            d_index(i) = find(simu_d>=exp_d(i),1);
        else
            d_index(i) = 1; %% the integration shut down
        end
    end
end

function s_error = target_func_val(y_simu,y_data)
    y_simu = reshape(y_simu, numel(y_simu), 1);
    y_data = reshape(y_data, numel(y_simu), 1);
    rel_error    = sum(abs((abs( 2.*y_data./(y_simu + y_data)) - 1))./length(y_data));
    abs_error    = sum(abs(y_data - y_simu)./length(y_data)); 
    s_error      = min(rel_error.^2, abs_error.^2);
end 

function err = simu_exp_pair(D_exp, cond, WT_met_reac_steady, WT_y_steady,...
                             Jgy, Jeh, cell, precursor, aa, nh4, protein, lipid, prot_sec, total_protein_con,...
                             D_chem, num_flux, num_y)
    % find D index
    papers = fieldnames(D_exp);
    for i = 1:length(papers)
        paper_i = papers{i};
        Dclim_idx.(paper_i) = find_d_index(D_chem, D_exp.(paper_i));
    end
    err = 0;
    
    % flux error
    if exist(['Jgy_' cond],'var') 
    papers = fieldnames(Jgy);
        for i = 1: length(papers)
            y_simu = WT_met_reac_steady.flux(Dclim_idx.(papers{i}),num_flux.gy);
            y_data = Jgy.(papers{i});
            err = err + target_func_val(y_simu,y_data);
        end    
    end
    if exist(['Jeh_' cond],'var') 
    papers = fieldnames(Jeh);
        for i = 1: length(papers)
            y_simu = WT_met_reac_steady.flux(Dclim_idx.(papers{i}),num_flux.fe) - WT_met_reac_steady.flux(Dclim_idx.(papers{i}),num_flux.gn);
            y_data = Jeh.(papers{i});
            err = err + target_func_val(y_simu,y_data);
        end 
    end 
    
    % met error
    if exist(['cell_' cond],'var') 
    papers = fieldnames(cell);
        for i = 1: length(papers)
            y_simu = WT_y_steady(Dclim_idx.(papers{i}),end);
            y_data = cell.(papers{i});
            err = err + target_func_val(y_simu,y_data);
        end 
    end 
    if exist(['precursor_' cond],'var') 
    precs = fieldnames(precursor);
        for k = 1: length(precs)
        prec = precs{k};
        papers = fieldnames(precursor.(prec)); 
            for i = 1: length(papers)
                y_simu = WT_y_steady(Dclim_idx.(papers{i}),num_y.pc);
                y_data = precursor.(prec).(papers{i})./precursor.(prec).(papers{i})(1);
                err = err + target_func_val(y_simu,y_data);
            end 
        end 
    end 
    if exist(['aa_' cond],'var') 
    papers = fieldnames(aa.glutamine); 
        for i = 1: length(papers)
            y_simu = WT_y_steady(Dclim_idx.(papers{i}),num_y.aa_in);
            y_data = aa.glutamine.(papers{i})./aa.glutamine.(papers{i})(1);
            err = err + target_func_val(y_simu,y_data);
        end 
    end 
    if exist(['nh4_' cond],'var') 
    papers = fieldnames(nh4);
        for i = 1: length(papers)
            y_simu = WT_y_steady(Dclim_idx.(papers{i}),num_y.nh4);
            y_data = nh4.(papers{i});
            err = err + target_func_val(y_simu,y_data);
        end 
    end
    if exist(['protein_' cond],'var')
    papers = fieldnames(protein); 
        for i = 1: length(papers)
            y_simu = 100.*(total_protein_con(Dclim_idx.(papers{i})).*pt_uM_to_gPerL)./par.rho_cell;
            y_data = protein.(papers{i});
            err = err + target_func_val(y_simu,y_data);
        end 
    end 
    if exist(['lipid_' cond],'var')
    papers = fieldnames(lipid); 
        for i = 1: length(papers)
            y_simu = 100.*WT_y_steady(Dclim_idx.(papers{i}),num_y.lp_e).*lp_uM_to_gPerL./par.rho_cell +lipid_constant_part;
            y_data = lipid.(papers{i});
            err = err + target_func_val(y_simu,y_data);
        end 
    end
    
    % protein error
    
    if exist(['prot_sec_' cond],'var') 
        % gy
        papers = fieldnames(prot_sec.gy); 
        for i = 1: length(papers)
            y_simu = WT_y_steady(Dclim_idx.(papers{i}),num_y.gy);
            y_data = log2(prot_sec.gy.(papers{i})./prot_sec.gy.(papers{i})(1));
            err = err + target_func_val(y_simu,y_data);
        end 
        % fe
        papers = fieldnames(prot_sec.fe); 
        for i = 1: length(papers)
            y_simu = WT_y_steady(Dclim_idx.(papers{i}),num_y.fe);
            y_data = log2(prot_sec.fe.(papers{i})./prot_sec.fe.(papers{i})(1));
            err = err + target_func_val(y_simu,y_data);
        end 
        % gn
        papers = fieldnames(prot_sec.gn); 
        for i = 1: length(papers)
            y_simu = WT_y_steady(Dclim_idx.(papers{i}),num_y.gn);
            y_data = log2(prot_sec.gn.(papers{i})./prot_sec.gn.(papers{i})(1));
            err = err + target_func_val(y_simu,y_data);
        end 
        % mt
        papers = fieldnames(prot_sec.mt); 
        for i = 1: length(papers)
            y_simu = WT_y_steady(Dclim_idx.(papers{i}),num_y.mt);
            y_data = log2(prot_sec.mt.(papers{i})./prot_sec.mt.(papers{i})(1));
            err = err + target_func_val(y_simu,y_data);
        end 
        % as
        papers = fieldnames(prot_sec.as); 
        for i = 1: length(papers)
            y_simu = WT_y_steady(Dclim_idx.(papers{i}),num_y.as);
            y_data = log2(prot_sec.as.(papers{i})./prot_sec.as.(papers{i})(1));
            err = err + target_func_val(y_simu,y_data);
        end 
        % at
        papers = fieldnames(prot_sec.at); 
        for i = 1: length(papers)
            y_simu = WT_y_steady(Dclim_idx.(papers{i}),num_y.at);
            y_data = log2(prot_sec.at.(papers{i})./prot_sec.at.(papers{i})(1));
            err = err + target_func_val(y_simu,y_data);
        end 
        % lp
        papers = fieldnames(prot_sec.lp); 
        for i = 1: length(papers)
            y_simu = WT_y_steady(Dclim_idx.(papers{i}),num_y.lp);
            y_data = log2(prot_sec.lp.(papers{i})./prot_sec.lp.(papers{i})(1));
            err = err + target_func_val(y_simu,y_data);
        end 
        % lo
        papers = fieldnames(prot_sec.lo); 
        for i = 1: length(papers)
            y_simu = WT_y_steady(Dclim_idx.(papers{i}),num_y.lo);
            y_data = log2(prot_sec.lo.(papers{i})./prot_sec.lo.(papers{i})(1));
            err = err + target_func_val(y_simu,y_data);
        end 
        % R
        papers = fieldnames(prot_sec.r); 
        for i = 1: length(papers)
            y_simu = WT_y_steady(Dclim_idx.(papers{i}),num_y.r);
            y_data = log2(prot_sec.r.(papers{i})./prot_sec.r.(papers{i})(1));
            err = err + target_func_val(y_simu,y_data);
        end 
        % Z
        papers = fieldnames(prot_sec.z); 
        for i = 1: length(papers)
            y_simu = WT_y_steady(Dclim_idx.(papers{i}),num_y.z);
            y_data = log2(prot_sec.z.(papers{i})./prot_sec.z.(papers{i})(1));
            err = err + target_func_val(y_simu,y_data);
        end 
    end
end
end
