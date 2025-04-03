%% ----------- load exp data ---------------------------------------------
[D_clim, cell_clim, glucose_clim, Jgy_clim, Jeh_clim,... 
precursor_clim, atp_clim, aa_clim,  prot_sec_clim, ...
lipid_clim, protein_clim,carbo_clim, RNA_clim, paperinfo_clim,...
color_clim, shape_clim, gas_clim] = SC_C_lim_chem();   

[D_nlim, cell_nlim, glucose_nlim, Jgy_nlim, Jeh_nlim, nh4_nlim, gas_nlim, ... 
precursor_nlim, atp_nlim, aa_nlim,  prot_sec_nlim, ...
lipid_nlim, protein_nlim,carbo_nlim, RNA_nlim, paperinfo_nlim,...                    
color_nlim, shape_nlim, Jnh4_nlim]  = SC_N_lim_chem(); 

[glucose_climb, ethanol_climb,  cell_climb, ...
precursor_climb, atp_climb, aa_climb, prot_sec_climb,  ...
tag_climb, pl_climb, paperinfo_climb,...
color_climb, shape_climb] =  SC_C_lim_batch(); 

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