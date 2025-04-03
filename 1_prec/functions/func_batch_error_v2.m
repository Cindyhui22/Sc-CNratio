function error = func_batch_error_v2(clim_WT_batch_min_t, clim_WT_batch_min_y, ...
                                     clim_WT_batch_YPD_t, clim_WT_batch_YPD_y,...
                                     glucose_climb, ethanol_climb,  cell_climb, prot_sec_climb,...
                                     par, num_y, num_flux, num_prot)
%
% This function calculate the error between experiment and simulation for chemostat
% Relative Error = abs(exp data/sim data - 1)
% Absolute Error = abs(exp data - sim data)
%
%% -------- load simu data and processing --------------------------------
% REZ fraction calculation before normalizing
% clim_WT_batch_min_t, nlim_WT_batch_min_t, clim_WT_batch_YPD_t,
clim_YPD_total_protein_con = sum(clim_WT_batch_YPD_y(:,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);
% nlim_total_protein_con   = sum(nlim_WT_batch_min_y(:,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);
clim_total_protein_con   = sum(clim_WT_batch_min_y(:,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);

% fraction
R_clim_YPD_sim = clim_WT_batch_YPD_y(:,num_y.r).*par.l(num_prot.r)'./clim_YPD_total_protein_con;
Z_clim_YPD_sim = clim_WT_batch_YPD_y(:,num_y.z).*par.l(num_prot.z)'./clim_YPD_total_protein_con;
E_clim_YPD_sim = sum(clim_WT_batch_YPD_y(:,num_y.gy:num_y.lo).*par.l(num_prot.gy:num_prot.lo)',2)./clim_YPD_total_protein_con;     

R_clim_sim = clim_WT_batch_min_y(:,num_y.r).*par.l(num_prot.r)'./clim_total_protein_con;
Z_clim_sim = clim_WT_batch_min_y(:,num_y.z).*par.l(num_prot.z)'./clim_total_protein_con;
E_clim_sim = sum(clim_WT_batch_min_y(:,num_y.gy:num_y.lo).*par.l(num_prot.gy:num_prot.lo)',2)./clim_total_protein_con;
% 
% R_nlim_sim = nlim_WT_batch_min_y(:,num_y.r).*par.l(num_prot.r)'./nlim_total_protein_con;
% Z_nlim_sim = nlim_WT_batch_min_y(:,num_y.z).*par.l(num_prot.z)'./nlim_total_protein_con;
% E_nlim_sim = sum(nlim_WT_batch_min_y(:,num_y.gy:num_y.lo).*par.l(num_prot.gy:num_prot.lo)',2)./nlim_total_protein_con;



%% normalizing  
clim_WT_batch_min_y(:,num_y.pc)   = clim_WT_batch_min_y(:,num_y.pc)   ./ clim_WT_batch_min_y(1,num_y.pc);
% nlim_WT_batch_min_y(:,num_y.pc)   = nlim_WT_batch_min_y(:,num_y.pc)   ./ nlim_WT_batch_min_y(1,num_y.pc);
%clim_WT_batch_YPD_y(:,num_y.pc) = clim_WT_batch_YPD_y(:,num_y.pc) ./ clim_WT_batch_YPD_y(1,num_y.pc);

clim_WT_batch_min_y(:,num_y.aa_in)   = clim_WT_batch_min_y(:,num_y.aa_in)   ./ clim_WT_batch_min_y(1,num_y.aa_in);
% nlim_WT_batch_min_y(:,num_y.aa_in)   = nlim_WT_batch_min_y(:,num_y.aa_in)   ./ nlim_WT_batch_min_y(1,num_y.aa_in);
%clim_WT_batch_YPD_y(:,num_y.aa_in) = clim_WT_batch_YPD_y(:,num_y.aa_in) ./ clim_WT_batch_YPD_y(1,num_y.aa_in);

clim_WT_batch_min_y(:,num_y.ae)   = clim_WT_batch_min_y(:,num_y.ae)   ./ clim_WT_batch_min_y(1,num_y.ae);
% nlim_WT_batch_min_y(:,num_y.ae)   = nlim_WT_batch_min_y(:,num_y.ae)   ./ nlim_WT_batch_min_y(1,num_y.ae);
%clim_WT_batch_YPD_y(:,num_y.ae) = clim_WT_batch_YPD_y(:,num_y.ae) ./ clim_WT_batch_YPD_y(1,num_y.ae);

% proteins frac (sim)
%_% clim_WT_batch_min_y(:,[num_y.r: num_y.sd]) = (clim_WT_batch_min_y(:,[num_y.r: num_y.sd]).*par.l')./ clim_total_protein_con;
%_% nlim_WT_batch_min_y(:,[num_y.r: num_y.sd]) = (nlim_WT_batch_min_y(:,[num_y.r: num_y.sd]).*par.l')./ nlim_total_protein_con;
%_% %clim_WT_batch_YPD_y(:,[num_y.r: num_y.sd]) = (clim_WT_batch_YPD_y(:,[num_y.r: num_y.sd]).*par.l')./ clim_total_protein_con;
clim_WT_batch_min_y(:,[num_y.r: num_y.lo]) = (clim_WT_batch_min_y(:,[num_y.r: num_y.lo]).*par.l')./ clim_total_protein_con;
% nlim_WT_batch_min_y(:,[num_y.r: num_y.lo]) = (nlim_WT_batch_min_y(:,[num_y.r: num_y.lo]).*par.l')./ nlim_total_protein_con;
%clim_WT_batch_YPD_y(:,[num_y.r: num_y.lo]) = (clim_WT_batch_YPD_y(:,[num_y.r: num_y.lo]).*par.l')./ clim_total_protein_con;

% protein log2fc 
%_% clim_WT_batch_min_y(:,[num_y.r: num_y.sd])   = log2(clim_WT_batch_min_y(:,[num_y.r: num_y.sd])   ./ clim_WT_batch_min_y(1,[num_y.r: num_y.sd]));
%_% nlim_WT_batch_min_y(:,[num_y.r: num_y.sd])   = log2(nlim_WT_batch_min_y(:,[num_y.r: num_y.sd])   ./ nlim_WT_batch_min_y(1,[num_y.r: num_y.sd]));
%_% %clim_WT_batch_YPD_y(:,[num_y.r: num_y.sd]) = log2(clim_WT_batch_YPD_y(:,[num_y.r: num_y.sd]) ./ clim_WT_batch_YPD_y(1,[num_y.r: num_y.sd]));
clim_WT_batch_min_y(:,[num_y.r: num_y.lo])   = log2(clim_WT_batch_min_y(:,[num_y.r: num_y.lo])   ./ clim_WT_batch_min_y(1,[num_y.r: num_y.lo]));
% nlim_WT_batch_min_y(:,[num_y.r: num_y.lo])   = log2(nlim_WT_batch_min_y(:,[num_y.r: num_y.lo])   ./ nlim_WT_batch_min_y(1,[num_y.r: num_y.lo]));
%clim_WT_batch_YPD_y(:,[num_y.r: num_y.lo]) = log2(clim_WT_batch_YPD_y(:,[num_y.r: num_y.lo]) ./ clim_WT_batch_YPD_y(1,[num_y.r: num_y.lo]));
%% ------------------ pair up C-lim simu and exp -----------------------------
err = zeros(2,1);
% bartolomeo
pro_time = [8.8, 12.5, 20];
batch_time = [0.00	
            2.01	
            2.97	
            4.01	
            5.02	
            5.99	
            6.99	
            7.99	
            9.00	
            10.00	
            11.00	
            11.97	
            13.01	
            13.98	
            14.98	
            15.95	
            16.99	
            17.99	
            19.00	
            20.00	
            ]; 
err(1) = simu_exp_pair(batch_time, pro_time, 'bartolomeo', clim_WT_batch_min_t, clim_WT_batch_min_y,...
                             glucose_climb, ethanol_climb,  cell_climb, prot_sec_climb,...
                             R_clim_sim, E_clim_sim, Z_clim_sim, ...
                             num_flux, num_y);
% murphy 
batch_time = [0             
            4.917355372	
            6.94214876	
            8.925619835	
            10.95041322	
            12.9338843	
            14.95867769	
            16.98347107	
            24.95867769	
            28.96694215	
            32.97520661
            ];
pro_time = [5
7
9
11
13
15
17
];
err(2) = simu_exp_pair(batch_time, pro_time, 'murphy', clim_WT_batch_YPD_t, clim_WT_batch_YPD_y,...
                             glucose_climb, ethanol_climb,  cell_climb, prot_sec_climb,...
                             R_clim_YPD_sim, E_clim_YPD_sim, Z_clim_YPD_sim, ...
                             num_flux, num_y);
error = sum(sum(err));

% function finds  find the corresponding dilution index for simulation 
function d_index = find_d_index(simu_d,exp_d)
    simu_d = reshape(simu_d, numel(simu_d), 1);
    exp_d = reshape(exp_d, numel(exp_d), 1);
    d_index = zeros(length(exp_d),1);
    for i = 1: length(exp_d)
        if exp_d(i) <= max(simu_d)
            d_index(i) = find(simu_d>=exp_d(i),1);
        else
            d_index(i) = length(simu_d); %% the integration shut down
        end
    end
end

function s_error = target_func_val(y_simu,y_data)
    y_simu = reshape(y_simu, numel(y_simu), 1);
    y_data = reshape(y_data, numel(y_simu), 1);
    rel_error    = sum(abs((abs( (2.*y_data + eps(1))./(y_simu + y_data + eps(1))) - 1))./length(y_data));
    abs_error    = sum(abs(y_data - y_simu)./length(y_data)); 
    s_error      = min(rel_error.^2, abs_error.^2);
end 

function err = simu_exp_pair(t_exp, t_pro, paper, simu_t, simu_y,...
                             glucose_climb, ethanol_climb,  cell_climb, prot_sec_climb,...
                             R_clim_YPD_sim, E_clim_YPD_sim, Z_clim_YPD_sim, ...
                             num_flux, num_y)
    err = 0;
    % find t index
    t_idx = find_d_index(simu_t, t_exp);
    tp_idx = find_d_index(simu_t, t_pro);
    

    % met error
    if isfield(glucose_climb, paper)
        y_simu = simu_y(t_idx, num_y.gl);
        y_data = glucose_climb.(paper)(:,2);
        err = err + target_func_val(y_simu,y_data);
    end
    if isfield(ethanol_climb, paper)
        y_simu = simu_y(t_idx, num_y.eh);
        y_data = ethanol_climb.(paper)(:,2);
        err = err + target_func_val(y_simu,y_data);
    end
    if isfield(cell_climb, paper)
        y_simu = simu_y(t_idx, end);
        if strcmp(paper, 'murphy')
            y_data = [cell_climb.(paper)(1,2), cell_climb.(paper)(3:end,2)'];
        else
            y_data = cell_climb.(paper)(:,2);
        end
        err = err + target_func_val(y_simu,y_data);
    end
   
    % protein error
    % REZ fraction
    y_simu = R_clim_YPD_sim(tp_idx);
    y_data = prot_sec_climb.r.(paper);
    err = err + target_func_val(y_simu,y_data);

    y_simu = E_clim_YPD_sim(tp_idx);
    y_data = prot_sec_climb.e.(paper);
    err = err + target_func_val(y_simu,y_data);

    y_simu = Z_clim_YPD_sim(tp_idx);
    y_data = prot_sec_climb.z.(paper);
    err = err + target_func_val(y_simu,y_data);

    % gy, fe, gn, mt, as, at, lp, lo, R, Z
    y_simu = simu_y(tp_idx,num_y.gy);
    y_data = log2(prot_sec_climb.gy.(paper)./prot_sec_climb.gy.(paper)(1));
    err = err + target_func_val(y_simu,y_data);

    y_simu = simu_y(tp_idx,num_y.fe);
    y_data = log2(prot_sec_climb.fe.(paper)./prot_sec_climb.fe.(paper)(1));
    err = err + target_func_val(y_simu,y_data);

    y_simu = simu_y(tp_idx,num_y.gn);
    y_data = log2(prot_sec_climb.gn.(paper)./prot_sec_climb.gn.(paper)(1));
    err = err + target_func_val(y_simu,y_data);
    
    y_simu = simu_y(tp_idx,num_y.mt);
    y_data = log2(prot_sec_climb.mt.(paper)./prot_sec_climb.mt.(paper)(1));
    err = err + target_func_val(y_simu,y_data);

    y_simu = simu_y(tp_idx,num_y.as);
    y_data = log2(prot_sec_climb.as.(paper)./prot_sec_climb.as.(paper)(1));
    err = err + target_func_val(y_simu,y_data);

    y_simu = simu_y(tp_idx,num_y.at);
    y_data = log2(prot_sec_climb.at.(paper)./prot_sec_climb.at.(paper)(1));
    err = err + target_func_val(y_simu,y_data);

    y_simu = simu_y(tp_idx,num_y.lp);
    y_data = log2(prot_sec_climb.lp.(paper)./prot_sec_climb.lp.(paper)(1));
    err = err + target_func_val(y_simu,y_data);

    y_simu = simu_y(tp_idx,num_y.lo);
    y_data = log2(prot_sec_climb.gy.(paper)./prot_sec_climb.lo.(paper)(1));
    err = err + target_func_val(y_simu,y_data);

    y_simu = simu_y(tp_idx,num_y.r);
    y_data = log2(prot_sec_climb.r.(paper)./prot_sec_climb.r.(paper)(1));
    err = err + target_func_val(y_simu,y_data);

    y_simu = simu_y(tp_idx,num_y.z);
    y_data = log2(prot_sec_climb.z.(paper)./prot_sec_climb.z.(paper)(1));
    err = err + target_func_val(y_simu,y_data);

end
end
