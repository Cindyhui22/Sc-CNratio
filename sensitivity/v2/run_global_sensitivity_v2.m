clc; close all; clear;  
addpath('../')
addpath(genpath('./'))


load_simulation_setting;
data_unit_conv; 
plot_num; 
plt_labels = plot_labels; 
plt_colors = plot_colors;

%% initialize by model para
par = read_parametersv2; 

% get parameters names
parName         = fieldnames(par);
parName(find(strcmp(parName,'K_3AT'))) = []; % not a base model parameters
parName(find(strcmp(parName,'c_3AT'))) = []; % not a base model parameters
parName(find(strcmp(parName,'K_mgl'))) = []; % not a base model parameters


% stack the SI sheets into a single table
par_type = {'ori_nonfree','ori_free','ori_adapt_nonfree','ori_adapt_free', 'new_nonfree', 'new_free'};
par_all_type_table = [];
for i = 1:length(par_type)
par_each_type       = readtable('SI_para_type.xlsx','Sheet', par_type{i});
par_each_type_table = [cell2table(par_each_type{:,1},'VariableNames',{'MATLAB_parname'}) ...
cell2table(par_each_type{:,end},'VariableNames',{'Paramter_category'})];
par_all_type_table  = [par_all_type_table; par_each_type_table];
par_type_size(i)    = numel(par_each_type{:,1}); % - 1 because of first row is empty 
end 
table([par_type';'total par'], [par_type_size'; sum(par_type_size)]) 

par_info       = par_all_type_table; 
parCategory    = cell(length(parName),1); 
par_info_row = 1; 
for i = 1:length(parName) % organize in MATLAB order
if ~isempty(find(strcmp(parName{i},par_info.MATLAB_parname)))
par_info_row       = find(strcmp(parName{i},par_info.MATLAB_parname)); 
parCategory{i}     = par_info.Paramter_category{par_info_row};
else 
parCategory{i}     = 'none';
disp('Did not foud this read_parametersv2 parameter in excel')
disp(parName{i})
end 
end 


none_indices = strcmp(parCategory, 'none');
if any(none_indices)
    none_params = parName(none_indices);  % Get the names of the parameters
    warning(['The following parameters have no category assigned and will be removed: ', strjoin(none_params, ', ')]);
    parCategory(none_indices) = [];
    parName(none_indices) = [];
end

params = {};
orif_idx = strcmp(parCategory, 'ori_free');
oaf_idx = strcmp(parCategory, 'ori_adapt_free');
newf_idx = strcmp(parCategory, 'new_free');
if any(orif_idx)
    % Get the names of the parameters
    param_orif = parName(orif_idx);  
    for i = 1:numel(param_orif)
        params{end+1} = {['par.' param_orif{i}], par.(param_orif{i}), 0};
    end
end
if any(oaf_idx)
    % Get the names of the parameters
    param_oaf = parName(oaf_idx);  
    for i = 1:numel(param_oaf)
        params{end+1} = {['par.' param_oaf{i}], par.(param_oaf{i}), 0};
    end
end
if any(newf_idx)
    % Get the names of the parameters
    param_newf = parName(newf_idx);  
    for i = 1:numel(param_newf)
        params{end+1} = {['par.' param_newf{i}], par.(param_newf{i}), 0};
    end
end
params = params';
para_f = [param_orif; param_oaf; param_newf];
%% -------------------------------------------------------------------------------------------------------  
%  Read in & Organize Exp for comparison;
%  ------------------------------------------------------------------------------------------------------- 
func_exp_load;

%% MCMC model 
model.ssfun = @(params, data) fun_mcmc_loss(params, data, para_f,...
                                            D_clim, Jgy_clim, Jeh_clim, cell_clim, precursor_clim, aa_clim, protein_clim, lipid_clim, prot_sec_clim, ...
                                            D_cnlim1, Jgy_cnlim1, Jeh_cnlim1, nh4_cnlim1, protein_cnlim1, prot_sec_cnlim1, ...
                                            D_cnlim2, Jgy_cnlim2, Jeh_cnlim2, cell_cnlim2, precursor_cnlim2, aa_cnlim2, nh4_cnlim2, protein_cnlim2, prot_sec_cnlim2, ...
                                            D_nlimf, Jgy_nlimf, Jeh_nlimf, cell_nlimf, precursor_nlimf, aa_nlimf, protein_nlimf, lipid_nlimf, prot_sec_nlimf, ...
                                            glucose_climb, ethanol_climb,  cell_climb, prot_sec_climb);   % fitting & loss calculation
model.sigma2 = 0.01^2;          % initialize N(0, sigma)
model.nbatch = 4+6+6+15+4+30+2+5+... % total number of dataset
               3+6+6+11+2+60+1+5+...
               2*16;    
model.N = (5+4+5+5)+...         % total observation points
          (5+4+5+5+9+1)+...
          (5+4+5+5+9+1)+...
          (7*5+4*4+2*5+2*5)+...
          (5+4+5+5)+...
          (10*(5+9+1))+...
          (5+7)+...
          (5+7+9+7+1)+...
            (5+5+6)+...
            (5+6+1+1+1+6)+...
            (5+6+1+1+1+6)+...
            (7*5+4*6)+...
            (5+6)+...
            (10*(5+1+1+1+6))+...
            (5)+...
            (5+1+1+1+6)+...
           (3*13+20*3)+(7*13+11*3); 


data = 0; %paased in data info by options, no need to give data
Markov_chain = [];
parfor chain_run= 1:4  
    dirname = sprintf('./my_mcmc_data%d/', chain_run); %%%Change this
    if ~exist(dirname, 'dir')
        % Directory does not exist, create it
        [status, msg] = mkdir(dirname);
        disp(['Directory "', dirname, '" created successfully.']);
    else
        % Directory already exists
        disp(['Directory "', dirname, '" already exists.']);
    end
    options = struct();
    options.chainfile = sprintf('my_mcmc_chain%d.bin', chain_run); %%%Change this
    options.sschainfile = sprintf('my_mcmc_sschain%d.bin', chain_run); %%%Change this
    options.s2chainfile = sprintf('my_mcmc_s2chain%d.bin', chain_run); %%%Change this
    options.savedir = dirname; %or other directory
    options.updatesigma = 1;        % update error variance N(0, sigma)
    
    options.verbosity   = 1;        % information printed option
    options.savesize = 1000;     %%%Change this
    options.nsimu = 10000;       %%%Change this
    options.burnintime = 0;
    
    [res,chain,s2chain] = mcmcrun(model,data,params,options); 
    % save(sprintf('result_mcmc%d.mat', chain_run), 'chain', 'res'); %%%Change this
    filename = sprintf('result_mcmc%d.mat', chain_run);
    parsave(filename, chain, res);
    % Markov_chain = [Markov_chain; chain];
end
% save('Markov_chains.mat', 'Markov_chain');
%% Read file
% addpath('./my_mcmc_data1/')
% chain = readbin(options.chainfile, 0, 1);
% sschain = readbin(options.sschainfile, 0, 1);
% if options.updatesigma
%     s2chain = readbin(options.s2chainfile, 0, 1);
% end

%% plotting
% %------------------ Chain plot------------------------------------------
% addpath('./sensitivity/loaded_functions/mcmcstat')
% figure();
% mcmcplot(Markov_chain,[],res,'chainpanel');
% res.names = para_f;
% 
% %------------------- Violin plot ---------------------------------------
% rmpath('./sensitivity/loaded_functions/mcmcstat')
% addpath('./sensitivity/loaded_functions/Violinplot-Matlab-master/');
% num_per_plot = 12; 
% figure();
% % box plot set up
% % y_lims = [10^-2 10^11; 10^-4 10^4; 10^-3 10^4];
% 
% for k = 1:ceil(numel(para_f)/num_per_plot)
%     if k < max(ceil(numel(para_f)/num_per_plot))
%     index_range = [(k-1)*num_per_plot+1:k*num_per_plot];
%     else 
%     index_range = [(k-1)*num_per_plot+1:numel(para_f)];
%     end 
%     chain_k = Markov_chain(:,index_range);
%     params_names_k = para_f(index_range);
%     params_k = params(index_range);
%     % y_lim = y_lims(k,:); 
%     yticks(10.^[log10(min(ylim)):2:log10(max(ylim))])
% 
%     subplot(ceil(numel(para_f)/num_per_plot),1, k) 
%     violinplot(chain_k,[],'ShowData',false);
%     ylabel('Values');
%     set(gca,'XTickLabel',params_names_k)  
%     set(gca, 'YMinorTick', 'off');
%     set(gca, 'XMinorTick', 'off');
%     bx = gca;
%     hold on; 
%     x_axis_space = (bx.XTick(2) - bx.XTick(1))/2;   
%     for i = 1: numel(index_range)
%         % plot initial value on top
%         plot(bx.XTick(i), params_k{i}{2},'o','MarkerFaceColor','r')
%     end
% end
% 
% 
% %% Statistics
% % mean
% marc_mean = mean(Markov_chain, 1);
% % median
% marc_median = median(Markov_chain, 1);
% % std
% marc_std = std(Markov_chain,0, 1);
% 
% % 95% equitailed calculation
% marc_sort = sort(Markov_chain, 1);
% marc_size = size(Markov_chain, 1);
% marc_CI2 = [marc_sort(marc_size*0.025, :); marc_sort(marc_size*0.975, :)];

%% Predictions
% rmpath('./sensitivity/loaded_functions/mcmcstat')
% addpath('./sensitivity/loaded_functions/Violinplot-Matlab-master/');
% data = [];
% n_sample=2000;  %%%Change this
% [protein_mcmc, protein_clim, protein_nlim] = mcmc_pred(Markov_chain, n_sample, para_f);
% % figure();
% % violinplot(protein_mcmc,[],'ShowData',false);
% % ylabel('C/N Protein fraction');
% % set(gca,'XTickLabel',{'r','z','gy','fe','gn','mt','as','at','lp','lo'});         
% % set(gca, 'YMinorTick', 'off');
% % set(gca, 'XMinorTick', 'off');
% % set(gca, 'YScale','log');
% % ylim([1e-3, 1e3])
% % axis square
% 
% % figure();%scatter plot
% % savefig(gcf, 'Protein_mcmc.fig');
% save('Protein_mcmc.mat', 'protein_mcmc','protein_clim', 'protein_nlim');



%%
function [protein_mcmc, protein_clim, protein_nlim] = mcmc_pred(Markov_chain, n_sample, para_f)
    protein_mcmc = zeros(n_sample, 10);
    protein_clim = zeros(n_sample, 10);
    protein_nlim = zeros(n_sample, 10);
    numRows = size(Markov_chain, 1);
    parfor i=1:n_sample %%%change this
        rand_idx = randi(numRows);
        para_i = Markov_chain(rand_idx,:);
        [pred_ratio, pred_clim, pred_nlim] = protein_predict(para_i,  para_f);
        protein_mcmc(i,:) = pred_ratio';
        protein_clim(i,:) = pred_clim;
        protein_nlim(i,:) = pred_nlim;
    end
end


function [pred_ratio, pred_clim, pred_nlim] = protein_predict(para_i,  para_f)
    load_simulation_setting;
    data_unit_conv; 
    plot_num; 
    par = read_par(para_i, para_f); 

    % run chemostat and calculate sim-exp error; with unperturbed reference parameter; use as reference later. 
    [WT_y_steady_clim,   ~, ...
     WT_y_steady_cnlim1, ~, ...
     WT_y_steady_cnlim2, ~, ...
     WT_y_steady_nlim,   ~,   ~] = func_chemostat_simulation(par);

    % proteins frac (sim)
    total_protein_con_clim = sum(WT_y_steady_clim(:,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);
    total_protein_con_nlim = sum(WT_y_steady_nlim(:,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);
    % total_protein_con_cnlim1 = sum(WT_y_steady_cnlim1(:,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);
    % total_protein_con_cnlim2 = sum(WT_y_steady_cnlim2(:,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);
    % 
    WT_y_steady_clim(:,[num_y.r: num_y.lo]) = (WT_y_steady_clim(:,[num_y.r: num_y.lo]).*par.l')./ total_protein_con_clim;
    WT_y_steady_nlim(:,[num_y.r: num_y.lo]) = (WT_y_steady_nlim(:,[num_y.r: num_y.lo]).*par.l')./ total_protein_con_nlim;
    % WT_y_steady_cnlim1(:,[num_y.r: num_y.lo]) = (WT_y_steady_cnlim1(:,[num_y.r: num_y.lo]).*par.l')./ total_protein_con_cnlim1;
    % WT_y_steady_cnlim2(:,[num_y.r: num_y.lo]) = (WT_y_steady_cnlim2(:,[num_y.r: num_y.lo]).*par.l')./ total_protein_con_cnlim2;
    
    WT_y_steady_clim_frac = WT_y_steady_clim;
    WT_y_steady_nlim_frac = WT_y_steady_nlim;
    % WT_y_steady_cnlim1_frac = WT_y_steady_cnlim1;
    % WT_y_steady_cnlim2_frac = WT_y_steady_cnlim2;


    % % REZ fraction calculation before normalizing
    % % fraction
    % R_chemo_sim_clim = WT_y_steady_clim(:,num_y.r).*par.l(num_prot.r)' ./ total_protein_con_clim;
    % Z_chemo_sim_clim = WT_y_steady_clim(:,num_y.z).*par.l(num_prot.z)' ./ total_protein_con_clim;
    % E_chemo_sim_clim = sum(WT_y_steady_clim(:,num_y.gy:num_y.lo).*par.l(num_prot.gy:num_prot.lo)',2) ./ total_protein_con_clim;
    % 
    % R_chemo_sim_nlim = WT_y_steady_nlim(:,num_y.r).*par.l(num_prot.r)' ;
    % Z_chemo_sim_nlim = WT_y_steady_nlim(:,num_y.z).*par.l(num_prot.z)' ;
    % E_chemo_sim_nlim = sum(WT_y_steady_nlim(:,num_y.gy:num_y.lo), 2);
    % % relative fold change
    % R_chemo_sim_fc = R_chemo_sim_clim./R_chemo_sim_clim(1); 
    % Z_chemo_sim_fc = Z_chemo_sim_clim./Z_chemo_sim_clim(1);
    % E_chemo_sim_fc = E_chemo_sim_clim./E_chemo_sim_clim(1);
    % 
    % R_chemo_sim_fc_nlim = R_chemo_sim_nlim./R_chemo_sim_nlim(1); 
    % Z_chemo_sim_fc_nlim = Z_chemo_sim_nlim./Z_chemo_sim_nlim(1);
    % E_chemo_sim_fc_nlim = E_chemo_sim_nlim./E_chemo_sim_nlim(1);
    % 
    % % normalize precursor, amino acid, atp, protein with respect to first data point 
    % D_exp = 0.05; % first D in exp
    % [norm_D, norm_indx] = min(abs(D_chem - D_exp)); 
    % % protein log2fc 
    % WT_y_steady_clim(:,[num_y.r: num_y.lo]) = log2(WT_y_steady_clim(:,[num_y.r: num_y.lo])./WT_y_steady_clim(norm_indx,[num_y.r: num_y.lo])); 
    % WT_y_steady_nlim(:,[num_y.r: num_y.lo]) = log2(WT_y_steady_nlim(:,[num_y.r: num_y.lo])./WT_y_steady_nlim(norm_indx,[num_y.r: num_y.lo])); 
    % WT_y_steady_cnlim1(:,[num_y.r: num_y.lo]) = log2(WT_y_steady_cnlim1(:,[num_y.r: num_y.lo])./WT_y_steady_cnlim1(norm_indx,[num_y.r: num_y.lo])); 
    % WT_y_steady_cnlim2(:,[num_y.r: num_y.lo]) = log2(WT_y_steady_cnlim2(:,[num_y.r: num_y.lo])./WT_y_steady_cnlim2(norm_indx,[num_y.r: num_y.lo])); 
    index = 12;
    pred_ratio = [100*WT_y_steady_clim_frac(index,num_y.r)./100*WT_y_steady_nlim_frac(index,num_y.r);
                  100*WT_y_steady_clim_frac(index,num_y.z)./100*WT_y_steady_nlim_frac(index,num_y.z);
                  100*WT_y_steady_clim_frac(index,num_y.gy)./100*WT_y_steady_nlim_frac(index,num_y.gy);
                  100*WT_y_steady_clim_frac(index,num_y.fe)./100*WT_y_steady_nlim_frac(index,num_y.fe);
                  100*WT_y_steady_clim_frac(index,num_y.gn)./100*WT_y_steady_nlim_frac(index,num_y.gn);
                  100*WT_y_steady_clim_frac(index,num_y.mt)./100*WT_y_steady_nlim_frac(index,num_y.mt);
                  100*WT_y_steady_clim_frac(index,num_y.as)./100*WT_y_steady_nlim_frac(index,num_y.as);
                  100*WT_y_steady_clim_frac(index,num_y.at)./100*WT_y_steady_nlim_frac(index,num_y.at);
                  100*WT_y_steady_clim_frac(index,num_y.lp)./100*WT_y_steady_nlim_frac(index,num_y.lp);
                  100*WT_y_steady_clim_frac(index,num_y.lo)./100*WT_y_steady_nlim_frac(index,num_y.lo);];
    pred_clim = 100*WT_y_steady_clim_frac(index,num_y.r:num_y.lo);
    pred_nlim = 100*WT_y_steady_nlim_frac(index,num_y.r:num_y.lo);
end


function error = fun_mcmc_loss(params, data, para_f,...
                               D_clim, Jgy_clim, Jeh_clim, cell_clim, precursor_clim, aa_clim, protein_clim, lipid_clim, prot_sec_clim, ...
                               D_cnlim1, Jgy_cnlim1, Jeh_cnlim1, nh4_cnlim1, protein_cnlim1, prot_sec_cnlim1, ...
                               D_cnlim2, Jgy_cnlim2, Jeh_cnlim2, cell_cnlim2, precursor_cnlim2, aa_cnlim2, nh4_cnlim2, protein_cnlim2, prot_sec_cnlim2, ...
                               D_nlimf, Jgy_nlimf, Jeh_nlimf, cell_nlimf, precursor_nlimf, aa_nlimf, protein_nlimf, lipid_nlimf, prot_sec_nlimf, ...
                               glucose_climb, ethanol_climb,  cell_climb, prot_sec_climb)
    load_simulation_setting;
    data_unit_conv; 
    plot_num; 
    par = read_par(params, para_f); 
    
    [WT_y_steady_clim,   WT_met_reac_steady_clim, ...
     WT_y_steady_cnlim1, WT_met_reac_steady_cnlim1, ...
     WT_y_steady_cnlim2, WT_met_reac_steady_cnlim2, ...
     WT_y_steady_nlim,   WT_met_reac_steady_nlim,   D_chem] = func_chemostat_simulation(par);
    
    error_chem = func_chemo_error_v2(WT_y_steady_clim,   WT_met_reac_steady_clim, ...
                                    WT_y_steady_cnlim1, WT_met_reac_steady_cnlim1, ...
                                    WT_y_steady_cnlim2, WT_met_reac_steady_cnlim2, ...
                                    WT_y_steady_nlim,   WT_met_reac_steady_nlim, ...
                                    D_clim, Jgy_clim, Jeh_clim, cell_clim, precursor_clim, aa_clim, protein_clim, lipid_clim, prot_sec_clim, ...
                                    D_cnlim1, Jgy_cnlim1, Jeh_cnlim1, nh4_cnlim1, protein_cnlim1, prot_sec_cnlim1, ...
                                    D_cnlim2, Jgy_cnlim2, Jeh_cnlim2, cell_cnlim2, precursor_cnlim2, aa_cnlim2, nh4_cnlim2, protein_cnlim2, prot_sec_cnlim2, ...
                                    D_nlimf, Jgy_nlimf, Jeh_nlimf, cell_nlimf, precursor_nlimf, aa_nlimf, protein_nlimf, lipid_nlimf, prot_sec_nlimf, ...
                                    num_y, num_flux, num_prot, ...
                                    pt_uM_to_gPerL, lp_uM_to_gPerL, lipid_sc_c,...
                                    D_chem, par);

    [clim_WT_batch_min_t, clim_WT_batch_min_y, ...
     clim_WT_batch_YPD_t, clim_WT_batch_YPD_y] = func_batch_simulation(par);
    
    error_batch = func_batch_error_v2(clim_WT_batch_min_t, clim_WT_batch_min_y, ...
                                         clim_WT_batch_YPD_t, clim_WT_batch_YPD_y,...
                                         glucose_climb, ethanol_climb,  cell_climb, prot_sec_climb,...
                                         par, num_y, num_flux, num_prot);
    
    
    error = error_chem + error_batch; 
end


function par = read_par(p, para_f)
%--------------------------------------------------
% running global sa for most sensitive fitted parameter
%--------------------------------------------------
par = read_parametersv2; 

for i = 1:numel(para_f)
    par.(para_f{i})= p(i);
end
par.k_tj = [par.K_t_r;
            par.K_t_z;
            par.K_t_gy;
            par.K_t_fe;
            par.K_t_gn;
            par.K_t_mt;
            par.K_t_as;
            par.K_t_at;
            par.K_t_lp;
            par.K_t_lo;];
end 

function parsave(filename, chain, res)
    save(filename, 'chain', 'res'); 
end



