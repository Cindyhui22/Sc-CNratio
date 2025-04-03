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

%% MCMC model 
model.ssfun = @(params, data) fun_mcmc_loss(params, data, para_f);   % fitting & loss calculation
model.sigma2 = 0.01^2;          % initialize N(0, sigma)
model.nbatch = 4+6+6+15+4+30+2+5+... % total number of dataset
               3+6+6+11+2+60+1+5;    
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
            (5+1+1+1+6); 


options.updatesigma = 1;        % update error variance N(0, sigma)
options.verbosity   = 1;        % information printed option
options.nsimu = 20000;       %%%Change this
options.burnintime = 1000;   %%%Change this

data = 0; %paased in data info by options, no need to give data
Markov_chain = [];
parfor chain_run= 1:4  
    [res,chain,s2chain] = mcmcrun(model,data,params,options); 
    % save(sprintf('result_mcmc%d.mat', chain_run), 'chain', 'res');
    filename = sprintf('result_mcmc%d.mat', chain_run);
    parsave(filename, chain, res);
    Markov_chain = [Markov_chain; chain];
end
save('Markov_chains.mat', 'Markov_chain');

%% plotting
%------------------ Chain plot------------------------------------------
addpath('./sensitivity/loaded_functions/mcmcstat')
figure();
mcmcplot(Markov_chain,[],res,'chainpanel');
res.names = para_f;

%------------------- Violin plot ---------------------------------------
rmpath('./sensitivity/loaded_functions/mcmcstat')
addpath('./sensitivity/loaded_functions/Violinplot-Matlab-master/');
num_per_plot = 15; 
figure();
% box plot set up
y_lims = [10^-4 10^11;10^-4 10^11; 10^-4 10^11; 10^-4 10^11];

for k = 1:ceil(numel(para_f)/num_per_plot)
    if k < max(ceil(numel(para_f)/num_per_plot))
    index_range = [(k-1)*num_per_plot+1:k*num_per_plot];
    else 
    index_range = [(k-1)*num_per_plot+1:numel(para_f)];
    end 
    chain_k = Markov_chain(:,index_range);
    params_names_k = para_f(index_range);
    params_k = params(index_range);
    y_lim = y_lims(k,:); 
    yticks(10.^[log10(min(y_lim)):2:log10(max(y_lim))])

    subplot(ceil(numel(para_f)/num_per_plot),1, k) 
    violinplot(chain_k,[],'ShowData',false);
    ylabel('Values');
    set(gca,'YScale','log')
    set(gca,'XTickLabel',params_names_k)  
    set(gca, 'YMinorTick', 'off');
    set(gca, 'XMinorTick', 'off');
    bx = gca;
    hold on; 
    x_axis_space = (bx.XTick(2) - bx.XTick(1))/2;   
    for i = 1: numel(index_range)
        % plot initial value on top
        plot(bx.XTick(i), params_k{i}{2},'o','MarkerFaceColor','r')
    end
end

%% Statistics
% mean
marc_mean = mean(Markov_chain, 1);
% median
marc_median = median(Markov_chain, 1);
% std
marc_std = std(Markov_chain,0, 1);

% 95% equitailed calculation
marc_sort = sort(Markov_chain, 1);
marc_size = size(Markov_chain, 1);
marc_CI2 = [marc_sort(marc_size*0.025, :); marc_sort(marc_size*0.975, :)];

%% Predictions
rmpath('./sensitivity/loaded_functions/mcmcstat')
addpath('./sensitivity/loaded_functions/Violinplot-Matlab-master/');
load('Markov_chains.mat');
load('result_mcmc1.mat');
data = [];
n_sample=2000;  %%%Change this
[protein_mcmc, protein_clim, protein_nlim] = mcmc_pred(Markov_chain, n_sample, para_f);
save('Protein_mcmc_pred.mat', 'protein_mcmc','protein_clim', 'protein_nlim');
%%
load('Protein_mcmc_pred.mat')
protein_colors = {'#FED304';
'#0FA64A';
'#0778A9';
'#3D57A7';
'#6E4EA0';
'#773A96';
'#A92179';
'#EE3324';
'#F47C20';
'#FBA919';
};
protein_shape = {"o", "square", "diamond", "^", "x", "v","hexagram", "*", "pentagram", "+"};
% protein_colors = {    [0.9961    0.8275    0.0157];  % #FED304
%                       [0.0627    0.6510    0.2902];  % #0FA64A
%                       [0.0275    0.4706    0.6235];  % #0778A9
%                       [0.2431    0.3412    0.6510];  % #3D57A7
%                       [0.4314    0.3255    0.6275]; % #6E4EA0
%                       [0.4667    0.2235    0.5882]; % #773A96
%                       [0.6627    0.1333    0.4745];  % #A92179
%                       [0.9333    0.2000    0.1412];  % #EE3324
%                       [0.9569    0.4863    0.1255];  % #F47C20
%                       [0.9843    0.6588    0.1020]; % #FBA919};
%                     };
pro_num = size(protein_clim);

figure();
hold on;
for i = 1:pro_num(2)
% scatter(protein_clim(:,i), protein_nlim(:,i),'MarkerEdgeColor',protein_colors{i}, 'MarkerFaceColor', protein_colors{i}, 'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.2,'LineWidth', 2); 
plot(protein_clim(:,i), protein_nlim(:,i), protein_shape{i}, 'Color', protein_colors{i}, 'MarkerFaceColor',protein_colors{i}, 'MarkerSize', 8, 'LineWidth', 2);
end
plot([0.01 100], [0.01 100], 'k-', 'LineWidth', 0.5)
plot([0.01 50], [0.02 100], 'k:', 'LineWidth', 0.5)
plot([0.01 200], [0.005 100], 'k:', 'LineWidth', 0.5)
xlim([0.01 100])
ylim([0.01 100])
xticks([0.01 1 100])
yticks([0.01 1 100])
set(gca, 'XScale','log')
set(gca, 'YScale','log')
set(gca, 'XMinorTick', 'off')
set(gca, 'YMinorTick', 'off')
xlabel('protein sector % carbon-lim')
ylabel('protein sector % nitrogen-lim')
axis square; 
box on;
hold off;



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



function error = fun_mcmc_loss(params, data, para_f)
    load_simulation_setting;
    data_unit_conv; 
    plot_num; 
    par = read_par(params, para_f); 
    
    [WT_y_steady_clim,   WT_met_reac_steady_clim, ...
     WT_y_steady_cnlim1, WT_met_reac_steady_cnlim1, ...
     WT_y_steady_cnlim2, WT_met_reac_steady_cnlim2, ...
     WT_y_steady_nlim,   WT_met_reac_steady_nlim,   D_chem] = func_chemostat_simulation(par);
    
    error = func_chemo_error(WT_y_steady_clim,   WT_met_reac_steady_clim, ...
                              WT_y_steady_cnlim1, WT_met_reac_steady_cnlim1, ...
                              WT_y_steady_cnlim2, WT_met_reac_steady_cnlim2, ...
                              WT_y_steady_nlim,   WT_met_reac_steady_nlim, ...
                              num_y, num_flux, num_prot, ...
                              pt_uM_to_gPerL, lp_uM_to_gPerL, lipid_sc_c,...
                              D_chem, par);
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



