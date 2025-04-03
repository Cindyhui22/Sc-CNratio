clc; close all; clear;  
addpath('../')
addpath(genpath('./'))


% read parameters 
par = read_parametersv2; % reference parameter
load_simulation_setting;
data_unit_conv; 
plot_num; 
plt_labels = plot_labels; 
plt_colors = plot_colors;

keySet_l = {'l_r','l_z','l_gy','l_fe','l_gn','l_mt','l_as','l_at', 'l_lp', 'l_lo'};
valueSet_l = 1: numel(keySet_l);
par_l    =  containers.Map(keySet_l,valueSet_l); 

keySet_ktj = {'K_t_r','K_t_z','K_t_gy','K_t_fe','K_t_gn','K_t_mt','K_t_as','K_t_at','K_t_lp','K_t_lo'};
valueSet_ktj = 1: numel(keySet_ktj);
par_k_tj =  containers.Map(keySet_ktj,valueSet_ktj);


% get parameters names
parName         = fieldnames(par);
% parName(find(strcmp(parName,'l')))    = []; % delete par.l not perturb par.l going to perturb each element in par.l
% parName(find(strcmp(parName,'k_tj'))) = [];


parName(find(strcmp(parName,'K_3AT'))) = []; % not a base model parameters
parName(find(strcmp(parName,'c_3AT'))) = []; % not a base model parameters
parName(find(strcmp(parName,'K_mgl'))) = []; % not a base model parameters

% -------------------------------------------------------------------------------------------------------  
%  Plotting stuff: parameter colors by parameter category 
%  ------------------------------------------------------------------------------------------------------- 

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
% writetable(par_all_type_table,'par_value_excle.xlsx')
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

% parameters perturbing ratio 
delta = 1;

% perturb every parameter by delta
par_sto = repmat(par, 1, length(parName)); % 
for i = 1:length(parName)
if  sum(strcmp(parName(i),keySet_l)) == 1
par_sto(i).l(par_l(parName{i})) = max(par_sto(i).l(par_l(parName{i})) + ...
delta * par_sto(i).l(par_l(parName{i})), 0.1);
end
if sum(strcmp(parName(i),keySet_ktj)) == 1
par_sto(i).k_tj(par_k_tj(parName{i})) = max(par_sto(i).k_tj(par_k_tj(parName{i})) + ...
delta * par_sto(i).k_tj(par_k_tj(parName{i})), 0.1);
end 
par_sto(i).(parName{i})= max(par_sto(i).(parName{i})+ delta * par_sto(i).(parName{i}), 0.1);  % perturbed parameters
end

none_indices = strcmp(parCategory, 'none');
if any(none_indices)
    none_params = parName(none_indices);  % Get the names of the parameters
    warning(['The following parameters have no category assigned and will be removed: ', strjoin(none_params, ', ')]);

    % Remove the corresponding entries from parCategory, par_sto, and parName
    parCategory(none_indices) = [];
    par_sto(none_indices) = [];
    parName(none_indices) = [];
end


% perturbed parameter colors list 
% set up colormap based on parameter type
color_map = zeros(length(parName),3);

par_cat = {'ori_nonfree','ori_free','ori_adapt_nonfree','ori_adapt_free', 'new_nonfree', 'new_free'};

% par_cat_cls =  ['#868686','#0096c5','#b6b6b6','6fccec','#b86f1c','#ffc21f']; 
par_cat_cls = [ 0.5255    0.5255    0.5255  ;% #868686
                0.0000    0.5882    0.7725  ;% #0096c5
                0.7137    0.7137    0.7137  ;% #b6b6b6
                0.4353    0.8039    0.9255  ;% #6fccec
                0.7216    0.4353    0.1137  ;% #b86f1c
                1.0000    0.7608    0.1216  ;% #ffc21f
               ];


par_cat_size = zeros(numel(par_cat),1); 
for i = 1: numel(par_cat)
indx = find(strcmp(parCategory ,par_cat{i}));
color_map(indx,:) = repmat(par_cat_cls(i,:), length(indx), 1);
par_cat_size(i) = numel(indx); 
end 
%% -------------------------------------------------------------------------------------------------------  
%  Read in & Organize Exp for comparison;
%  ------------------------------------------------------------------------------------------------------- 
func_exp_load;

%% -------------------------------------------------------------------------------------------------------  
%  Simulate with reference parameters
%  ------------------------------------------------------------------------------------------------------- 
%
% run chemostat and calculate sim-exp error; with unperturbed reference parameter; use as reference later. 
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
%
% run batch and calculate sim-exp error; with unperturbed reference parameter; use as reference later. 
[clim_WT_batch_min_t, clim_WT_batch_min_y, ...
 clim_WT_batch_YPD_t, clim_WT_batch_YPD_y] = func_batch_simulation(par);

error_batch = func_batch_error_v2(clim_WT_batch_min_t, clim_WT_batch_min_y, ...
                                     clim_WT_batch_YPD_t, clim_WT_batch_YPD_y,...
                                     glucose_climb, ethanol_climb,  cell_climb, prot_sec_climb,...
                                     par, num_y, num_flux, num_prot);


error_ref = error_chem + error_batch; 

%% -------------------------------------------------------------------------------------------------------  
%  Simulate with purturbed parameters and calculate sim-exp error
%  ------------------------------------------------------------------------------------------------------- 

% n = 3; %run a small sample
n = length(parName);

error_preturb     = zeros(n,1); % 

% run simulation and calculate sim-exp error with perturbed parameter  
% chemostat
for i = 1:n
[WT_y_steady_clim,   WT_met_reac_steady_clim, ...
 WT_y_steady_cnlim1, WT_met_reac_steady_cnlim1, ...
 WT_y_steady_cnlim2, WT_met_reac_steady_cnlim2, ...
 WT_y_steady_nlim,   WT_met_reac_steady_nlim,   D_chem] = func_chemostat_simulation(par_sto(i));

error_prechem = func_chemo_error_v2(WT_y_steady_clim,   WT_met_reac_steady_clim, ...
                                WT_y_steady_cnlim1, WT_met_reac_steady_cnlim1, ...
                                WT_y_steady_cnlim2, WT_met_reac_steady_cnlim2, ...
                                WT_y_steady_nlim,   WT_met_reac_steady_nlim, ...
                                D_clim, Jgy_clim, Jeh_clim, cell_clim, precursor_clim, aa_clim, protein_clim, lipid_clim, prot_sec_clim, ...
                                D_cnlim1, Jgy_cnlim1, Jeh_cnlim1, nh4_cnlim1, protein_cnlim1, prot_sec_cnlim1, ...
                                D_cnlim2, Jgy_cnlim2, Jeh_cnlim2, cell_cnlim2, precursor_cnlim2, aa_cnlim2, nh4_cnlim2, protein_cnlim2, prot_sec_cnlim2, ...
                                D_nlimf, Jgy_nlimf, Jeh_nlimf, cell_nlimf, precursor_nlimf, aa_nlimf, protein_nlimf, lipid_nlimf, prot_sec_nlimf, ...
                                num_y, num_flux, num_prot, ...
                                pt_uM_to_gPerL, lp_uM_to_gPerL, lipid_sc_c,...
                                D_chem, par_sto(i));

[clim_WT_batch_min_t, clim_WT_batch_min_y, ...
 clim_WT_batch_YPD_t, clim_WT_batch_YPD_y] = func_batch_simulation(par_sto(i));

error_prebatch = func_batch_error_v2(clim_WT_batch_min_t, clim_WT_batch_min_y, ...
                                     clim_WT_batch_YPD_t, clim_WT_batch_YPD_y,...
                                     glucose_climb, ethanol_climb,  cell_climb, prot_sec_climb,...
                                     par_sto(i), num_y, num_flux, num_prot);
 

error_preturb(i) = error_prechem + error_prebatch; 
end

% -------------------------------------------------------------------------------------------------------  
%%  Calculate sensitivity and sort the sensitivity 
%  ------------------------------------------------------------------------------------------------------- 
% Local sensitivity (defined as the first-order partial derivative of the averaged relative squared
% error of predictions)  (same as Chen's SI equation 63)
S = (abs(log(error_preturb)-log(error_ref)))/log(1+delta);
[sorted_S, ind_order] = sort(S,'descend');

%%
% -------------------------------------------------------------------------------------------------------  
%  plotting and save plot
%  ------------------------------------------------------------------------------------------------------- 

% find the corresponding parameter name and color after sort local sensitivity from high to low
sorted_parName       = parName(ind_order);
sorted_color_map     = color_map(ind_order,:);

% bar plot setup
bar_width       = 0.4;
bar_num         = 26; % bar per plot
subplot_num     = 6;  % subplot per figure
% Y_ticks         = [10^-7,10^-4,10^-1]; 

% plot
count = subplot_num + 1;  
for i = 1: floor(n/bar_num)
if mod(count, subplot_num) == 1 % subplot_num of subplot per figure
figure; 
count = count - subplot_num; 
end 

subplot(subplot_num,1,count)
hold on; 
for k = 1:bar_num % bar_num of bar per subplot
indx = (i-1) * bar_num + k; 
h = bar(k,sorted_S(indx),bar_width);
set(h,'FaceColor',sorted_color_map(indx,:)); % change color for each bar
end
hold off;
% ylim([min(Y_ticks), max(Y_ticks)*100])
ylim([10^-4, 10^1])
set(gca,'Ytick',[10^-4, 10^-2, 10^0])
set(gca,'Xtick',1:bar_num)
xticklabels(sorted_parName(indx - bar_num + 1:indx))
set(gca,'TickLabelInterpreter','none')
set(gca,'fontsize',12)
set(gca,'fontname','Arial')    
% set(gca,'Ytick',Y_ticks)
set(gca,'YScale','log')
set(gca,'YMinorTick','off')
set(gca,'XTickLabelRotation',45)
box on;
xlim([0 bar_num+1])
count = count + 1;
end 

subplot(subplot_num,1,count)
hold on; 
lastbar_num = (length(parName)- indx -1);
for k = 1: lastbar_num % bar_num of bar per subplot
indx = indx + 1; 
h = bar(k,sorted_S(indx),bar_width);
set(h,'FaceColor',sorted_color_map(indx,:)); % change color for each bar
end
hold off;
ylim([10^-4, 10^1])
set(gca,'Ytick',[10^-4, 10^-2, 10^0])
set(gca,'Xtick',1:lastbar_num)
xticklabels(sorted_parName(indx - lastbar_num + 1:indx+1))
set(gca,'TickLabelInterpreter','none')
set(gca,'fontsize',12)
set(gca,'fontname','Arial')    
set(gca,'YScale','log')
set(gca,'YMinorTick','off')
set(gca,'XTickLabelRotation',45)
box on;
xlim([0 bar_num+1])
set(gcf,'Position',[1000 432 488 905])

%%
% save data
% save('local_sensi_v2.mat', "parName","sorted_parName", "S", "sorted_S", "ind_order", "error_preturb","error_ref","delta", "par_sto");
%
