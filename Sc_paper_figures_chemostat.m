%% SC paper figure 1, 2, 3.


% GREEN is ref par1 
clc; close all; clear; 
addpath('../')
addpath(genpath('./'))

par = read_parametersv2;

% unit conversions
cell_fac = par.V_c*par.rho_cell; % cells in g/L
% cell_fac = 1; % cells in number
media_type  = 'clim_nlim'; % clim: green, nlim: blue
yeast_type  = 'sc'; 


data_unit_conv; 
plot_num; 
plt_labels = plot_labels;
plt_colors = plot_colors; 
load_simulation_setting;
% set(0,'DefaultLineLineWidth',2);



%% parameters and simulation setup

par_sc = read_parametersv2; %_%

gy_ratio = 1;%_% 1.4; 
amino_acid_ex = 0;%_% 0*1*2*10^5;
batch   = 0; % 1: plot batch 



%% chemostat simulation - WT, minimal media --------------------------------
% chemostat simulation ----------------------------------------------------
% chemostat experiment conditions
D_chem = linspace(0.05, 0.42, 30);
D = D_chem; 

% experiment conditions
% batch simulation --------------------------------------------------------
chem_t_start = 0; 
chem_t_final = 6000;   


snf1_vals = 0; % does not do anything - just need some value for 'module_signaling.m' 
jgy_vals  = 0; % does not do anything - just need some value for 'module_signaling.m' 
mgl       = 0; % does not do anything - just need some value for 'module_met.m' 

% initial conditions - generic, minimal media 
aa_ex       = 0; 
gl          = 100; 
eh          = 0; 
cells       = 1.0E10 * cell_fac; 
nh_ex       = 0.001 * nh4_gPerL_to_uM;
mutant      = 'none'; % WT

%% Clim chemostat SC
%_% C/N ratio = 2.2*10/62.5 = 0.352 
par = par_sc; 
%_% par.k_gy    = par.k_gy   * gy_ratio; 
par.gl_in   = 10 * gl_gPerL_to_uM;  % inlet glucose concentration
par.nh4_in  = 62.5 * nh4_gPerL_to_uM;  %_% 5 * nh4_gPerL_to_uM; % inlet NH4 concentration

[WT_y_steady_clim, WT_met_reac_steady_clim, WT_other_met_reac_steady_clim, WT_prot_syn_steady_clim,...
WT_sig_steady_clim, WT_rib_steady_clim, WT_g_rate_steady_clim, WT_ribo_rate_steady_clim]...
= function_chemostat_all_D(par, D_chem, t_chem_start, t_chem_final, aa_ex,...
nh_ex, gl, eh, cells, mutant, snf1_vals, jgy_vals, mgl, ...
plt_labels, num_y, num_flux, num_prot);



disp('WT chemostat simulation running - glucose limited condition...')


%% CNlim1 chemostat SC
%_% C/N ratio = 2.2*10/(5/3.5) = 15.4; 2.2*10/(5/2)=8.8
par = par_sc; 
%_% par.k_gy    = par.k_gy   * gy_ratio; % 
par.gl_in   = 7.5/1.5* gl_gPerL_to_uM;%_%10 * gl_gPerL_to_uM;  % inlet glucose concentration
%par.gl_in   = 7.5 * gl_gPerL_to_uM;  % inlet glucose concentration
par.nh4_in  = 1.875/1.5* nh4_gPerL_to_uM;%_%5 * nh4_gPerL_to_uM/2; %_% 5 * nh4_gPerL_to_uM/3.5; % inlet NH4 concentration

[WT_y_steady_cnlim1, WT_met_reac_steady_cnlim1, WT_other_met_reac_steady_cnlim1, WT_prot_syn_steady_cnlim1,...
WT_sig_steady_cnlim1, WT_rib_steady_cnlim1, WT_g_rate_steady_cnlim1, WT_ribo_rate_steady_cnlim1]...
= function_chemostat_all_D(par, D_chem, t_chem_start, t_chem_final, aa_ex,...
nh_ex, gl, eh, cells, mutant, snf1_vals, jgy_vals, mgl, ...
plt_labels, num_y, num_flux, num_prot);



disp('WT chemostat simulation running - glucose limited condition...')



%% CNlim2 chemostat SC
%_% C/N ratio = 2.2*10/(5/5) = 22; 2.2*10/(5/3)=13.2
par = par_sc; 
%_% par.k_gy    = par.k_gy   * gy_ratio; % 
par.gl_in   = 7.5/1.5* gl_gPerL_to_uM;%_% 10 * gl_gPerL_to_uM;  % inlet glucose concentration
%par.gl_in   = 7.5 * gl_gPerL_to_uM;  % inlet glucose concentration

par.nh4_in  = 1.25/1.5 * nh4_gPerL_to_uM;%_% 5 * nh4_gPerL_to_uM/3; %_% 5 * nh4_gPerL_to_uM/5; % inlet NH4 concentration
%par.nh4_in  = 5 * nh4_gPerL_to_uM/5; % inlet NH4 concentration

[WT_y_steady_cnlim2, WT_met_reac_steady_cnlim2, WT_other_met_reac_steady_cnlim2, WT_prot_syn_steady_cnlim2,...
WT_sig_steady_cnlim2, WT_rib_steady_cnlim2, WT_g_rate_steady_cnlim2, WT_ribo_rate_steady_cnlim2]...
= function_chemostat_all_D(par, D_chem, t_chem_start, t_chem_final, aa_ex,...
nh_ex, gl, eh, cells, mutant, snf1_vals, jgy_vals, mgl, ...
plt_labels, num_y, num_flux, num_prot);



disp('WT chemostat simulation running - glucose limited condition...')


%% Nlim chemostat SC
%_% C/N ratio = 2.2*10/(5/10) = 44; 2.2*10/(5/15)=66
par = par_sc; 
%_% par.k_gy    = par.k_gy   * gy_ratio; % 
par.gl_in   = 7.5 * gl_gPerL_to_uM; % inlet glucose concentration
par.nh4_in  = 3.75 * nh4_gPerL_to_uM/15; %_% 5 * nh4_gPerL_to_uM/10; % inlet NH4 concentration
[WT_y_steady_nlim, WT_met_reac_steady_nlim, WT_other_met_reac_steady_nlim, WT_prot_syn_steady_nlim,...
WT_sig_steady_nlim, WT_rib_steady_nlim, WT_g_rate_steady_nlim, WT_ribo_rate_steady_nlim]...
= function_chemostat_all_D(par, D_chem, t_chem_start, t_chem_final, aa_ex,...
nh_ex, gl, eh, cells, mutant, snf1_vals, jgy_vals, mgl, ...
plt_labels, num_y, num_flux, num_prot);

disp('WT chemostat simulation running - NH4 limited condition...')





%% Fig2: plot Sc Clim vs Nlim chemo (metabolic pathway)

figure_output_format  = ...
{'NH4'    , 'glucose' , 'cells'     ,'protein' , 'snf1'  ; 
'aa_in'   , 'ethanol' ,  'precursor', 'lipid'  , 'tor'  ;
'NA'      , 'NA'      ,  'NA'       , 'NA'     , 'JehJgy',   ;
'NA'      , 'NA'      , 'NA'        , 'NA'     , 'NA'   ; 
'NA'      , 'NA'      , 'NA'        , 'NA'     ,  'NA'}'; 


figure_position = [101 67 668 738];

paper_plot_sc_met_lerp(D, D_crit_simu, ...
WT_y_steady_clim, WT_y_steady_nlim, WT_y_steady_cnlim1, WT_y_steady_cnlim2, ...
WT_met_reac_steady_clim, WT_met_reac_steady_nlim, WT_met_reac_steady_cnlim1, WT_met_reac_steady_cnlim2,...
WT_sig_steady_clim, WT_sig_steady_nlim, WT_sig_steady_cnlim1, WT_sig_steady_cnlim2, ...
par, num_y, num_flux, num_prot, figure_output_format, figure_position)




%% Fig3: plot Sc Clim vs Nlim chemo (gene expression pathway)


figure_output_format  = ...
{
'R'       , 'Egy'     ,  'Efe'   , 'REZ'     , 'NA',   ;
'Emt'     , 'Eas'     , 'Egn'    , 'REZ_nlim', 'NA'   ; 
'Elp'     , 'Elo'     , 'Z'      , 'NA'      , 'NA'; 
'Eat'      , 'NA'      , 'NA'     ,'NA'       , 'NA'  ; 
'NA'      , 'NA'      ,  'NA'    , 'NA'      , 'NA'  ;}'; 
figure_position = [101 67 668 738];

paper_plot_sc_protein_lerp(D, D_crit_simu, ...
WT_y_steady_clim, WT_y_steady_nlim, WT_y_steady_cnlim1, WT_y_steady_cnlim2, ...
WT_met_reac_steady_clim, WT_met_reac_steady_nlim, WT_met_reac_steady_cnlim1, WT_met_reac_steady_cnlim2,...
WT_sig_steady_clim, WT_sig_steady_nlim, WT_sig_steady_cnlim1, WT_sig_steady_cnlim2, ...
par, num_y, num_flux, num_prot, figure_output_format, figure_position)




%% Fig4: plot Sc Clim vs Nlim chemo (ATP)

figure_output_format = ''; 
paper_plot_sc_atp_lerp(D, D_crit_simu, ...
WT_y_steady_clim, WT_y_steady_nlim, WT_y_steady_cnlim1, WT_y_steady_cnlim2, ...
WT_met_reac_steady_clim, WT_met_reac_steady_nlim, WT_met_reac_steady_cnlim1, WT_met_reac_steady_cnlim2,...
WT_sig_steady_clim, WT_sig_steady_nlim, WT_sig_steady_cnlim1, WT_sig_steady_cnlim2, ...
par, num_y, num_flux, num_prot, figure_output_format, figure_position)




%% plot Sc Clim vs Nlim chemo (SI)

figure_output_format = {'aa_in', 'precursor', 'NA', 'NA'; 
'atp', 'NA', 'NA', 'NA'; 
'REZ', 'REZ_nlim', 'NA', 'NA';
'Egn', 'R', 'NA', 'NA'; }';

figure_position = [440 248 726 550];

paper_plot_sc_SI_lerp(D, D_crit_simu, ...
WT_y_steady_clim, WT_y_steady_nlim, WT_y_steady_cnlim1, WT_y_steady_cnlim2, ...
WT_met_reac_steady_clim, WT_met_reac_steady_nlim, WT_met_reac_steady_cnlim1, WT_met_reac_steady_cnlim2,...
WT_sig_steady_clim, WT_sig_steady_nlim, WT_sig_steady_cnlim1, WT_sig_steady_cnlim2, ...
par, num_y, num_flux, num_prot, figure_output_format, figure_position)


%% Fig5: plot Sc Clim vs Nlim batch 

