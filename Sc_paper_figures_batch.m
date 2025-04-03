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



%% change all


par_sc = read_parametersv2; %_%


% green ref par1 
% blue perturbed par2

gy_ratio = 1;%_% 1.4; 
amino_acid_ex = 0; %_% 0*1*2*10^5;
batch   = 1 ; % 1: plot batch 
chemo   = 1; 


batch_t_start = 0;  
batch_t_final = 60;  

% chemostat simulation ----------------------------------------------------
% chemostat experiment conditions
D_chem = linspace(0.05, 0.42, 30);
D = D_chem; 




%% experiment conditions

% batch simulation --------------------------------------------------------
chem_t_start = 0; 
chem_t_final = 6000;   

%% WT batch and chemostat simulations -------------------------------------


disp(' ')
disp('------------------------- WT batch simulations -------------------------')

% simulation conditions
mutant     = 'none'; % WT

snf1_vals = 0; % does not do anything - just need some value for 'module_signaling.m' 
jgy_vals  = 0; % does not do anything - just need some value for 'module_signaling.m' 
mgl       = 0; % does not do anything - just need some value for 'module_met.m' 


%% Clim batch SC
%
par        = par_sc; 
% initial conditions - bartolomeo
clim_gl    = 23 * gl_gPerL_to_uM; 
clim_eh    = 0; 
clim_cells = 6.89E9*cell_fac; 
clim_nh4   = 5*nh4_gPerL_to_uM*3; 
clim_aa_ex = 0; 

disp('WT batch simulation running - Clim...')

% run simulation 
[clim_WT_batch_min_t, clim_WT_batch_min_y, clim_WT_batch_min_met_reac, clim_WT_batch_min_other_met_reac, clim_WT_batch_min_sig, clim_WT_batch_min_prot_syn, clim_WT_batch_min_rib, clim_WT_batch_min_g_rate, clim_WT_batch_min_ribo_rate] = function_batch(par, batch_t_start, batch_t_final, clim_aa_ex, clim_nh4, clim_gl, clim_eh, clim_cells, mutant, snf1_vals, jgy_vals, mgl);



%% Nlim batch SC
par       = par_sc; 
% initial conditions - zampar
nlim_gl    = 23 * gl_gPerL_to_uM; 
nlim_eh    = 0; 
nlim_cells = 6.89E9*cell_fac; 
nlim_nh4   = 5*nh4_gPerL_to_uM/5;
nlim_aa_ex = 0; 

disp('WT batch simulation running - Nlim...')

% run simulation 
[nlim_WT_batch_min_t, nlim_WT_batch_min_y, nlim_WT_batch_min_met_reac, nlim_WT_batch_min_other_met_reac, nlim_WT_batch_min_sig, nlim_WT_batch_min_prot_syn, nlim_WT_batch_min_rib, nlim_WT_batch_min_g_rate, nlim_WT_batch_min_ribo_rate] = function_batch(par, batch_t_start, batch_t_final, nlim_aa_ex, nlim_nh4, nlim_gl, nlim_eh, nlim_cells, mutant, snf1_vals, jgy_vals, mgl);



%% complex SC
% complex ------------------------------------------------------------------

par       = par_sc; 
%initial conditions - complex 
clim_aa_ex = 2*10^6; 
clim_gl    = 18 * gl_gPerL_to_uM; 
clim_eh    = 0; 
clim_cells = 3.76E09*cell_fac; %7.52E8; 
clim_nh4   = 5 * nh4_gPerL_to_uM;

disp('WT batch simulation running - Complex...')

% run simulation 
[clim_WT_batch_YPD_t, clim_WT_batch_YPD_y, clim_WT_batch_YPD_met_reac, clim_WT_batch_YPD_other_met_reac, clim_WT_batch_YPD_sig, clim_WT_batch_YPD_prot_syn, clim_WT_batch_YPD_rib, clim_WT_batch_YPD_g_rate, clim_WT_batch_YPD_ribo_rate] = function_batch(par, batch_t_start, batch_t_final, clim_aa_ex, clim_nh4, clim_gl, clim_eh, clim_cells, mutant, snf1_vals, jgy_vals, mgl);

%}


%%  ---------------------  plot batch

%% Sc batch clim vs nlim
clc; close all
% figure_output_format = ...
%   {'glucose', 'precursor', 'ethanol', 'aa_in' , 'atp'  ,  'cells'   , 'NH4'   ; 
%    'protein',  'lipid'   , 'snf1'   , 'tor'   , 'R'    ,  'Z'       , 'Egy'   ;
%     'Efe'   , 'Eas'      ,  'Emt'   , 'Egn'   , 'Elp'  , 'Elo'      ,  ''     ;
%     'NA'    , 'NA'       ,  'REZ'   , 'NA'   , 'NA',  'NA'       ,  'NA'   ; }'; 
% figure_position = [143 571 993 490];

% figure_output_format  = ...
%  {'NH4'    , 'glucose' , 'cells'     ,'snf1' , 'protein'  ; 
%  'aa_in'   , 'ethanol' ,  'precursor', 'tor'  , 'lipid'  ;
%  'R'    ,  'Z'       , 'Egy'      , 'Efe'     , 'REZ',   ;
%  'Eas'      ,  'Emt'   , 'Egn'   , 'Elp'  , 'Elo'   ; 
%   'NA' , 'atp'     , 'NA'        , 'NA'     ,  'NA'}'; 

figure_output_format  = ...
{'NH4'    , 'aa_in' , 'protein'     ,'cells' , 'snf1'  ; 
'glucose'   , 'ethanol' ,  'lipid', 'tor'  , 'precursor'  ;
'R'    ,  'Z'       , 'Egy'      , 'Efe'     , 'REZ',   ;
'Eas'      ,  'Emt'   , 'Egn'   , 'Elp'  , 'Elo'   ; 
'NA' , 'atp'     , 'NA'        , 'NA'     ,  'NA'}'; 

figure_position = [101 67 668 738];


paper_plot_sc_batch_clim_v_nlim(clim_WT_batch_min_t, nlim_WT_batch_min_t,  clim_WT_batch_YPD_t,...
clim_WT_batch_min_y, nlim_WT_batch_min_y, clim_WT_batch_YPD_y,  ...
clim_WT_batch_min_sig, nlim_WT_batch_min_sig, clim_WT_batch_YPD_sig, ...
clim_WT_batch_min_met_reac, nlim_WT_batch_min_met_reac, clim_WT_batch_YPD_met_reac, ...
par, num_y, num_prot, num_flux, figure_output_format, figure_position) 

