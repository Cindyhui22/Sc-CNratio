% chemostat simu function - adpated from Sc_paper_figures_chemostat.m
function [WT_y_steady_clim,   WT_met_reac_steady_clim, ...
          WT_y_steady_cnlim1, WT_met_reac_steady_cnlim1, ...
          WT_y_steady_cnlim2, WT_met_reac_steady_cnlim2, ...
          WT_y_steady_nlim,   WT_met_reac_steady_nlim, D_chem] = func_chemostat_simulation(par_sto) 


%%% parameters and simulation setup
par_sc = par_sto; 
gy_ratio = 1;
amino_acid_ex = 0;
batch   = 0; 
% unit conversions
load_simulation_setting;
data_unit_conv; 
plot_num; 
plt_labels = plot_labels; 
plt_colors = plot_colors;
cell_fac = par_sc.V_c*par_sc.rho_cell; % cells in g/L
media_type  = 'clim_nlim'; % clim: green, nlim: blue
yeast_type  = 'sc'; 

D_chem = linspace(0.05, 0.42, 30);
D = D_chem; 
chem_t_start = 0; 
chem_t_final = 6000;   

% initial conditions - generic, minimal media 
aa_ex       = 0; 
gl          = 100; 
eh          = 0; 
cells       = 1.0E10 * cell_fac; 
nh_ex       = 0.001 * nh4_gPerL_to_uM;
mutant      = 'none'; % WT

snf1_vals = 0;
jgy_vals  = 0; 
mgl       = 0; 

%%% Clim chemostat SC
par = par_sc; 
par.gl_in   = 10 * gl_gPerL_to_uM; 
par.nh4_in  = 62.5 * nh4_gPerL_to_uM;  
[WT_y_steady_clim, WT_met_reac_steady_clim, ~]...
= function_chemostat_all_D(par, D_chem, t_chem_start, t_chem_final, aa_ex,...
nh_ex, gl, eh, cells, mutant, snf1_vals, jgy_vals, mgl, ...
plt_labels, num_y, num_flux, num_prot);



%%% CNlim1 chemostat SC
par = par_sc; 
par.gl_in   = 7.5/1.5* gl_gPerL_to_uM;
par.nh4_in  = 1.875/1.5* nh4_gPerL_to_uM;
[WT_y_steady_cnlim1, WT_met_reac_steady_cnlim1, ~]...
= function_chemostat_all_D(par, D_chem, t_chem_start, t_chem_final, aa_ex,...
nh_ex, gl, eh, cells, mutant, snf1_vals, jgy_vals, mgl, ...
plt_labels, num_y, num_flux, num_prot);




%%% CNlim2 chemostat SC
par = par_sc; 
par.gl_in   = 7.5/1.5* gl_gPerL_to_uM;
par.nh4_in  = 1.25/1.5 * nh4_gPerL_to_uM;
[WT_y_steady_cnlim2, WT_met_reac_steady_cnlim2,~]...
= function_chemostat_all_D(par, D_chem, t_chem_start, t_chem_final, aa_ex,...
nh_ex, gl, eh, cells, mutant, snf1_vals, jgy_vals, mgl, ...
plt_labels, num_y, num_flux, num_prot);



%%% Nlim chemostat SC
par = par_sc; 
par.gl_in   = 7.5 * gl_gPerL_to_uM; % inlet glucose concentration
par.nh4_in  = 3.75 * nh4_gPerL_to_uM/15; %_% 5 * nh4_gPerL_to_uM/10; % inlet NH4 concentration
[WT_y_steady_nlim, WT_met_reac_steady_nlim, ~]...
= function_chemostat_all_D(par, D_chem, t_chem_start, t_chem_final, aa_ex,...
nh_ex, gl, eh, cells, mutant, snf1_vals, jgy_vals, mgl, ...
plt_labels, num_y, num_flux, num_prot);


end