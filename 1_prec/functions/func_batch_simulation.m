function [clim_WT_batch_min_t, clim_WT_batch_min_y, ...
          clim_WT_batch_YPD_t, clim_WT_batch_YPD_y] = func_batch_simulation(par_sto)     

load_simulation_setting
data_unit_conv; 
plot_num; 
% disp(' ')
par = par_sto;
cell_fac = par.V_c*par.rho_cell; % cells in g/L
batch_t_start = 0;  
batch_t_final = 40;  
%% WT batch simulations --------------------------------------------------
%
% disp(' ')
% disp('------------------------- WT simulations -------------------------')

% simulation conditions
mutant     = 'none'; % WT

snf1_vals = 0; % does not do anything - just need some value for 'module_signaling.m' 
jgy_vals  = 0; % does not do anything - just need some value for 'module_signaling.m' 
mgl       = 0; % does not do anything - just need some value for 'module_met.m' 

% bartolomeo --------------------------------------------------------------
% disp('WT batch simulation running - bartolomeo...')

% initial conditions - bartolomeo
clim_gl    = 23 * gl_gPerL_to_uM; 
clim_eh    = 0; 
clim_cells = 6.89E9*cell_fac; 
clim_nh4   = 5*nh4_gPerL_to_uM*3; 
clim_aa_ex = 0; 
% run simulation 
[clim_WT_batch_min_t, clim_WT_batch_min_y, ~] = function_batch(par, batch_t_start, batch_t_final, clim_aa_ex, clim_nh4, clim_gl, clim_eh, clim_cells, mutant, snf1_vals, jgy_vals, mgl);


% murphy ------------------------------------------------------------------
% disp('WT batch simulation running - murphy...')

%initial conditions - complex 
clim_aa_ex = 2*10^6; 
clim_gl    = 18 * gl_gPerL_to_uM; 
clim_eh    = 0; 
clim_cells = 3.76E09*cell_fac; %7.52E8; 
clim_nh4   = 5 * nh4_gPerL_to_uM;
% run simulation 
[clim_WT_batch_YPD_t, clim_WT_batch_YPD_y,~] = function_batch(par, batch_t_start, batch_t_final, clim_aa_ex, clim_nh4, clim_gl, clim_eh, clim_cells, mutant, snf1_vals, jgy_vals, mgl);

end