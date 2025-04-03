
% line for critical dilution rate 
% D_crit_simu = {'set_by_hand', 0.19}; % manually set critical dilution rate in simulation (read by eye when running unnormalized simulation)
D_crit_simu = {'auto_calculation'}; % automatically determine dilution rate as where the highest ethanol prod. fold change happens in simulation 

% plotting stuff
x_label_batch = {'time (h)'};
x_label_chem  = {'dilution rate (h^{-1})'};



par_type   = 'change_batch_chemo';
mutant     = 'none';
%par_type = 'same_all_paper';
%par_type = 'diff_all_paper';
%snf1_vals = 0; % does not do anything - just need some value for 'module_signaling.m' 
%jgy_vals  = 0; % does not do anything - just need some value for 'module_signaling.m' 
%mgl       = 0; % does not do anything - just need some value for 'module_met.m' 

%% experiment conditions

% chemost1at simulation ----------------------------------------------------
%D_chem = linspace(0.027, 0.43, 50*2); % 
D_chem = linspace(0.027, 0.41, 12); % start with lowest D in exp.

t_chem_start = 0; 
t_chem_final = 6000;   

% hxt and snf1 mutant simulations -----------------------------------------

D_hxt_eh = 0.1; % dilution rate for ethanol chemostat
D_hxt_gl = 0.1; % dilution rate for glucose consumption phase

D_snf1 = linspace(0.05, 0.35, 2);

% batch simulation --------------------------------------------------------

D_batch = 0;
batch_t_start = 0; 
batch_t_final =  35;
%batch_t_final = 35; 20; %35;   


