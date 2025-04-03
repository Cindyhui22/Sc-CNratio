function error = function_chemo_sim_exp_error(D_chem, WT_y_steady_10gL, WT_met_reac_steady_10gL, WT_sig_steady_10gL, D_crit_simu, par, num_y, num_flux, num_prot)
%
% This function calculate the error between experiment and simulation for chemostat
% Relative Error = (abs(exp data/sim data - 1))^2
% Absolute Error = (abs(exp data - sim data))^2
%
% Example: 
% 
% 
% 
%
%% -------------------------------------------------------------------------------------------------------  
% experiment, simulation matrix index 
% 
num_exp.Jgl   = 1;
num_exp.Jeh   = 2; 
num_exp.cells = 3;
num_exp.aa    = 4;
num_exp.ae    = 5; 
num_exp.r     = 6; 
num_exp.z     = 7;
num_exp.gy    = 8; 
num_exp.fe    = 9; 
num_exp.gn    = 10; 
num_exp.mt    = 11; 
num_exp.as    = 12; 
num_exp.at    = 13; 


%% -------------------------------------------------------------------------------------------------------  
%  load exp data
%  ------------------------------------------------------------------------------------------------------- 
%% load exp data
data_chemostat; 
data_hackett_prot; 
data_4_prec;
data_different_carbon; 
data_xia_new; 


%% -------------------------------------------------------------------------------------------------------  
%  normalize precursor, amino acid, atp, protein with respect to first data
%  point
%  ------------------------------------------------------------------------------------------------------- 
norm_indx = 1;

% normalize precursor, amino acid, atp to fold change
WT_y_steady_10gL(:,[num_y.pc, num_y.aa_in, num_y.ae]) = WT_y_steady_10gL(:,[num_y.pc,num_y.aa_in, num_y.ae]) ./ WT_y_steady_10gL(norm_indx,[num_y.pc,num_y.aa_in, num_y.ae]);

kumar_lg.glutamine = kumar_lg.glutamine ./ kumar_lg.glutamine(1);
kumar_hg.glutamine = kumar_hg.glutamine ./ kumar_hg.glutamine(1);
hackett.glutamine  = hackett.glutamine  ./ hackett.glutamine(1);
kumar_lg.atp       = kumar_lg.atp ./ kumar_lg.atp(1);
kumar_hg.atp       = kumar_hg.atp ./ kumar_hg.atp(1);
hackett.atp        = hackett.atp  ./ hackett.atp(1);

% normalize protein to log2 fold change
WT_y_steady_10gL(:,[num_y.r: num_y.at]) = log2(WT_y_steady_10gL(:,[num_y.r: num_y.at]) ./ WT_y_steady_10gL(norm_indx,[num_y.r: num_y.at]));

hackett.r  = log2(hackett.r ./ hackett.r(1));
hackett.z  = log2(hackett.z ./ hackett.z(1));
hackett.gy = log2(hackett.gy ./ hackett.gy(1));
hackett.fe = log2(hackett.fe ./ hackett.fe(1));
hackett.gn = log2(hackett.gn ./ hackett.gn(1));
hackett.mt = log2(hackett.mt ./ hackett.mt(1));
hackett.sp = log2(hackett.sp ./ hackett.sp(1));
hackett.sd = log2(hackett.sd ./ hackett.sd(1));
hackett.as = log2(hackett.as ./ hackett.as(1));
hackett.at = log2(hackett.at ./ hackett.at(1));

%}

%% -------------------------------------------------------------------------------------------------------  
%  Organize simulation result in to matrix 
%  ------------------------------------------------------------------------------------------------------- 
% find corresponding simulation dilution rate for each experiment dilution rate; 
% find if simulation variable is measured in each experiment variable; if not measured in experiment set as 0; 

% get simulation dilution rate index of corresponding experiment data dilution rate  
D_temp = D_chem.*ones(numel(boer.dr),1); 
[~, sim_dr_indx.boer]     = min(abs(D_temp - boer.dr    ),[],2);
D_temp = D_chem.*ones(numel(kumar_lg.dr),1); 
[~, sim_dr_indx.kumar_lg] = min(abs(D_temp - kumar_lg.dr),[],2);
D_temp = D_chem.*ones(numel(kumar_hg.dr),1); 
[~, sim_dr_indx.kumar_hg] = min(abs(D_temp - kumar_hg.dr),[],2);

% onehot encoder for variable measured in each experiment [variable numbers by 1]; 
% measured = 1; not measured = 0 
sim_have_exp.boer     = ones(13,1);  % boer and hacekett 
sim_have_exp.kumar_lg = ones(13,1); 
sim_have_exp.kumar_hg = ones(13,1); 
sim_have_exp.boer(num_exp.cells)  = 0; % no exp measured
sim_have_exp.kumar_hg([num_exp.r:num_exp.at]) = 0; % no exp measured
sim_have_exp.kumar_lg([num_exp.cells,num_exp.r:num_exp.at]) = 0; % no exp measured

num_exp_per_variable = sim_have_exp.boer + sim_have_exp.kumar_lg + sim_have_exp.kumar_hg;

% pass sim_have_exp, sim_dr_indx in function chemo_sim_result_mat 
% to get simulation result vector
sim_result.boer     = chemo_sim_result_mat(sim_dr_indx.boer,     sim_have_exp.boer,     WT_y_steady_10gL, WT_met_reac_steady_10gL, num_y, num_flux);
sim_result.kumar_lg = chemo_sim_result_mat(sim_dr_indx.kumar_lg, sim_have_exp.kumar_lg, WT_y_steady_10gL, WT_met_reac_steady_10gL, num_y, num_flux);
sim_result.kumar_hg = chemo_sim_result_mat(sim_dr_indx.kumar_hg, sim_have_exp.kumar_hg, WT_y_steady_10gL, WT_met_reac_steady_10gL, num_y, num_flux);


%% -------------------------------------------------------------------------------------------------------  
%  Organize experiment result in to matrix 
%  ------------------------------------------------------------------------------------------------------- 

% initialize experiment result vector
exp_result.boer     =  zeros(size(sim_result.boer));
exp_result.kumar_lg =  zeros(size(sim_result.kumar_lg));
exp_result.kumar_hg =  zeros(size(sim_result.kumar_hg));

% boer result vector
exp_result.boer(num_exp.Jgl,:) = boer.Jgy';
exp_result.boer(num_exp.Jeh,:) = boer.Jfe_minus_Jgo';
exp_result.boer(num_exp.aa,:)  = hackett.glutamine';
exp_result.boer(num_exp.ae,:)  = hackett.atp';
exp_result.boer(num_exp.r,:)   = hackett.r';
exp_result.boer(num_exp.z,:)   = hackett.z';
exp_result.boer(num_exp.gy,:)  = hackett.gy';
exp_result.boer(num_exp.fe,:)  = hackett.fe';
exp_result.boer(num_exp.gn,:)  = hackett.gn';
exp_result.boer(num_exp.mt,:)  = hackett.mt';
exp_result.boer(num_exp.as,:)  = hackett.as';
exp_result.boer(num_exp.at,:)  = hackett.at';

% kumar lg result vector
exp_result.kumar_lg(num_exp.Jgl,:) = kumar_lg.Jgy';
exp_result.kumar_lg(num_exp.Jeh,:) = kumar_lg.Jfe_minus_Jgo';
exp_result.kumar_lg(num_exp.aa,:)  = kumar_lg.glutamine';
exp_result.kumar_lg(num_exp.ae,:)  = kumar_lg.atp';

% kumar hg result vector
exp_result.kumar_hg(num_exp.Jgl,:) = kumar_hg.Jgy';
exp_result.kumar_hg(num_exp.Jeh,:) = kumar_hg.Jfe_minus_Jgo';
exp_result.kumar_hg(num_exp.aa,:)  = kumar_hg.glutamine';
exp_result.kumar_hg(num_exp.ae,:)  = kumar_hg.atp';
exp_result.kumar_hg(num_exp.cells,:) = kumar_hg.cells';
exp_result.kumar_hg(num_exp.cells,end) = 0; % not compare the last one 

%% -------------------------------------------------------------------------------------------------------  
%  Data processing before calculating error 
%  ------------------------------------------------------------------------------------------------------- 

% lambda reg coefficient/weight for protein data 
% set as 1 --> no reg effect; or same weight for protein log2 data and metabolite data
lambda = 1; 

% add small number to prevent dividing 0; 
exp_result.boer     = exp_result.boer + eps(1); 
exp_result.kumar_lg = exp_result.kumar_lg + eps(1); 
exp_result.kumar_hg = exp_result.kumar_hg + eps(1); 

sim_result.boer     = sim_result.boer + eps(1); 
sim_result.kumar_lg = sim_result.kumar_lg + eps(1); 
sim_result.kumar_hg = sim_result.kumar_hg + eps(1); 

% Special treatment for Ethanol production data; Up to discussion
sim_result.boer(num_exp.Jeh,sim_result.boer(num_exp.Jeh,:)<0.3*10^4)         = eps(1); 
sim_result.kumar_lg(num_exp.Jeh,sim_result.kumar_lg(num_exp.Jeh,:)<0.3*10^4) = eps(1); 
sim_result.kumar_hg(num_exp.Jeh,sim_result.kumar_hg(num_exp.Jeh,:)<0.3*10^4) = eps(1); 

%% -------------------------------------------------------------------------------------------------------  
%  calculate errors 
%  ------------------------------------------------------------------------------------------------------- 

% calculate errors 
rel_error = sum(abs((abs(exp_result.boer./sim_result.boer) - 1).^2),2);
abs_error = sum(lambda * abs(exp_result.boer - sim_result.boer).^2,2); 
error.boer  = min(rel_error,abs_error);

rel_error = sum(abs((abs(exp_result.kumar_lg./sim_result.kumar_lg) - 1).^2),2);
abs_error = sum(lambda * abs(exp_result.kumar_lg - sim_result.kumar_lg).^2,2); 
error.kumar_lg  = min(rel_error,abs_error);

rel_error = sum(abs((abs(exp_result.kumar_hg./sim_result.kumar_hg) - 1).^2),2);
abs_error = sum(lambda * abs(exp_result.kumar_hg - sim_result.kumar_hg).^2,2); 
error.kumar_hg  = min(rel_error,abs_error);

error.chemo_per_variable = (error.boer+error.kumar_lg +error.kumar_hg)./num_exp_per_variable; 
error.chemo_sum = sum (error.chemo_per_variable);


end 


function sim_result = chemo_sim_result_mat(idx,sim_have_exp,WT_y_steady_10gL, WT_met_reac_steady_10gL, num_y, num_flux)
% variable output
sim_result = [ WT_met_reac_steady_10gL.flux(idx,num_flux.gy)';...
WT_met_reac_steady_10gL.flux(idx,num_flux.fe)' - WT_met_reac_steady_10gL.flux(idx,num_flux.gn)';
WT_y_steady_10gL(idx,end)';
WT_y_steady_10gL(idx,num_y.aa_in)';
WT_y_steady_10gL(idx,num_y.ae)';
WT_y_steady_10gL(idx,num_y.r)';
WT_y_steady_10gL(idx,num_y.z)';
WT_y_steady_10gL(idx,num_y.gy)';
WT_y_steady_10gL(idx,num_y.fe)';
WT_y_steady_10gL(idx,num_y.gn)';
WT_y_steady_10gL(idx,num_y.mt)';
WT_y_steady_10gL(idx,num_y.as)';
WT_y_steady_10gL(idx,num_y.at)';];

% set variable as 0 if not measured in exp.             
sim_result = sim_result.*sim_have_exp;            
end                
