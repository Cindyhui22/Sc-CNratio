function error = function_batch_sim_exp_error(bart_WT_batch_min_t,   zamp_WT_batch_min_t,   murphy_WT_batch_YPD_t, ...
bart_WT_batch_min_y,   zamp_WT_batch_min_y,   murphy_WT_batch_YPD_y, ...
bart_WT_batch_min_sig, zamp_WT_batch_min_sig, murphy_WT_batch_YPD_sig,...
par, num_y, num_flux, num_prot)
%
% This function calculate the error between experiment and simulation for chemostat
% Relative Error = abs(exp data/sim data - 1)
% Absolute Error = abs(exp data - sim data)
%
% Example: 
% 
%
% 
%
%% -------------------------------------------------------------------------------------------------------  
% experiment, simulation matrix index 
%
num_exp.gl    = 1;
num_exp.eh    = 2; 
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

num_exp.bart  = 1; 
num_exp.zamp  = 2; 
num_exp.murphy= 3; 

%% -------------------------------------------------------------------------------------------------------  
%  load exp data
%  ------------------------------------------------------------------------------------------------------- 
%
data_bartolomeo; 
data_murphy; 
data_zampar; 
data_batch;
data_4_prec; 

%% -------------------------------------------------------------------------------------------------------  
%  normalize precursor, amino acid, atp, protein with respect to first data
%  point
%  ------------------------------------------------------------------------------------------------------- 
% precursor, amino acid, atp
bart_WT_batch_min_y(:,num_y.pc)   = bart_WT_batch_min_y(:,num_y.pc)   ./ bart_WT_batch_min_y(1,num_y.pc);
zamp_WT_batch_min_y(:,num_y.pc)   = zamp_WT_batch_min_y(:,num_y.pc)   ./ zamp_WT_batch_min_y(1,num_y.pc);
murphy_WT_batch_YPD_y(:,num_y.pc) = murphy_WT_batch_YPD_y(:,num_y.pc) ./ murphy_WT_batch_YPD_y(1,num_y.pc);

bart_WT_batch_min_y(:,num_y.ae)   = bart_WT_batch_min_y(:,num_y.ae)   ./ bart_WT_batch_min_y(1,num_y.ae);
zamp_WT_batch_min_y(:,num_y.ae)   = zamp_WT_batch_min_y(:,num_y.ae)   ./ zamp_WT_batch_min_y(1,num_y.ae);
murphy_WT_batch_YPD_y(:,num_y.ae) = murphy_WT_batch_YPD_y(:,num_y.ae) ./ murphy_WT_batch_YPD_y(1,num_y.ae);

bart_WT_batch_min_y(:,num_y.aa_in)   = bart_WT_batch_min_y(:,num_y.aa_in)   ./ bart_WT_batch_min_y(1,num_y.aa_in);
zamp_WT_batch_min_y(:,num_y.aa_in)   = zamp_WT_batch_min_y(:,num_y.aa_in)   ./ zamp_WT_batch_min_y(1,num_y.aa_in);
murphy_WT_batch_YPD_y(:,num_y.aa_in) = murphy_WT_batch_YPD_y(:,num_y.aa_in) ./ murphy_WT_batch_YPD_y(1,num_y.aa_in);

zampar_4prec.g6p               = zampar_4prec.g6p              ./ zampar_4prec.g6p(1);
zampar_4prec.f6p_g1p           = zampar_4prec.f6p_g1p          ./ zampar_4prec.f6p_g1p(1);
zampar_4prec.f16bp             = zampar_4prec.f16bp            ./ zampar_4prec.f16bp(1);
solis_escalante_4prec.g6p      = solis_escalante_4prec.g6p     ./ solis_escalante_4prec.g6p(1);
solis_escalante_4prec.f6p      = solis_escalante_4prec.f6p     ./ solis_escalante_4prec.f6p(1);
solis_escalante_4prec.f16bp    = solis_escalante_4prec.f16bp   ./ solis_escalante_4prec.f16bp(1);
solis_escalante_4prec.g3p      = solis_escalante_4prec.g3p     ./ solis_escalante_4prec.g3p(1);
zampar_4prec.dhap              = zampar_4prec.dhap             ./ zampar_4prec.dhap(1);
zampar_4prec.up2pg             = zampar_4prec.up2pg            ./ zampar_4prec.up2pg(1);
zampar_4prec.pep               = zampar_4prec.pep              ./ zampar_4prec.pep(1);
solis_escalante_4prec.dhap     = solis_escalante_4prec.dhap    ./ solis_escalante_4prec.dhap(1);
solis_escalante_4prec.up2pg    = solis_escalante_4prec.up2pg   ./ solis_escalante_4prec.up2pg(1);
solis_escalante_4prec.up3pg    = solis_escalante_4prec.up3pg   ./ solis_escalante_4prec.up3pg(1);
solis_escalante_4prec.pep      = solis_escalante_4prec.pep     ./ solis_escalante_4prec.pep(1);
zampar_4prec.pyruvate          = zampar_4prec.pyruvate         ./ zampar_4prec.pyruvate(1);
solis_escalante_4prec.pyruvate = solis_escalante_4prec.pyruvate./ solis_escalante_4prec.pyruvate(1);

zampar.glutamine(:,2) = zampar.glutamine(:,2) ./ zampar.glutamine(1,2);
bisschops.atp(1:4,1)  = bisschops.atp(1:4,1) + 13;                    % only use data point 1-4
bisschops.atp(1:4,2) = bisschops.atp(1:4,2) ./ bisschops.atp(1,2);
zampar.atp(:,2)    = zampar.atp(:,2)    ./ zampar.atp(1,2);

murphy.glucose(5:end,2) =0; % do not consider none-zeros after gl. dep
% protein 
bart_WT_batch_min_y(:,[num_y.r: num_y.at])   = log2(bart_WT_batch_min_y(:,[num_y.r: num_y.at])   ./ bart_WT_batch_min_y(1,[num_y.r: num_y.at]));
zamp_WT_batch_min_y(:,[num_y.r: num_y.at])   = log2(zamp_WT_batch_min_y(:,[num_y.r: num_y.at])   ./ zamp_WT_batch_min_y(1,[num_y.r: num_y.at]));
murphy_WT_batch_YPD_y(:,[num_y.r: num_y.at]) = log2(murphy_WT_batch_YPD_y(:,[num_y.r: num_y.at]) ./ murphy_WT_batch_YPD_y(1,[num_y.r: num_y.at]));

bartolomeo.r  = log2(bartolomeo.r ./ bartolomeo.r(1));
bartolomeo.z  = log2(bartolomeo.z ./ bartolomeo.z(1));
bartolomeo.gy = log2(bartolomeo.gy ./ bartolomeo.gy(1));
bartolomeo.fe = log2(bartolomeo.fe ./ bartolomeo.fe(1));
bartolomeo.gn = log2(bartolomeo.gn ./ bartolomeo.gn(1));
bartolomeo.mt = log2(bartolomeo.mt ./ bartolomeo.mt(1));
bartolomeo.sp = log2(bartolomeo.sp ./ bartolomeo.sp(1));
bartolomeo.sd = log2(bartolomeo.sd ./ bartolomeo.sd(1));
bartolomeo.as = log2(bartolomeo.as ./ bartolomeo.as(1));
bartolomeo.at = log2(bartolomeo.at ./ bartolomeo.at(1));

murphy.r  = log2(murphy.r ./ murphy.r(1));
murphy.z  = log2(murphy.z ./ murphy.z(1));
murphy.gy = log2(murphy.gy ./ murphy.gy(1));
murphy.fe = log2(murphy.fe ./ murphy.fe(1));
murphy.gn = log2(murphy.gn ./ murphy.gn(1));
murphy.mt = log2(murphy.mt ./ murphy.mt(1));
murphy.sp = log2(murphy.sp ./ murphy.sp(1));
murphy.sd = log2(murphy.sd ./ murphy.sd(1));
murphy.as = log2(murphy.as ./ murphy.as(1));
murphy.at = log2(murphy.at ./ murphy.at(1));


%% -------------------------------------------------------------------------------------------------------  
%  Organize simulation protein result in to matrix 
%  ------------------------------------------------------------------------------------------------------- 
% find corresponding simulation time for each experiment time; 

% initialize simulation protein result vector
sim_prot_result.bart   = zeros(10,numel(bartolomeo.time_prot)); % 10 proteins by bart time point  
sim_prot_result.murphy = zeros(10,numel(murphy.time_prot));     % 10 proteins by murhpy time point


% get simulation time index of corresponding experiment time 
t_temp                      = bart_WT_batch_min_t'.*ones(numel(bartolomeo.time_prot),1); 
[~, sim_prot_indx.bart]     = min(abs(t_temp - bartolomeo.time_prot),[],2);
t_temp                      = murphy_WT_batch_YPD_t'.*ones(numel(murphy.time_prot),1); 
[~, sim_prot_indx.murphy]   = min(abs(t_temp - murphy.time_prot),[],2);

% onehot encoder for variable measured in each experiment [variable numbers by 1]; 
% measured = 1; not measured = 0 
sim_have_exp.bart     = ones(13-num_exp.ae,1);  
sim_have_exp.murphy   = ones(13-num_exp.ae,1); 

% pass sim_prot_indx, sim_have_exp in function batch_sim_result_mat 
% to get simulation protein result vector
sim_prot_result.bart     = batch_sim_result_mat(sim_prot_indx.bart,    sim_have_exp.bart,     bart_WT_batch_min_y, num_y, num_flux);
sim_prot_result.murphy   = batch_sim_result_mat(sim_prot_indx.murphy,  sim_have_exp.murphy,    murphy_WT_batch_YPD_y, num_y, num_flux);


%% -------------------------------------------------------------------------------------------------------  
%  Organize experiment protein result in to matrix 
%  ------------------------------------------------------------------------------------------------------- 

% initializ eexperiment result vector
exp_prot_result.bart     =  zeros(size(sim_prot_result.bart));
exp_prot_result.murphy   =  zeros(size(sim_prot_result.murphy));

% put exp result in matrix 
% bart result vector
exp_prot_result.bart(num_exp.r-num_exp.ae,:)   = bartolomeo.r';
exp_prot_result.bart(num_exp.z-num_exp.ae,:)   = bartolomeo.z';
exp_prot_result.bart(num_exp.gy-num_exp.ae,:)  = bartolomeo.gy';
exp_prot_result.bart(num_exp.fe-num_exp.ae,:)  = bartolomeo.fe';
exp_prot_result.bart(num_exp.gn-num_exp.ae,:)  = bartolomeo.gn';
exp_prot_result.bart(num_exp.mt-num_exp.ae,:)  = bartolomeo.mt';
exp_prot_result.bart(num_exp.as-num_exp.ae,:)  = bartolomeo.as';
exp_prot_result.bart(num_exp.at-num_exp.ae,:)  = bartolomeo.at';

% murphy result vector
exp_prot_result.murphy(num_exp.r-num_exp.ae,:)   = murphy.r';
exp_prot_result.murphy(num_exp.z-num_exp.ae,:)   = murphy.z';
exp_prot_result.murphy(num_exp.gy-num_exp.ae,:)  = murphy.gy';
exp_prot_result.murphy(num_exp.fe-num_exp.ae,:)  = murphy.fe';
exp_prot_result.murphy(num_exp.gn-num_exp.ae,:)  = murphy.gn';
exp_prot_result.murphy(num_exp.mt-num_exp.ae,:)  = murphy.mt';
exp_prot_result.murphy(num_exp.as-num_exp.ae,:)  = murphy.as';
exp_prot_result.murphy(num_exp.at-num_exp.ae,:)  = murphy.at';


%% -------------------------------------------------------------------------------------------------------  
%  Data processing before calculating protein error 
%  ------------------------------------------------------------------------------------------------------- 
% lambda reg coefficient/weight for protein data 
% set as 1 --> no reg effect; or same weight for protein log2 data and metabolite data
lambda = 1; 

% add small number to prevent dividing 0; 
exp_prot_result.bart     = exp_prot_result.bart   + eps(1); 
exp_prot_result.murphy   = exp_prot_result.murphy + eps(1); 
sim_prot_result.bart     = sim_prot_result.bart   + eps(1); 
sim_prot_result.murphy   = sim_prot_result.murphy + eps(1);

%% -------------------------------------------------------------------------------------------------------  
%  Compare sim and exp protein matrix and calculate error 
%  ------------------------------------------------------------------------------------------------------- 
% initialize 
% error between sim and exp computes from exp_prot_result, sim_prot_result
error.prot.all          = zeros(8,3); % 10 proteins by 3 different exp  

% protein error
rel_prot_error.bart     = abs((abs(exp_prot_result.bart./sim_prot_result.bart) - 1));
rel_prot_error.murphy   = abs((abs(exp_prot_result.murphy./sim_prot_result.murphy) - 1));
abs_prot_error.bart     = lambda * abs(exp_prot_result.bart - sim_prot_result.bart); 
abs_prot_error.murphy   = lambda * abs(exp_prot_result.murphy - sim_prot_result.murphy); 

error.prot.bart         = sum(min(rel_prot_error.bart ,abs_prot_error.bart ),2);
error.prot.murphy       = sum(min(rel_prot_error.murphy ,abs_prot_error.murphy ),2);
error.prot.zamp         = zeros(size(error.prot.murphy));  % did not measure protein; 


error.prot.all          = [error.prot.bart, error.prot.zamp, error.prot.murphy]; 


%% -------------------------------------------------------------------------------------------------------  
%  Organize simulation and experiment metabolite result in to matrix 
%  ------------------------------------------------------------------------------------------------------- 

% initialize error between sim and exp
error.met.all    = zeros(5,3); % 6 met by 3 different exp  

% pass in simulation ys,ts and experiment ye,te to batch_sim_met_at_exp_t
% to find simulation ys' at corresponding te, compute error 
[sim_met_result.bart_gl,     exp_met_result.bart_gl,     rel_prot_error.bart_gl,  abs_prot_error.bart_gl,     error.met.all(num_exp.gl, num_exp.bart)]      = batch_sim_met_at_exp_t(bart_WT_batch_min_t, bart_WT_batch_min_y, num_y.gl, bartolomeo.glucose, lambda);
[sim_met_result.bart_eh,     exp_met_result.bart_eh,     rel_prot_error.bart_eh,  abs_prot_error.bart_eh,     error.met.all(num_exp.eh, num_exp.bart)]      = batch_sim_met_at_exp_t(bart_WT_batch_min_t, bart_WT_batch_min_y, num_y.eh, bartolomeo.ethanol, lambda);
[sim_met_result.bart_cells,  exp_met_result.bart_cells,  rel_prot_error.bart_cells,abs_prot_error.bart_cells, error.met.all(num_exp.cells, num_exp.bart)]   = batch_sim_met_at_exp_t(bart_WT_batch_min_t, bart_WT_batch_min_y, num_y.cells, bartolomeo.cells,lambda);
[sim_met_result.bart_ae,     exp_met_result.bart_ae,     rel_prot_error.bart_ae,  abs_prot_error.bart_ae,     error.met.all(num_exp.ae, num_exp.bart)]      = batch_sim_met_at_exp_t(bart_WT_batch_min_t, bart_WT_batch_min_y, num_y.ae, bisschops.atp(1:4,:),lambda);

[sim_met_result.zamp_gl,     exp_met_result.zamp_gl,     rel_prot_error.zamp_gl,  abs_prot_error.zamp_gl,     error.met.all(num_exp.gl, num_exp.zamp)]      = batch_sim_met_at_exp_t(zamp_WT_batch_min_t, zamp_WT_batch_min_y, num_y.gl, zampar.glucose, lambda);
[sim_met_result.zamp_eh,     exp_met_result.zamp_eh,     rel_prot_error.zamp_eh,  abs_prot_error.zamp_eh,     error.met.all(num_exp.eh, num_exp.zamp)]      = batch_sim_met_at_exp_t(zamp_WT_batch_min_t, zamp_WT_batch_min_y, num_y.eh, zampar.ethanol, lambda);
[sim_met_result.zamp_cells,  exp_met_result.zamp_cells,  rel_prot_error.zamp_cells, abs_prot_error.zamp_cells,error.met.all(num_exp.cells, num_exp.zamp)]   = batch_sim_met_at_exp_t(zamp_WT_batch_min_t, zamp_WT_batch_min_y, num_y.cells, zampar.cells,lambda);
[sim_met_result.zamp_ae,     exp_met_result.zamp_ae,     rel_prot_error.zamp_ae,  abs_prot_error.zamp_ae,     error.met.all(num_exp.ae, num_exp.zamp)]      = batch_sim_met_at_exp_t(zamp_WT_batch_min_t, zamp_WT_batch_min_y, num_y.ae, zampar.atp, lambda);
[sim_met_result.zamp_aa,     exp_met_result.zamp_aa,     rel_prot_error.zamp_aa,  abs_prot_error.zamp_aa,     error.met.all(num_exp.aa, num_exp.zamp)]      = batch_sim_met_at_exp_t(zamp_WT_batch_min_t, zamp_WT_batch_min_y, num_y.aa_in, zampar.glutamine,lambda);

[sim_met_result.murphy_gl,   exp_met_result.murphy_gl,   rel_prot_error.murphy_gl, abs_prot_error.murphy_gl,  error.met.all(num_exp.gl,    num_exp.murphy)] = batch_sim_met_at_exp_t(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y, num_y.gl,    murphy.glucose,lambda);
[sim_met_result.murphy_cells,exp_met_result.murphy_cells,rel_prot_error.murphy_cells, abs_prot_error.bart_gl, error.met.all(num_exp.cells, num_exp.murphy)] = batch_sim_met_at_exp_t(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y, num_y.cells, murphy.cells, lambda);


%% -------------------------------------------------------------------------------------------------------  
%  Combine protein and metabolite error matrix 
%  ------------------------------------------------------------------------------------------------------- 

sim_have_exp.bart     = ones(13,1);  
sim_have_exp.zamp     = ones(13,1); 
sim_have_exp.murphy   = ones(13,1); 
sim_have_exp.bart([num_exp.aa])                                      = 0;  % no exp measured
sim_have_exp.zamp([num_exp.r:num_exp.at])                            = 0;  % no exp measured
sim_have_exp.murphy([num_exp.eh,num_exp.aa,num_exp.ae])  = 0;  % no exp measured

num_exp_per_variable = sim_have_exp.bart + sim_have_exp.zamp + sim_have_exp.murphy;


error.batch_per_variable = sum(([error.met.all; error.prot.all].^2)./num_exp_per_variable,2);  
error.batch_sum = sum(error.batch_per_variable);




end 


function sim_result = batch_sim_result_mat(idx,sim_have_exp, y, num_y, num_flux)
% variable output
sim_result = [  y(idx,num_y.r)';
y(idx,num_y.z)';
y(idx,num_y.gy)';
y(idx,num_y.fe)';
y(idx,num_y.gn)';
y(idx,num_y.mt)';
y(idx,num_y.as)';
y(idx,num_y.at)';];

% set variable as 0 if not measured in exp.             
sim_result = sim_result.*sim_have_exp;            
end                



function [sim_met_result,exp_met_result, rel_prot_error,  abs_prot_error, error] = batch_sim_met_at_exp_t(sim_t, sim_y, num, exp, lambda)

t_temp                   = sim_t'.*ones(numel(exp(:,1)),1); 
[~, sim_met_indx]        = min(abs(t_temp - exp(:,1)),[],2);
sim_met_result.t         = sim_t(sim_met_indx);
sim_met_result.y         = sim_y(sim_met_indx,num) + eps(1);
exp_met_result.y         = exp(:,2) + eps(1);
exp_met_result.t         = exp(:,1);

rel_prot_error          = sum(abs((abs(exp_met_result.y./sim_met_result.y) - 1))./numel(exp_met_result.t));
abs_prot_error          = sum(lambda * abs(exp_met_result.y - sim_met_result.y)./numel(exp_met_result.t)); 
error                   = min(rel_prot_error, abs_prot_error);


end
