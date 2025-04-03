function [t, y, met_reac, other_met_reac, sig, prot_syn, rib, g_rate, ribo_rate] = function_ode_solver(par, t_start, t_final, met, prot, R0, cells, mutant, snf1_vals, jgy_vals, mgl)

plt_labels = plot_labels;
%
tic
disp('Start batch simulation...')

% cell state
%{
cell_state = [table2array(struct2table(met))';  % 1 - 8
table2array(struct2table(prot))'; % 9 - 20 
R0;                               % 21
cells;                            % 22
];
%}

cell_state = [met';  % 1 - 8 %_% not yet modified
prot'; % 9 - 20              %_% not yet modified
R0;                               % 21
cells;                            % 22
];          

% set ode options
%Rel_tol  = 1.0E-09; 
%Abs_tol  = 1.0E-12; 
Rel_tol  = 1.0E-03; 
Abs_tol  = 1.0E-06; 
%max_step = 100;
tstart = tic;
options  = odeset('RelTol',Rel_tol, ...
'AbsTol',Abs_tol, ...
'NonNegative',(1:length(cell_state)),...
'Events',@(t,y) myevent(t,y,tstart));       

disp('Running batch simulation...')    
[t, y,te,ye,ie] = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals, jgy_vals, mgl), ...
[t_start, t_final], ...
cell_state, ...
options); 

y = real(y);                         

% calculate intermediates
met_reac.prot      = ones(numel(t), 6); 
met_reac.substrate = ones(numel(t), 6); 
met_reac.atp       = ones(numel(t), 6); 
met_reac.sig       = ones(numel(t), 6); 
met_reac.flux      = ones(numel(t), 6); 

other_met_reac.flux.po  = ones(numel(t), 1); 
other_met_reac.flux.eo  = ones(numel(t), 1);
other_met_reac.flux.ao  = ones(numel(t), 1);

prot_syn.alpha         = ones(numel(t), 8); 
prot_syn.tc            = ones(numel(t), 1);  
prot_syn.tu            = ones(numel(t), 1);  
prot_syn.t0            = ones(numel(t), 1);  
prot_syn.beta          = ones(numel(t), 8);
prot_syn.prot_syn_rate = ones(numel(t), 1); 
prot_syn.eIF_a_s       = ones(numel(t), 1); 
prot_syn.eIF_a_tau     = ones(numel(t), 1); 
prot_syn.prot_syn_rate = ones(numel(t), 1); 
prot_syn.rib.others    = ones(numel(t), length(plt_labels.others));

sig.snf1 = ones(numel(t), 1);
sig.tor  = ones(numel(t), 1); 

rib.rf    = ones(numel(t), 1);
rib.rat_j = ones(numel(t), 8);
rib.rat   = ones(numel(t), 1);
rib.ras_j = ones(numel(t), 8);
rib.ras   = ones(numel(t), 1);

g_rate    = ones(numel(t), 1);
ribo_rate = ones(numel(t), 1); 

% get intermediate values  
%
disp('Calculating intermediate values...')    
for k = 1:length(t)
%[dydt_t, sig_t, met_reac_t, prot_syn_rate_t, beta_t, alpha_t, ra_t, rat_t, r0_t, tRNA_t, eIF_a_t, eIF_a_s_t, eIF_a_tau_t, other_met_reac_t] = yeast_model(t_batch(k), y_batch(k,:)', par, mutant_type);
%[dydt_t, sig_t, met_reac_t, prot_syn_rate_t, beta_t, alpha_t, rib_t, tRNA_t, eIF_a_s_t, eIF_a_tau_t, other_met_reac_t, g_rate_t, ribo_rate_t] = yeast_model(t_batch(k), y_batch(k,:)', par, mutant_type);
[~, sig_t, met_reac_t, prot_syn_rate_t, beta_t, alpha_t, rib_t, tRNA_t, eIF_a_s_t, eIF_a_tau_t, other_met_reac_t, g_rate_t, ribo_rate_t] = yeast_model_update_ribosome(t(k), y(k,:)', par, mutant, snf1_vals, jgy_vals, mgl);

met_reac.prot(k,:)      = real(met_reac_t.prot)';
met_reac.substrate(k,:) = real(met_reac_t.substrate)';
met_reac.atp(k,:)       = real(met_reac_t.atp)';
met_reac.sig(k,:)       = real(met_reac_t.sig)';
met_reac.flux(k,:)      = real(met_reac_t.flux)';

other_met_reac.flux.po(k,:)  = real(other_met_reac_t.flux.po)';
other_met_reac.flux.eo(k,:)  = real(other_met_reac_t.flux.eo)';
%    other_met_reac.flux.ao(k,:)  = real(other_met_reac_t.flux.ao)';

prot_syn.alpha(k,:)       = real(table2array(struct2table(alpha_t))); 
prot_syn.tc(k)            = real(tRNA_t.tc); 
prot_syn.tu(k)            = real(tRNA_t.tu);
prot_syn.t0(k)            = real(tRNA_t.t0); 
%prot_syn.eIF_a(k)         = real(eIF_a_t)'; 
prot_syn.beta(k,:)        = real(beta_t)';
prot_syn.prot_syn_rate(k) = real(prot_syn_rate_t); 
prot_syn.eIF_a_s(k)       = real(eIF_a_s_t);
prot_syn.eIF_a_tau(k)     = real(eIF_a_tau_t);
prot_syn.prot_syn_rate(k) = real(prot_syn_rate_t);
prot_syn.rib.others(k,:)     = real(rib_t.others);    

sig.snf1(k) = real(sig_t.snf1);
sig.tor(k)  = real(sig_t.tor); 

rib.rf(k)      = real(rib_t.rf); 
%rib.raf(k)     = real(rib_t.raf);
%rib.ri(k)      = real(rib_t.ri);
rib.rat_j(k,:) = real(rib_t.rat_j)';
rib.rat(k)     = real(rib_t.rat);
rib.ras_j(k,:) = real(rib_t.ras_j)';
rib.ras(k)     = real(rib_t.ras); 

g_rate(k)      = real(g_rate_t); 
ribo_rate(k)   = real(ribo_rate_t); 
end 
%}                         
end 

function [values, isterminal, direction] = myevent(t, y,tstart)
%  Don't let t cross zero...a dummy "event" to illustrate how 
%  one might handle other events in conjunction with the time
%  constraint.  Not necessary, but I put it in just in case.
%  Don't let integration go for more than 5 seconds.
values = toc(tstart) < 5; % The value that we want to be zero
isterminal = 1;           % Halt integration 
direction  = 0;           % The zero can be approached from either direction
end 
