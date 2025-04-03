function [t_batch, y_batch, met_reac, other_met_reac, sig, prot_syn, rib, g_rate, ribo_rate] = function_pulse(cell_state, par, t_batch_start, t_batch_final, aa_ex, nh_ex,  gl, eh, cells, mutant, snf1_vals, jgy_vals, mgl)

%


plot_num;
plt_labels = plot_labels; 
Rel_tol  = 1.0E-03; 
Abs_tol  = 1.0E-06; 
tstart = tic;
options  = odeset('RelTol',Rel_tol, ...
'AbsTol',Abs_tol, ...
'NonNegative',(1:length(cell_state)),...
'Events',@(t,y) myevent(t,y,tstart));          


cell_state(num_y.gl)     = gl;
cell_state(num_y.eh)     = eh;
cell_state(num_y.aaex)   = aa_ex;
cell_state(num_y.nh4)    = nh_ex;
cell_state(num_y.cells)  = cells; 
% get initial cell state
% cell_state = get_cell_state_update_ribosome(cell_state_temp, par, mutant, snf1_vals, jgy_vals, mgl);
% cell_state(8) = 0; 
% cell_state(9) = 0; 

disp('Calculate initial conditions for batch simulation: done')  

disp('Running batch simulation...')    
[t_batch, y_batch,te,ye,ie] = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals, jgy_vals, mgl), ...
[t_batch_start t_batch_final], ...
cell_state, ...
options); 

y_batch = real(y_batch);                         

% calculate intermediates
met_reac.prot      = ones(numel(t_batch), length(fieldnames(num_flux))); 
met_reac.substrate = ones(numel(t_batch), length(fieldnames(num_flux))); 
met_reac.atp       = ones(numel(t_batch), length(fieldnames(num_flux))); 
met_reac.sig       = ones(numel(t_batch), length(fieldnames(num_flux))); 
met_reac.flux      = ones(numel(t_batch), length(fieldnames(num_flux))); 
met_reac.snf1      = ones(numel(t_batch), length(fieldnames(num_flux))); 
met_reac.pc        = ones(numel(t_batch), length(fieldnames(num_flux))); 

met_reac.k_val_glu       = ones(numel(t_batch), length(fieldnames(num_flux)));
met_reac.k_val_glu_fc    = ones(numel(t_batch), length(fieldnames(num_flux)));
met_reac.k_val_eth       = ones(numel(t_batch), length(fieldnames(num_flux)));
met_reac.k_val_eth_fc    = ones(numel(t_batch), length(fieldnames(num_flux)));

other_met_reac.flux.po  = ones(numel(t_batch), 1); 
other_met_reac.flux.eo  = ones(numel(t_batch), 1);
other_met_reac.flux.no  = ones(numel(t_batch), 1);

prot_syn.alpha         = ones(numel(t_batch), length(fieldnames(num_prot))); 
prot_syn.tc            = ones(numel(t_batch), 1);  
prot_syn.tu            = ones(numel(t_batch), 1);  
prot_syn.t0            = ones(numel(t_batch), 1);  
prot_syn.beta          = ones(numel(t_batch), length(fieldnames(num_prot)));
prot_syn.prot_syn_rate = ones(numel(t_batch), 1); 
prot_syn.eIF_a_s       = ones(numel(t_batch), 1); 
prot_syn.eIF_a_tau     = ones(numel(t_batch), 1); 
prot_syn.prot_syn_rate = ones(numel(t_batch), 1); 
prot_syn.rib.others    = ones(numel(t_batch), length(plt_labels.others));

sig.snf1 = ones(numel(t_batch), 1);
sig.tor  = ones(numel(t_batch), 1); 

rib.rf    = ones(numel(t_batch), 1);
rib.rat_j = ones(numel(t_batch), length(fieldnames(num_prot)));
rib.rat   = ones(numel(t_batch), 1);
rib.ras_j = ones(numel(t_batch), length(fieldnames(num_prot)));
rib.ras   = ones(numel(t_batch), 1);

g_rate    = ones(numel(t_batch), 1);
ribo_rate = ones(numel(t_batch), 1); 

% get intermediate values      
disp('Calculating intermediate values...')    
for k = 1:length(t_batch)
%[dydt_t, sig_t, met_reac_t, prot_syn_rate_t, beta_t, alpha_t, ra_t, rat_t, r0_t, tRNA_t, eIF_a_t, eIF_a_s_t, eIF_a_tau_t, other_met_reac_t] = yeast_model(t_batch(k), y_batch(k,:)', par, mutant_type);
%[dydt_t, sig_t, met_reac_t, prot_syn_rate_t, beta_t, alpha_t, rib_t, tRNA_t, eIF_a_s_t, eIF_a_tau_t, other_met_reac_t, g_rate_t, ribo_rate_t] = yeast_model(t_batch(k), y_batch(k,:)', par, mutant_type);
[~, sig_t, met_reac_t, prot_syn_rate_t, beta_t, alpha_t, rib_t, tRNA_t, eIF_a_s_t, eIF_a_tau_t, other_met_reac_t, g_rate_t, ribo_rate_t] = yeast_model_update_ribosome(t_batch(k), y_batch(k,:)', par, mutant, snf1_vals, jgy_vals, mgl);

met_reac.prot(k,:)      = real(met_reac_t.prot)';
met_reac.substrate(k,:) = real(met_reac_t.substrate)';
met_reac.atp(k,:)       = real(met_reac_t.atp)';
met_reac.sig(k,:)       = real(met_reac_t.sig)';
met_reac.flux(k,:)      = real(met_reac_t.flux)';
met_reac.snf1(k,:)      = real(met_reac_t.snf1)'; 
met_reac.pc(k,:)        = real(met_reac_t.pc)'; 
%met_reac.k_val_glu(k,:)       = real(met_reac_t.k_val_glu)';
%met_reac.k_val_glu_fc(k,:)    = real(met_reac_t.k_val_glu_fc)';
%met_reac.k_val_eth(k,:)       = real(met_reac_t.k_val_eth)';
%met_reac.k_val_eth_fc(k,:)    = real(met_reac_t.k_val_eth_fc)';

other_met_reac.flux.po(k,:)  = real(other_met_reac_t.flux.po)';
other_met_reac.flux.eo(k,:)  = real(other_met_reac_t.flux.eo)';
other_met_reac.flux.no(k,:)  = real(other_met_reac_t.flux.no)';

prot_syn.alpha(k,:)       = real(table2array(struct2table(alpha_t)));
prot_syn.rib.others(k,:)  = real(rib_t.others);    
prot_syn.tc(k)            = real(tRNA_t.tc); 
prot_syn.tu(k)            = real(tRNA_t.tu);
prot_syn.t0(k)            = real(tRNA_t.t0); 
%prot_syn.eIF_a(k)         = real(eIF_a_t)'; 
prot_syn.beta(k,:)        = real(beta_t)';
prot_syn.prot_syn_rate(k) = real(prot_syn_rate_t); 
prot_syn.eIF_a_s(k)       = real(eIF_a_s_t);
prot_syn.eIF_a_tau(k)     = real(eIF_a_tau_t);
prot_syn.prot_syn_rate(k) = real(prot_syn_rate_t);

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


end 


function [values, isterminal, direction] = myevent(t, y,tstart)
%  Don't let t cross zero...a dummy "event" to illustrate how 
%  one might handle other events in conjunction with the time
%  constraint.  Not necessary, but I put it in just in case.
%  Don't let integration go for more than 5 seconds.
values = toc(tstart) < 15; % The value that we want to be zero
isterminal = 1;           % Halt integration 
direction  = 0;           % The zero can be approached from either direction
end 
