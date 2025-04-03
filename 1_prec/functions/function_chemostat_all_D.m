function [WT_y_steady, WT_met_reac_steady, WT_other_met_reac_steady, WT_prot_syn_steady,...
WT_sig_steady, WT_rib_steady, WT_g_rate_steady, WT_ribo_rate_steady]...
= function_chemostat_all_D(par, D_chem,  t_chem_start, t_chem_final, aa_ex,...
nh_ex, gl, eh, cells, mutant, snf1_vals, jgy_vals, mgl, ...
plt_labels, num_y, num_flux, num_prot)







WT_y_steady = ones(numel(D_chem), numel(fieldnames(num_y)));

WT_met_reac_steady.prot      = ones(numel(D_chem), length(fieldnames(num_flux))); 
WT_met_reac_steady.substrate = ones(numel(D_chem), length(fieldnames(num_flux))); 
WT_met_reac_steady.atp       = ones(numel(D_chem), length(fieldnames(num_flux))); 
WT_met_reac_steady.sig       = ones(numel(D_chem), length(fieldnames(num_flux))); 
WT_met_reac_steady.flux      = ones(numel(D_chem), length(fieldnames(num_flux))); 
WT_met_reac_steady.snf1      = ones(numel(D_chem), length(fieldnames(num_flux))); 
WT_met_reac_steady.pc        = ones(numel(D_chem), length(fieldnames(num_flux))); 

WT_other_met_reac_steady.flux.po  = ones(numel(D_chem), 1); 
WT_other_met_reac_steady.flux.eo  = ones(numel(D_chem), 1);
WT_other_met_reac_steady.flux.no  = ones(numel(D_chem), 1);

WT_prot_syn_steady.alpha         = ones(numel(D_chem), length(fieldnames(num_prot))); 
WT_prot_syn_steady.tc            = ones(numel(D_chem), 1);  
WT_prot_syn_steady.tu            = ones(numel(D_chem), 1);  
WT_prot_syn_steady.t0            = ones(numel(D_chem), 1);  
WT_prot_syn_steady.beta          = ones(numel(D_chem), length(fieldnames(num_prot)));
WT_prot_syn_steady.prot_syn_rate = ones(numel(D_chem), 1); 
WT_prot_syn_steady.eIF_a_s       = ones(numel(D_chem), 1); 
WT_prot_syn_steady.eIF_a_tau     = ones(numel(D_chem), 1); 
WT_prot_syn_steady.prot_syn_rate = ones(numel(D_chem), 1); 
WT_prot_syn_steady.rib.others    = ones(numel(D_chem), length(plt_labels.others));

WT_sig_steady.snf1 = ones(numel(D_chem), 1);
WT_sig_steady.tor  = ones(numel(D_chem), 1); 

WT_rib_steady.rf    = ones(numel(D_chem), 1);
WT_rib_steady.rat_j = ones(numel(D_chem), length(fieldnames(num_prot)));
WT_rib_steady.rat   = ones(numel(D_chem), 1);
WT_rib_steady.ras_j = ones(numel(D_chem), length(fieldnames(num_prot)));
WT_rib_steady.ras   = ones(numel(D_chem), 1);

WT_g_rate_steady    = ones(numel(D_chem), 1);
WT_ribo_rate_steady = ones(numel(D_chem), 1); 


new_cell_state.overwrite = 'False';

for k = 1:numel(D_chem)
% run simulation
% fprintf('D = %g \n', D_chem(k))
par.D = D_chem(k); 
[chem_min_y, chem_min_met_reac, chem_min_other_met_reac, chem_min_sig, chem_min_prot_syn, chem_min_rib, chem_min_g_rate, chem_min_ribo_rate] = function_chemostat(par, t_chem_start, t_chem_final, aa_ex, nh_ex, gl, eh, cells, mutant, snf1_vals, jgy_vals, mgl, new_cell_state);

WT_y_steady(k,:) = chem_min_y; 
new_cell_state.overwrite = 'True'; % Use the steady-state solution as the initial condition
new_cell_state.values = chem_min_y; 

WT_met_reac_steady.prot(k,:)      = chem_min_met_reac.prot; 
WT_met_reac_steady.substrate(k,:) = chem_min_met_reac.substrate; 
WT_met_reac_steady.atp(k,:)       = chem_min_met_reac.atp; 
WT_met_reac_steady.sig(k,:)       = chem_min_met_reac.sig; 
WT_met_reac_steady.flux(k,:)      = chem_min_met_reac.flux; 
WT_met_reac_steady.snf1(k,:)      = chem_min_met_reac.snf1; 
WT_met_reac_steady.pc(k,:)        = chem_min_met_reac.pc;

WT_other_met_reac_steady.flux.po(k,:)  = chem_min_other_met_reac.flux.po; 
WT_other_met_reac_steady.flux.eo(k,:)  = chem_min_other_met_reac.flux.eo;
WT_other_met_reac_steady.flux.no(k,:)  = chem_min_other_met_reac.flux.no;

WT_prot_syn_steady.alpha(k,:)         = chem_min_prot_syn.alpha; 
WT_prot_syn_steady.tc(k,:)            = chem_min_prot_syn.tc;  
WT_prot_syn_steady.tu(k,:)            = chem_min_prot_syn.tu;  
WT_prot_syn_steady.t0(k,:)            = chem_min_prot_syn.t0;  
WT_prot_syn_steady.beta(k,:)          = chem_min_prot_syn.beta;
WT_prot_syn_steady.prot_syn_rate(k,:) = chem_min_prot_syn.prot_syn_rate; 
WT_prot_syn_steady.eIF_a_s(k,:)       = chem_min_prot_syn.eIF_a_s; 
WT_prot_syn_steady.eIF_a_tau(k,:)     = chem_min_prot_syn.eIF_a_tau; 
WT_prot_syn_steady.prot_syn_rate(k,:) = chem_min_prot_syn.prot_syn_rate; 
WT_prot_syn_steady.rib.others(k,:)    = chem_min_rib.others; 

WT_sig_steady.snf1(k,:) = chem_min_sig.snf1;
WT_sig_steady.tor(k,:)  = chem_min_sig.tor; 

WT_rib_steady.rf(k,:)    = chem_min_rib.rf;
WT_rib_steady.rat_j(k,:) = chem_min_rib.rat_j;
WT_rib_steady.rat(k,:)   = chem_min_rib.rat;
WT_rib_steady.ras_j(k,:) = chem_min_rib.ras_j;
WT_rib_steady.ras(k,:)   = chem_min_rib.ras;

WT_g_rate_steady(k,:) = chem_min_g_rate;
WT_ribo_rate_steady(k,:) = chem_min_ribo_rate; 

end 

end 
