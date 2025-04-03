function [beta, alpha, rib, tRNA, eIF_a_s, eIF_a_tau] = module_gene_expression_update_ribosome(r0, met, prot, snf1, tor, par, mutant)


% useful functions
hill_transc_s_tau = @(s, t, ep, xi_s, xi_t, w_s, w_t, theta_s, theta_t)(ep + (xi_s .* ((s/w_s).^theta_s)) ...
+ (xi_t .* ((t/w_t).^theta_t))) ...
./ (1 + ((s/w_s).^theta_s) + ((t/w_t).^theta_t)); 

%{                                                                      
hill_transc_s_tau_aain = @(s, t, aa, ep, xi_s, xi_t, xi_aa, w_s, w_t, w_aa, theta_s, theta_t, theta_aa) ...
(ep + (xi_s .* ((s/w_s).^theta_s)) + (xi_t .* ((t/w_t).^theta_t)) + (xi_aa .* ((t/w_aa).^theta_aa))) ...
./ (1 + ((s/w_s).^theta_s) + ((t/w_t) .^theta_t) + ((s/w_aa).^theta_aa)); 
%}                                                                                                                                          

hill_transc_s = @(s, ep, xi_s, w_s, theta_s)(ep + (xi_s .* ((s/w_s).^theta_s)) ...
./ (1 + ((s/w_s).^theta_s)));   

hill_transc_tau = @(tau, ep, xi_t, w_t, theta_t)(ep + (xi_t .* ((tau/w_t).^theta_t)) ...
./ (1 + ((tau/w_t).^theta_t))); %_% 
%_% hill_transc_tau = @(tau, ep, xi_t, w_t, theta_t)(ep + (xi_t .* ((s/w_t).^theta_s)) ...
%_% ./ (1 + ((s/w_t).^theta_s))); 

%hill_transc_met  = @(s, ep, xi, w, theta)(ep + xi.*(s/w).^theta)./(1 + (s/w).^theta);

%hill_transc_tau_aaex = @(t, aa, ep, xi_t, xi_aa, w_t, w_aa, theta_t, theta_aa)(ep + xi_t.*(t/w_t).^theta_t + xi_aa.*(aa/w_aa).^theta_aa)./(1 + (t/w_t).^theta_t + (aa/w_aa).^theta_aa);

%% metabolites 

aa_ex = met(1); 
%gl    = met(2); 
aa_in = met(5); 
ae    = met(6);

%% transcription factors 

%
if strcmp(mutant, 'hap4_underexpress') || strcmp(mutant, 'const_gl_eh_aaex_cell_hap_under')
% wild type    
% hap 4 knock out
par.xi_gn_s = 0.35;    
par.xi_mt_s = 0.35; 
%par.xi_mt_s =0.000001;  
end
%}
%
if strcmp(mutant, 'hap4_overexpress') || strcmp(mutant, 'const_gl_eh_aaex_cell_hap_over')
% hap 4 OE
par.xi_gn_s = 2; 
par.xi_mt_s = 2; 
end
%}


alpha.r  = par.a0_r + (1 - par.a0_r)*hill_transc_s_tau(snf1, tor, par.ep_r, ...
par.xi_r_s, par.xi_r_tau, ...
par.w_r_s, par.w_r_tau, ...
par.theta_r_s, par.theta_r_tau);
alpha.z  = 1;
alpha.gy = 1;
alpha.fe = 1;
alpha.gn = par.a0_gn + (1 - par.a0_gn)*hill_transc_s(snf1, par.ep_gn, par.xi_gn_s, par.w_gn_s, par.theta_gn_s);
alpha.mt = par.a0_mt + (1 - par.a0_mt)*hill_transc_s(snf1, par.ep_mt, par.xi_mt_s, par.w_mt_s, par.theta_mt_s);
%alpha.as = hill_transc_s_tau_aain(snf1, tor, aa_in, ...
%par.ep_as, ...
%par.xi_as_s, par.xi_as_tau, par.xi_as_aa_in, ...
%par.w_as_s, par.w_as_tau, par.w_as_aa_in, ...
%par.theta_as_s, par.theta_as_tau, par.theta_as_aa_in);


%
alpha.as = par.a0_as + (1 - par.a0_as)*hill_transc_s_tau(snf1, tor, par.ep_as, ...
par.xi_as_s, par.xi_as_tau, ...
par.w_as_s, par.w_as_tau, ...
par.theta_as_s, par.theta_as_tau);  
%}
%{                         
alpha.as = hill_transc_s_tau_aain(snf1, tor, aa_in, ...
par.ep_as, ...
par.xi_as_s, par.xi_as_tau, par.xi_as_aa_in, ...
par.w_as_s, par.w_as_tau, par.w_as_aa_in, ...
par.theta_as_s, par.theta_as_tau, par.theta_as_aa_in); 
%}                              

% alpha.at = par.a0_at + (1 - par.a0_at)*hill_transc_tau_aaex(tor, aa_ex, par.ep_at, ...
%                                                             par.xi_at_tau, par.xi_at_aaex, ...
%                                                             par.w_at_tau, par.w_at_aaex, ...
%                                                             par.theta_at_tau, par.theta_at_aaex);

alpha.at = 1; 

 

%_% alpha.lp = par.a0_lp + (1 - par.a0_lp)*hill_transc_s_tau(snf1, tor, par.ep_lp, ...
%_% par.xi_lp_s, par.xi_lp_tau, ...
%_% par.w_lp_s, par.w_lp_tau, ...
%_% par.theta_lp_s, par.theta_lp_tau);
%_% 
alpha.lp = par.a0_lp + (1 - par.a0_lp)*hill_transc_tau(tor, par.ep_lp, par.xi_lp_tau, par.w_lp_tau, par.theta_lp_tau);

alpha.lo = par.a0_lo + (1 - par.a0_lo)*hill_transc_s_tau(snf1, tor, par.ep_lo, ...
par.xi_lo_s, par.xi_lo_tau, ...
par.w_lo_s, par.w_lo_tau, ...
par.theta_lo_s, par.theta_lo_tau);

% 
% alpha.lp = 1; 
% alpha.lo = 1; 

% alpha.at = par.alpha_at_0*(0.5^8/(tor^8+0.5^8)) + (1- par.alpha_at_0)*(aa_ex^2/(aa_ex^2+5000^2));

%% rRNA 
%rRNA = par.mr * p_r; 

%% tRNA 
%
tRNA.t0 = par.mt  * r0;                                                                         % total
tRNA.tc = tRNA.t0 * hill(aa_in, par.K_tc_aa, par.n_tc_aa) * hill(ae, par.K_tc_ae, par.n_tc_ae); % charged
tRNA.tu = tRNA.t0 - tRNA.tc;                                                                    % uncharged

%% ribosome partitioning

% term for tRNA binding to ribosomes 
trna_tc = tRNA.tc./par.K_rib_c;  
trna_tu = tRNA.tu./par.K_rib_u;
trna_contri_tc = trna_tc./(1 + trna_tc + trna_tu); 
trna_contri_tu = trna_tu./(1 + trna_tc + trna_tu); 
%}

eIF_a_s   = hill(par.I_rib_s, snf1, par.theta_rib_s);
eIF_a_tau = hill(tor, par.D_rib_tau, par.theta_rib_tau);

f_to_at   = par.k_tj .* table2array(struct2table(alpha))' .* trna_contri_tc .* hill(ae, par.K_rib_ae, par.n_rib_ae) * eIF_a_s * eIF_a_tau;
f_to_as   = par.k_tj .* table2array(struct2table(alpha))' .* trna_contri_tu .* hill(ae, par.K_rib_ae, par.n_rib_ae) * eIF_a_s * eIF_a_tau;


rib.rf    = r0./(1 + sum(f_to_at) + sum(f_to_as));

rib.rat_j = rib.rf .* f_to_at; 
rib.rat   = sum(rib.rat_j);
rib.ras_j = rib.rf .* f_to_as;
rib.ras   = sum(rib.ras_j);

fat = sum(f_to_at)/(1 + sum(f_to_at));

%% protein synthesis rate



beta = ((par.k_e .* rib.rat_j)./(par.l)); 

prot_syn_rate = sum(par.l .* beta); % aa units 

rib.others = [trna_tc trna_tu ...
trna_contri_tc trna_contri_tu...
hill(ae, par.K_rib_ae, par.n_rib_ae) sum(table2array(struct2table(alpha))')...
eIF_a_s eIF_a_tau...
rib.rf sum(f_to_at) rib.rat sum(f_to_as)...
rib.ras prot_syn_rate 0 fat ...
tRNA.t0  aa_in...
hill(aa_in, par.K_tc_aa, par.n_tc_aa)  1]; 



end 
