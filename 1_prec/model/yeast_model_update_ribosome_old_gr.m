function [dydt, sig, met_reac, prot_syn_rate, beta, alpha, rib, tRNA, eIF_a_s, eIF_a_tau, other_met_reac, g_rate, ribo_rate] = yeast_model_update_ribosome(t, cell_state, par, mutant_type, snf1_vals, jgy_vals, mgl)

%% ------------------------------------------------------------------------
%                               Variables 
%--------------------------------------------------------------------------
%disp(t) 

% metabolites 
met   = cell_state(1:8); 
aa_ex = met(1);  % 1
gl    = met(2);  % 2
pc    = met(3);  % 3
eh    = met(4);  % 4 
aa_in = met(5);  % 5 
ae    = met(6);  % 6 
nh4   = met(7);  % 7 
lp    = met(8);  % 8 

% proteins
prot = cell_state(9:18);
p_r  = prot(1);  % 9  % 1
%{
p_z  = prot(2);  % 10 % 2
p_gy = prot(3);  % 11 % 3
p_fe = prot(4);  % 12 % 4
p_gn = prot(5);  % 13 % 5
p_mt = prot(6);  % 14 % 6
p_as = prot(7);  % 15 % 7
p_at = prot(8);  % 16 % 8
p_lp = prot(9);  % 17 % 9
p_lo = prot(10); % 18 % 10
%}

% free ribosomes 
r0 = cell_state(19);

% number of cells 
n_cells = cell_state(20);

%% ------------------------------------------------------------------------
%                            Signaling network
%--------------------------------------------------------------------------

[sig] = module_signaling(met, par, mutant_type, snf1_vals);

%% ------------------------------------------------------------------------
%                            Metabolic network
%--------------------------------------------------------------------------

[met_reac] = module_met(met, prot, sig.snf1, par, mutant_type, jgy_vals, mgl); 

% fluxes
J_gy     = met_reac.flux(1);
J_fe     = met_reac.flux(2);
J_gn     = met_reac.flux(3);
J_mt     = met_reac.flux(4);
J_as     = met_reac.flux(5);
J_at     = met_reac.flux(6);
J_lp_fe  = met_reac.flux(7); 
J_lp_cit = met_reac.flux(8); 
J_lo     = met_reac.flux(9); 


J_lp     = J_lp_fe + J_lp_cit; 
%% ------------------------------------------------------------------------
%                           Gene expression network
%--------------------------------------------------------------------------

[beta, alpha, rib, tRNA, eIF_a_s, eIF_a_tau] = module_gene_expression_update_ribosome(r0, met, prot, sig.snf1, sig.tor, par, mutant_type);

%% ------------------------------------------------------------------------
%                               growth rate
%--------------------------------------------------------------------------

prot_syn_rate = sum(par.l .* beta); % aa units 
g_rate = prot_syn_rate./par.pro_den;
%g_rate = prot_syn_rate./par.pro_den;

%% ------------------------------------------------------------------------
%                              other fluxes
%--------------------------------------------------------------------------

% J_po = 0; %(par.q_po * g_rate) + (par.d_pm * pc);
J_po = (par.q_po * g_rate) - (J_lp -  J_lo) - J_as;
J_po = max(0,J_po); 
J_eo = (par.q_eo * g_rate) ; %  + (par.d_em * ae) + (par.d_em * J_gn); 
J_no = par.k_no * par.l(1) .* beta(1); 

if nh4 < 1000
J_eo = J_eo*1+ 0.05*par.q_eo; 
end 

other_met_reac.flux.po = J_po; 
other_met_reac.flux.eo = J_eo;  
other_met_reac.flux.no = J_no;  % NH4 to others
%% ------------------------------------------------------------------------
%                              mass balances
%--------------------------------------------------------------------------
J_ATP = (par.q_gy * J_gy) ...
+ (par.q_mt * J_mt) ...
+ (par.q_lo * J_lo) ...
- (par.q_gn * J_gn) ...
- (par.q_as * J_as) ...
- (par.q_lp * J_lp) ...      
- (par.q_p * prot_syn_rate) ...
- (par.q_at * J_at) ...
- J_eo ...
- (g_rate * ae);

if strcmp(mutant_type,'none') ...
|| strcmp(mutant_type,'const_snf1') ...
|| strcmp(mutant_type,'const_tor') ...   
|| strcmp(mutant_type,'const_jgy') ...
|| strcmp(mutant_type,'hxt') ... 
|| strcmp(mutant_type,'hap4_overexpress') ...
|| strcmp(mutant_type,'hap4_underexpress')...
|| strcmp(mutant_type,'addH')...
|| strcmp(mutant_type,'low_gn')...
|| strcmp(mutant_type,'low_mt')...
||strcmp(mutant_type,'low_mt_gn')   

if t > 5
h = 1; 
end 
% lipid syn precursor and ATP and glucose
% lipid oxi precursor and ATP and glucose

%(1.8*J_gy) + (J_gn/2) - J_fe - J_mt - (2*J_as) - J_po - (g_rate * pc) 
met_rate = [ (par.D * (par.aaex_in - aa_ex)) - ((par.V_c * n_cells/par.V_e) * J_at);                             % Aa_e:  extracellular amino acid                                                                                                % Ma
(par.D * (par.gl_in - gl)) - ((par.V_c * n_cells/par.V_e) * (J_gy+ par.n_lp_na*J_lp));              % Gl_ex: glucose                      
(par.n_gy * J_gy) + (par.n_gn * J_gn) + J_lo - J_lp - J_fe - J_mt - (par.n_as * J_as) - J_po - (g_rate * pc);         % Pc:    precursor 
(par.D * (par.eh_in - eh)) + (par.V_c * n_cells/par.V_e * (J_fe - J_gn));                           % Eh:    extracelluar EtOH                                       
par.n_aa * J_at + par.n_aa * J_as - par.n_aa * prot_syn_rate  - (g_rate * aa_in);                      % Aa_i:  amino acid; Jao = 0
J_ATP;                                                                                                 % Ae:    ATP
(par.D * (par.nh4_in - nh4)) - ((par.V_c * n_cells/par.V_e) * (par.n_nh * J_as + J_no));           % nh:    NH4   
(par.D * (par.lp_in - lp)) + (par.V_c * n_cells/par.V_e * (J_lp - J_lo));                          % lp:    lipid    
];       

elseif strcmp(mutant_type,'const_snf1_gl') ...
|| strcmp(mutant_type,'const_tor_gl') ...
|| strcmp(mutant_type,'const_gl') ...
|| strcmp(mutant_type,'const_jgy_gl')

met_rate = [ (par.D * (par.aaex_in - aa_ex)) - ((par.V_c * n_cells/par.V_e) * J_at);                             % Aa_e:  extracellular amino acid                                                                                               
0;                                                                                                     % Gl_ex: glucose                      
(par.n_gy * J_gy) + (par.n_gn * J_gn)  + J_lo - J_lp - J_fe - J_mt - (par.n_as * J_as)  - J_po - (g_rate * pc);        % Pc:    precursor
(par.D * (par.eh_in - eh)) + (par.V_c * n_cells/par.V_e * (J_fe - J_gn));                           % Eh:   extracelluar EtOH                                       
par.n_aa * J_at + par.n_aa * J_as - par.n_aa * prot_syn_rate  - (g_rate * aa_in);                      % Aa_i: amino acid
J_ATP;                                                                                                 % Ae:   ATP
(par.D * (par.nh4_in - nh4)) - ((par.V_c * n_cells/par.V_e) * (par.n_nh * J_as + J_no));           % nh:    NH4   
(par.D * (par.lp_in - lp)) + (par.V_c * n_cells/par.V_e * (J_lp - J_lo));                          % lp:    lipid                    
];      

elseif strcmp(mutant_type,'const_snf1_eh') ...
|| strcmp(mutant_type,'const_tor_eh') ...
|| strcmp(mutant_type,'const_eh') ...
|| strcmp(mutant_type,'const_jgy_eh')

met_rate = [ (par.D * (par.aaex_in - aa_ex)) - ((par.V_c * n_cells/par.V_e) * J_at);                             % Aa_e:  extracellular amino acid                                                                                 
(par.D * (par.gl_in - gl)) - ((par.V_c * n_cells/par.V_e) * (J_gy+ par.n_lp_na*J_lp));              % Gl_ex: glucose                      
(par.n_gy * J_gy) + (par.n_gn * J_gn)  + J_lo - J_lp - J_fe - J_mt - (par.n_as * J_as) - J_po - (g_rate * pc); % Pc:    precursor 
0;                                                                                                     % Eh:    extracelluar EtOH                                       
par.n_aa * J_at + par.n_aa * J_as - par.n_aa * prot_syn_rate  - (g_rate * aa_in);                      % Aa_i:  amino acid; Jao = 0
J_ATP;                                                                                                 % Ae:    ATP
(par.D * (par.nh4_in - nh4)) - ((par.V_c * n_cells/par.V_e) * (par.n_nh * J_as + J_no));           % nh:    NH4   
(par.D * (par.lp_in - lp)) + (par.V_c * n_cells/par.V_e * (J_lp - J_lo));                          % lp:    lipid                   
];         

elseif strcmp(mutant_type,'const_gl_eh_aaex_cell') ...
|| strcmp(mutant_type,'const_snf1_gl_eh_aaex_cell')...
|| strcmp(mutant_type,'const_tor_gl_eh_aaex_cell')...
|| strcmp(mutant_type,'const_gl_eh_aaex_cell_hap_over')...
|| strcmp(mutant_type,'const_gl_eh_aaex_cell_hap_under')...   
|| strcmp(mutant_type,'const_gl_eh_aaex_cell_addH')   

met_rate = [ 0;                                                                                                      % Aa_e:  extracellular amino acid                                                                                                % Ma
0;                                                                                                      % Gl_ex: glucose                      
(par.n_gy * J_gy) + (par.n_gn * J_gn)  + J_lo - J_lp - J_fe - J_mt - (par.n_as * J_as)  - J_po - (g_rate * pc);          % Pc:    precursor
0;                                                                                                      % Eh:    extracelluar EtOH                                       
par.n_aa * J_at + par.n_aa * J_as - par.n_aa * prot_syn_rate - (g_rate * aa_in);                      % Aa_i:  amino acid
J_ATP;                                                                                                  % Ae:    ATP
0;                                                                                                      % nh:    NH4   
(par.D * (par.lp_in - lp)) + 0*(par.V_c * n_cells/par.V_e * (J_lp - J_lo));                           % lp:    lipid                   
];                    
end 

prot_rate = beta - (g_rate * prot); % protein units 
ribo_rate = par.k_ro_plus * (p_r - (r0 * par.n_ro)) - (par.k_ro_minus * r0) - (g_rate * r0); % par.n_ro = 1

if strcmp(mutant_type,'none') ...
|| strcmp(mutant_type,'const_snf1') ...
|| strcmp(mutant_type,'const_snf1_gl') ...
|| strcmp(mutant_type,'const_tor') ...
|| strcmp(mutant_type,'const_tor_gl') ...        
|| strcmp(mutant_type,'const_gl') ...
|| strcmp(mutant_type,'const_jgy') ...
|| strcmp(mutant_type,'hxt') ...
|| strcmp(mutant_type,'hap4_overexpress') ...
|| strcmp(mutant_type,'hap4_underexpress')...
|| strcmp(mutant_type,'addH')...
|| strcmp(mutant_type,'const_snf1_eh')...
|| strcmp(mutant_type,'const_tor_eh') ...
|| strcmp(mutant_type,'const_eh') ...
|| strcmp(mutant_type,'const_jgy_eh')...    
|| strcmp(mutant_type,'low_gn')...
|| strcmp(mutant_type,'low_mt')...
||strcmp(mutant_type,'low_mt_gn')          

cell_rate = - (par.D * n_cells) + (g_rate * n_cells);

elseif strcmp(mutant_type,'const_gl_eh_aaex_cell') ...
|| strcmp(mutant_type,'const_snf1_gl_eh_aaex_cell')...
|| strcmp(mutant_type,'const_tor_gl_eh_aaex_cell')...
|| strcmp(mutant_type,'const_gl_eh_aaex_cell_hap_over')...
|| strcmp(mutant_type,'const_gl_eh_aaex_cell_hap_under')...  
|| strcmp(mutant_type,'const_gl_eh_aaex_cell_addH')     
cell_rate = 0; 

end 

%% ------------------------------------------------------------------------
%                                  outputs
%--------------------------------------------------------------------------

if  strcmp(mutant_type,'low_gn')...
||strcmp(mutant_type,'low_mt_gn')
prot_rate(5) = 0; % Egn concentration does not change
end 


if  strcmp(mutant_type,'low_mt')...
||strcmp(mutant_type,'low_mt_gn')    
prot_rate(6) = 0; % Emt concentration does not change
end 

dydt = [met_rate;
prot_rate;
ribo_rate;
cell_rate];







end 
