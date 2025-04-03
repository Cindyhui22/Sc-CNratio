function [dydt, sig, met_reac, prot_syn_rate, beta, alpha, rib, tRNA, eIF_a_s, eIF_a_tau, other_met_reac, g_rate, ribo_rate] = yeast_model_update_ribosome(t, cell_state, par, mutant_type, snf1_vals, jgy_vals, mgl)

%% ------------------------------------------------------------------------
%                               Variables 
%--------------------------------------------------------------------------
%disp(t) 

% metabolites 
met   = cell_state(1:8); %_% cell_state(1:9); 
aa_ex = met(1);  % 1
gl    = met(2);  % 2
pc    = met(3);  % 3
eh    = met(4);  % 4 
aa_in = met(5);  % 5 
ae    = met(6);  % 6 
nh4   = met(7);  % 7 
lp    = met(8);  % 8 
%_% sc    = met(9); 

% proteins
prot = cell_state(9:18); %_%cell_state(10:21);
p_r  = prot(1);  % 10  % 1
%{
p_z  = prot(2);  % 11 % 2
p_gy = prot(3);  % 12 % 3
p_fe = prot(4);  % 13 % 4
p_gn = prot(5);  % 14 % 5
p_mt = prot(6);  % 15 % 6
p_as = prot(7);  % 16 % 7
p_at = prot(8);  % 17 % 8
p_lp = prot(9);  % 18 % 9
p_lo = prot(10); % 19 % 10
p_sp = prot(11); % 20 % 11
p_sd = prot(12); % 21 % 12
%}

% free ribosomes 
r0 = cell_state(19); %_% cell_state(22);

% number of cells 
n_cells = cell_state(20); %_% cell_state(23);

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
%_% J_lp_cit = met_reac.flux(8); 
J_lo     = met_reac.flux(8); %_% met_reac.flux(9); 
%_% J_sp     = met_reac.flux(10); 
%_% J_sd     = met_reac.flux(11); 

J_lp     = J_lp_fe; %_% + J_lp_cit; 




%% ------------------------------------------------------------------------
%                           Gene expression network
%--------------------------------------------------------------------------

[beta, alpha, rib, tRNA, eIF_a_s, eIF_a_tau] = module_gene_expression_update_ribosome(r0, met, prot, sig.snf1, sig.tor, par, mutant_type);

%% ------------------------------------------------------------------------
%                               growth rate
%--------------------------------------------------------------------------



prot_syn_rate = sum(par.l .* beta); % aa units 



% no sc version
%_% g_rate = (10^-6)*(110* par.a_n *prot_syn_rate + 32*par.a_c*(J_lp-J_lo))/par.rho_cell;
% g_rate = (par.a_n *prot_syn_rate + par.a_c*(J_lp-J_lo))/par.rho_cell; %_% 
g_rate = (par.a_n *prot_syn_rate + par.a_c*(par.n_lp * J_lp - par.n_lo * J_lo))/par.rho_cell; %-% 
%g_rate = (10^-6)*(110*1.81*prot_syn_rate + par.fake_sc*31*4.66*(J_lp-J_lo))/par.rho_cell;


% has sc version
%g_rate = (10^-6)*(110*par.a_n*prot_syn_rate + par.fake_sc*181*par.a_c*(J_sp-J_sd))/par.rho_cell;
v_n_cells = n_cells/par.rho_cell; % n_cells is DW(g/L) for now 

% hold protein constant
%g_rate = prot_syn_rate./par.pro_den;
%g_rate = prot_syn_rate./par.pro_den;

%% ------------------------------------------------------------------------
%                              other fluxes
%--------------------------------------------------------------------------

% J_po = 0; %(par.q_po * g_rate) + (par.d_pm * pc);


%_% J_po = (par.q_po * g_rate) - (J_lp -  J_lo) - J_as; 
%_% J_po = max(0, J_po); 
J_po = par.q_po * g_rate;
J_eo = (par.q_eo * g_rate) + (par.d_em * ae) + (par.d_em * J_gn); 
%_% J_no = par.k_no * par.l(1) .* beta(1); 
J_no = par.q_no * par.q_po * g_rate;

%_% J_en = 40000*(100000/(nh4+1000));  % bigger too smallgf 
J_en = par.k_en * (par.I_en/(nh4+par.I_en)) * g_rate; %_%


J_eo = J_eo + J_en;  % sc %_%


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

met_rate = [ (par.D * (par.aaex_in - aa_ex)) - (v_n_cells * J_at);                % Aa_e:  extracellular amino acid                                                                                                % Ma
(par.D * (par.gl_in - gl)) - (v_n_cells * (J_gy+ par.n_lp_na*J_lp));              % Gl_ex: glucose                      
(par.n_gy * J_gy) + (par.n_gn * J_gn) + J_lo - J_lp - J_fe - J_mt - (par.n_as * J_as) - J_po - (g_rate * pc); % Pc:    precursor 
(par.D * (par.eh_in - eh)) + (v_n_cells * (J_fe - J_gn));                           % Eh:    extracelluar EtOH                                       
par.n_aa * J_at + par.n_aa * J_as - par.n_aa * prot_syn_rate  - (g_rate * aa_in);  % Aa_i:  amino acid; Jao = 0
J_ATP;                                                                             % Ae:    ATP
(par.D * (par.nh4_in - nh4)) - (v_n_cells * (par.n_nh * J_as + J_no));             % nh:    NH4 
(par.n_lp * J_lp - par.n_lo * J_lo)- (g_rate * lp);                                % lp:    lipid  %_%
%_% (par.D * (par.lp_in - lp)) + (v_n_cells * (J_lp - J_lo));                      % lp:    lipid  
];       

elseif strcmp(mutant_type,'const_snf1_gl') ...
|| strcmp(mutant_type,'const_tor_gl') ...
|| strcmp(mutant_type,'const_gl') ...
|| strcmp(mutant_type,'const_jgy_gl')

met_rate = [ (par.D * (par.aaex_in - aa_ex)) - (v_n_cells * J_at);                               % Aa_e:  extracellular amino acid                                                                                               
0;                                                                                 % Gl_ex: glucose                      
(par.n_gy * J_gy) + (par.n_gn * J_gn)  + J_lo - J_lp - J_fe - J_mt - (par.n_as * J_as)  - J_po - (g_rate * pc);  % Pc:    precursor
(par.D * (par.eh_in - eh)) + (v_n_cells * (J_fe - J_gn));                           % Eh:   extracelluar EtOH                                       
par.n_aa * J_at + par.n_aa * J_as - par.n_aa * prot_syn_rate  - (g_rate * aa_in);  % Aa_i: amino acid
J_ATP;                                                                             % Ae:   ATP
(par.D * (par.nh4_in - nh4)) - (v_n_cells * (par.n_nh * J_as + J_no));             % nh:    NH4   
(par.n_lp * J_lp - par.n_lo * J_lo)- (g_rate * lp);                                                      % lp:    lipid  %_%
%_% (par.D * (par.lp_in - lp)) + (v_n_cells * (J_lp - J_lo));                      % lp:    lipid                        
];      

elseif strcmp(mutant_type,'const_snf1_eh') ...
|| strcmp(mutant_type,'const_tor_eh') ...
|| strcmp(mutant_type,'const_eh') ...
|| strcmp(mutant_type,'const_jgy_eh')

met_rate = [ (par.D * (par.aaex_in - aa_ex)) - (v_n_cells * J_at);                              % Aa_e:  extracellular amino acid                                                                                 
(par.D * (par.gl_in - gl)) - (v_n_cells * (J_gy+ par.n_lp_na*J_lp));               % Gl_ex: glucose                      
(par.n_gy * J_gy) + (par.n_gn * J_gn)  + J_lo - J_lp - J_fe - J_mt - (par.n_as * J_as) - J_po - (g_rate * pc); % Pc:    precursor 
0;                                                                                 % Eh:    extracelluar EtOH                                       
par.n_aa * J_at + par.n_aa * J_as - par.n_aa * prot_syn_rate  - (g_rate * aa_in);  % Aa_i:  amino acid; Jao = 0
J_ATP;                                                                             % Ae:    ATP
(par.D * (par.nh4_in - nh4)) - (v_n_cells * (par.n_nh * J_as + J_no));             % nh:    NH4   
(par.n_lp * J_lp - par.n_lo * J_lo)- (g_rate * lp);                                                      % lp:    lipid  %_%
%_% (par.D * (par.lp_in - lp)) + (v_n_cells * (J_lp - J_lo));                      % lp:    lipid  
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
par.n_aa * J_at + par.n_aa * J_as - par.n_aa * prot_syn_rate - (g_rate * aa_in);                        % Aa_i:  amino acid
J_ATP;                                                                                                  % Ae:    ATP
0;                                                                                                      % nh:    NH4   
(par.n_lp * J_lp - par.n_lo * J_lo)- (g_rate * lp);                                                      % lp:    lipid  %_%
%_% (par.D * (par.lp_in - lp)) + (v_n_cells * (J_lp - J_lo));                      % lp:    lipid 
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
