function [label] = plot_labels

label.met = {{'amino'; 'acids_{ex}'};
'glucose_{ex}';
'precursor';
'ethanol';
{'amino'; 'acids_{in}'};
'atp';
'nh4';
'lipid'; 
%_% 'sc'
};

label.prot = {'R';
'Z';
'gy';
'fe';
'gn';
'mt';
'as';
'at';
'lp';
'lo';
%_% 'sp';
%_% 'sd'
};   

label.rib_cells = {'R_{0}';
'cells (g/L)'};

label.rib = {'R_{af}';
'R_{i}';
'R_{at_j}';
'R_{at}';
'R_{as_j}';
'R_{as}';
'R_{f}'};              

label.met_reac.flux = {'J_{gy}';
'J_{fe}';
'J_{gn}';
'J_{mt}';
'J_{as}';
'J_{at}';
'J_{lp,fe}';
%_% 'J_{lp,ct}';
'J_{lo}';
%_% 'J_{sp}';
%_% 'J_{sd}';
};   

label.met_reac.prot = {'E_{gy}';
'E_{fe}';
'E_{gn}';
'E_{mt}';
'E_{as}';
'E_{at}';
'E_{lp}'; % J_{lp,fe};
%_% 'E_{lp}'; % J_{lp,cit};
'E_{lo}';
%_% 'E_{sp}';
%_% 'E_{sd}';
};        
% flux:
%         gy: 1
%         fe: 2
%         gn: 3
%         mt: 4
%         as: 5
%         at: 6
%      lp_fe: 7
%     lp_cit: 8
%         lo: 9   

label.met_reac.substrate = {{'k*[G_{l}]/' , '(K+[G_{l}])'};
{'k*[P_{c}]/' , '(K+[P_{c}])'};
{'k*[E_{h}]/' , '(K+[E_{h}])'};
{'k*[P_{c}]/' , '(K+[P_{c}])'};
{'[P_{c}][NH4]/' , '(K+[P_{c}])(K+[NH4])'};
{'k*[A_{a}^{e}]/' , '(K+[A_{a}^{e}])'}
{'k*[P_{c}]/' , '(K+[P_{c}])'};
%_% {'k*[P_{c}]/' , '(K+[P_{c}])'};
{'k*[L_{p}]/' , '(K+[L_{p}])'}; 
%_% {'k*[P_{c}]/' , '(K+[P_{c}])'};
%_% {'k*[S_{c}]/' , '(K+[S_{c}])'};
};  

label.met_reac.atp = {{'(I^{n})/' , '(I^{n} + [A_{e}]^{n})'};
'none';
{'([A_{e}]^{n})/' , '(I^{n} + [A_{e}]^{n})'};
{'(I^{n})/' , '(I^{n} + [A_{e}]^{n})'};
{'([A_{e}]^{n})/' , '(I^{n} + [A_{e}]^{n})'};
{'([A_{e}]^{n})/' , '(I^{n} + [A_{e}]^{n})'};
{'([A_{e}]^{n})/' , '(I^{n} + [A_{e}]^{n})'};
'none';
%_% 'none';
'none';
%_% 'Ae'; 
%_% 'Ae'
};    

%{                        
label.met_reac.atp = {{'(\eta_{gy,ae} + I_{gy,ae}^{n_{gy,ae}})/' , '(I_{gy,ae}^{n_{gy,ae}} + [A_{e}]^{n_{gy,ae}})'};
'none';
{'(\eta_{gn,ae} + [A_{e}]^{n_{gn,ae}})/'   , '(I_{gn,ae}^{n_{gn,ae}} + [A_{e}]^{n_{gn,ae}})'};
{'(\eta_{mt,ae} + [A_{e}]^{n_{mt,ae}})/'   , '(I_{mt,ae}^{n_{mt,ae}} + [A_{e}]^{n_{mt,ae}})'};
{'(\eta_{as,ae} + [A_{e}]^{n_{as,ae}})/'   , '(I_{as,ae}^{n_{as,ae}} + [A_{e}]^{n_{as,ae}})'};
{'(\eta_{sp,ae} + [A_{e}]^{n_{sp,ae}})/'   , '(I_{sp,ae}^{n_{sp,ae}} + [A_{e}]^{n_{sp,ae}})'};
'none';
{'(\eta_{at,ae} + [A_{e}]^{n_{at,ae}})/'   , '(I_{at,ae}^{n_{at,ae}} + [A_{e}]^{n_{at,ae}})'}};                         
%}

% label.met_reac.sig = {{'(I^{n})/' , '(I_{gy,s}^{n_{gy,s}} + [s^{*}]^{n_{gy,s}})'};
%                       'none';
%                       {'([s^{*}]^{n})/' , '(I^{n} + [s^{*}]^{n})'};
%                       {'([s^{*}]^{n})/' , '(I^{n} + [s^{*}]^{n})'};
%                       {'(I^{n} )/' , '(I^{n} + [A_{a}^{i}]^{n})'};
%                       {'(I^{n} )/' , '(I^{n} + [A_{a}^{i}]^{n})'};
%                       {'(I^{n} )/' , '(I^{n} + [s^{*}]^{n})'};
%                       {'(I^{n} )/' , '(I^{n} + [s^{*}]^{n})AA'};
%                       'none'}; 
% 
% 
% label.met_reac.sig = {{'(I_{gy,s}^{n_{gy,s}})/'    , '(I_{gy,s}^{n_{gy,s}} + [s^{*}]^{n_{gy,s}})'};
%                       'none';
%                       {'([s^{*}]^{n_{gn,s}})/'     , '(I_{gn,s}^{n_{gn,s}} + [s^{*}]^{n_{gn,s}})'};
%                       {'( [s^{*}]^{n_{mt,s}})/'     , '(I_{mt,s}^{n_{mt,s}} + [s^{*}]^{n_{mt,s}})'};
%                       {'(I_{as,ai}^{n_{as,ai}})/' , '(I_{as,ai}^{n_{as,ai}} + [A_{a}^{i}]^{n_{as,s}})'};
%                       {'(I_{at,ai}^{n_{at,ai}})/' , '(I_{at,ai}^{n_{at,ai}} + [A_{a}^{i}]^{n_{at,s}})'};
%                       {'(I^{n} )/' , '(I^{n} + [s^{*}]^{n})'};
%                       {'(I^{n} )/' , '(I^{n} + [s^{*}]^{n})AA'};
%                       'none'}; 


label.met_reac.sig = {{'[Pc],[s^{*}] term'};
'none';
{'[s^{*}] term'};
{'[Pc],[s^{*}] term'};
{'[Pc],[Aa] term'};
{'(I_{at,ai}^{n_{at,ai}})/' , '(I_{at,ai}^{n_{at,ai}} + [A_{a}^{i}]^{n_{at,s}})'},
{'[s^{*}]'};
%_% {'[s^{*}], [AA] term'};
{'none'}; 
%_% {'[s^{*}]'};
%_% {'[s^{*}]'};
}; 

%%
label.met_reac.snf1 = {{'[s^{*}] term'};
'none';
{'[s^{*}] term'};
{'none'};
{'none'};
{'(I_{at,ai}^{n_{at,ai}})/' , '(I_{at,ai}^{n_{at,ai}} + [A_{a}^{i}]^{n_{at,s}})'};
{'[s^{*}] term'};
{'[s^{*}] term'};
{'none'};}; 

label.met_reac.pc = {{'[Pc] term'};
'none';
{'none'};
{'[Pc] term'};
{'[Aa] term'};
{'(I_{at,ai}^{n_{at,ai}})/' , '(I_{at,ai}^{n_{at,ai}} + [A_{a}^{i}]^{n_{at,s}})'};
{'none'};
{'[Aa] term'}; 
{'none'}; }; 

%}

label.prot_syn.alpha = {'\alpha_{r}';
'\alpha_{z}';
'\alpha_{gy}';
'\alpha_{fe}';
'\alpha_{gn}';
'\alpha_{mt}';
'\alpha_{as}';
'\alpha_{at}'};   

label.prot_syn.beta  = {'\beta_{r}';
'\beta_{z}';
'\beta_{gy}';
'\beta_{fe}';
'\beta_{gn}';
'\beta_{mt}';
'\beta_{as}';
'\beta_{at}'};                      

label.prot_syn.prot_syn_rate = {'l_{j} \beta_{j}'};   
label.prot_syn.eIF_a_s       = {'K^{n}/','(K^{n} + [s^{*}]^{n})'}; 
label.prot_syn.eIF_a_tau     = {'[\tau^{*}]^{n})/','(K^{n} + [\tau^{*}]^{n})'};

label.prot_syn.tc = {'t_{c}'};
%label.prot_syn.eIF_a = {'eIF_{a}'};

label.rat_j = {'R_{at,r}';
'R_{at,z}';
'R_{at,gy}';
'R_{at,fe}';
'R_{at,gn}';
'R_{at,mt}';
'R_{at,as}';
'R_{at,at}'};

label.sig.snf1 = {'SNF1 activity'};
label.sig.tor  = {'TORC1 activity'};

label.others = {'trna_tc';
'trna_tu';
'trna_contri_tc';
'trna_contri_tu';
'atp_contri';
'sum \alpha';
'eIF_a_s'; 
'eIF_a_tau';
'Rf';
'f_to_at';
'Rat';
'f_to_as';                    
'Ras';
'prot_sys_rate'; 
'g_rate';
'fat';
't_0';
'aa'; 
'aa term';
'\rho term'};

label.atp_flux = {'q_{gy} J_{gy}';
'q_{mt} J_{mt}';
'q_{gn} J_{gn}';
'q_{as} J_{as}';
'q_{p} l_{j}\beta_{j}';
'J_{eo}';
'q_{lo} J_{lo}'; 
'q_{lp} J_{lp}'; 
};  

label.atp_flux_chem = {'q_{gy} J_{gy}';
'q_{mt} J_{mt}';
'q_{gn} J_{gn}';
'q_{as} J_{as}';
'q_{p} l_{j}\beta_{j}';
'J_{eo}';
'q_{lo} J_{lo}'; 
'q_{lp} J_{lp}'; 
};                

end 
