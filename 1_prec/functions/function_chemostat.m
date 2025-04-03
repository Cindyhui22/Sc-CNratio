function [y_steady, met_reac_steady, other_met_reac_steady, sig_steady, prot_syn_steady, rib_steady, g_rate_steady, ribo_rate_steady] = function_chemostat(par, t_start, t_final, aa_ex,nh_ex, gl, eh, cells, mutant, snf1_vals, jgy_vals, mgl, new_cell_state)

%
tic
%disp('Start chemostat simulation...')

%disp('Calculating initiation conditions...')

% initial conditions 
% metabolites
met.ex_amino_acids = aa_ex;       % 1
met.glucose        = gl;          % 2
met.precursor      = 1.00E+02;    % 3
met.ethanol        = eh;          % 4  
met.in_amino_acids = 5*10^2;      % 5  
met.atp            = 3.00E+03;    % 6 
met.nh             = nh_ex;       % 7
met.lipid          = 0;           % 8 
%num.met = numel(fieldnames(met));

% proteins 
prot.r  = 0.15      * par.pro_den / par.l(1);     % 9   1 % converted amino acid fraction to protein fraction
prot.z  = 0.524     * par.pro_den / par.l(2);    % 10  2
prot.gy = 0.0567    * par.pro_den / par.l(3);    % 11  3
prot.fe = 0.0200    * par.pro_den / par.l(4);    % 12  4
prot.gn = 0.0540    * par.pro_den / par.l(5);    % 13  5
prot.mt = 0.1137    * par.pro_den / par.l(6);    % 14  6
prot.as = 0.0707    * par.pro_den / par.l(7);    % 15  7
prot.at = 0.0066    * par.pro_den / par.l(8);    % 16  8
prot.lp = 0.0101    * par.pro_den / par.l(9);    % 17  9 %_% par.pro_den / par.l(7);
prot.lo = 0.0102    * par.pro_den / par.l(10);    % 18  10 %_% par.pro_den / par.l(7);

% total ribosomes 
R0 = 7.8; % 19

% temp cell state

cell_state = [table2array(struct2table(met))';  % 1 - 8
table2array(struct2table(prot))'; % 9 - 18
R0;                               % 19
cells;                            % 20
];


if strcmp( new_cell_state.overwrite, 'True' )
cell_state = new_cell_state.values; 
end 

% load WT_y_steady_10gL.mat; 
% cell_state = WT_y_steady_10gL(1,:);
%WT_y_steady_10gL =[];


%load ss_cell_state.mat 
%cell_state = cell_state_steady;

% set ode options
%Rel_tol  = 1.0E-09; 
%Abs_tol  = 1.0E-12; 
Rel_tol  = 1.0E-03; 
Abs_tol  = 1.0E-06; 
%max_step = 100;
tstart = tic;
%
options  = odeset('RelTol',Rel_tol, ...
'AbsTol',Abs_tol, ...
'NonNegative',(1:length(cell_state)),...
'Events',@(t,y) myevent(t,y,tstart));    
%}

% get initial cell state
%{
cell_state = get_cell_state_update_ribosome(cell_state_temp, par, mutant, snf1_vals, jgy_vals, mgl);
disp('Calculate initial conditions for chemostat simulation: done')  
%}
%disp('Running chemostat simulation...')    
[t, y,te,ye,ie] = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals, jgy_vals, mgl), ...
[t_start:t_final], ...
cell_state, ...
options); 

y = real(y);           
y_steady = y(end,:); 

% get intermediate values      
%disp('Calculating intermediate values...')    
%for k = 1:length(t)
%[dydt_t, sig_t, met_reac_t, prot_syn_rate_t, beta_t, alpha_t, ra_t, rat_t, r0_t, tRNA_t, eIF_a_t, eIF_a_s_t, eIF_a_tau_t, other_met_reac_t] = yeast_model(t_batch, y_batch', par, mutant_type);
%[dydt_t, sig_t, met_reac_t, prot_syn_rate_t, beta_t, alpha_t, rib_t, tRNA_t, eIF_a_s_t, eIF_a_tau_t, other_met_reac_t, g_rate_t, ribo_rate_t] = yeast_model(t_batch, y_batch', par, mutant_type);
[~, sig_t, met_reac_t, prot_syn_rate_t, beta_t, alpha_t, rib_t, tRNA_t, eIF_a_s_t, eIF_a_tau_t, other_met_reac_t, g_rate_t, ribo_rate_t] = yeast_model_update_ribosome(t(end), y(end,:)', par, mutant, snf1_vals, jgy_vals, mgl);

met_reac_steady.prot      = real(met_reac_t.prot)';
met_reac_steady.substrate = real(met_reac_t.substrate)';
met_reac_steady.atp       = real(met_reac_t.atp)';
met_reac_steady.sig       = real(met_reac_t.sig)';
met_reac_steady.flux      = real(met_reac_t.flux)';
met_reac_steady.snf1      = real(met_reac_t.snf1)';
met_reac_steady.pc        = real(met_reac_t.pc)';

other_met_reac_steady.flux.po  = real(other_met_reac_t.flux.po)';
other_met_reac_steady.flux.eo  = real(other_met_reac_t.flux.eo)';
other_met_reac_steady.flux.no  = real(other_met_reac_t.flux.no)';

prot_syn_steady.alpha         = real(table2array(struct2table(alpha_t))); 
prot_syn_steady.tc            = real(tRNA_t.tc); 
prot_syn_steady.tu            = real(tRNA_t.tu);
prot_syn_steady.t0            = real(tRNA_t.t0); 
%prot_syn.eIF_a                = real(eIF_a_t)'; 
prot_syn_steady.beta          = real(beta_t)';
prot_syn_steady.prot_syn_rate = real(prot_syn_rate_t); 
prot_syn_steady.eIF_a_s       = real(eIF_a_s_t);
prot_syn_steady.eIF_a_tau     = real(eIF_a_tau_t);
prot_syn_steady.prot_syn_rate = real(prot_syn_rate_t);
rib_steady.others             = real(rib_t.others); 

sig_steady.snf1 = real(sig_t.snf1);
sig_steady.tor  = real(sig_t.tor); 

rib_steady.rf      = real(rib_t.rf); 
%rib.raf            = real(rib_t.raf);
%rib.ri             = real(rib_t.ri);
rib_steady.rat_j   = real(rib_t.rat_j)';
rib_steady.rat     = real(rib_t.rat);
rib_steady.ras_j   = real(rib_t.ras_j)';
rib_steady.ras     = real(rib_t.ras); 

g_rate_steady      = real(g_rate_t); 
ribo_rate_steady   = real(ribo_rate_t); 
%end 

end 
%
function [values, isterminal, direction] = myevent(t, y,tstart)
%  Don't let t cross zero...a dummy "event" to illustrate how 
%  one might handle other events in conjunction with the time
%  constraint.  Not necessary, but I put it in just in case.
%  Don't let integration go for more than 5 seconds.
values = toc(tstart) < 20; % The value that we want to be zero
isterminal = 1;           % Halt integration 
direction  = 0;           % The zero can be approached from either direction
end 
%}
