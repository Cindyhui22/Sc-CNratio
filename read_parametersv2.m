function par = read_parametersv2
% --------------------------------------------------

% --------------------------------------------------
%% metabolic network
%  total snf1 and tor (reorganize later)
par.s_tot            = 0.22; %_% 
par.tau_tot          = 0.09; %_% 

% rate                 
par.k_gy         =  8.5279e+04 * 1.1; %_% newly modified
par.k_fe         =  2.7555e+05; 
par.k_as         =  5.0255e+03;
par.k_mt         =  3.2010e+04;
par.k_gn         =  3.9203e+04; 
par.k_at         =  3.5254e+09;

% substrate contribution 
par.K_mgl        = 3200;
par.K_gy         = 4300; 
par.K_fe         = 6000; 
par.n_fe_p       = 4; 
par.K_as_p       = 390;
par.K_mt         = 390; % * 2
par.K_gn_eh      = 13000;
par.K_at_aa      = 59;

%% added regulation
par.K_gy_gl     = 5.6; 

par.n_mt_p      = 2;  
par.n_as_pc      = 2; 

par.I_gy_pc      = 20000; % high Ki for T6P
par.n_gy_pc      = 1; 


% let PkA and pc repression contribute to 30% of Emt activity 
par.I_mt_pc    = 10000;
par.n_mt_pc    = 1;


%%
% ATP contribution 
par.I_gy_ae         = 3000;
par.n_gy_ae         = 1.00E+00;


par.K_as_ae         = 1.00E+02;

par.I_mt_ae         = 1100;
par.n_mt_ae         = 2; 

par.D_at_ae         = 6.6825e+03; % I think K_at_ae, n_at_ae value should be samilar as K_rib_ae.
par.n_at_ae         = 4; 


% signaling contribution 
par.I_gy_s          = 1.2*par.s_tot; 0.95*par.s_tot;
par.n_gy_s          = 1; 

par.I_as_ai         = 2000;
par.n_as_ac         = 1; 

par.D_gn_s          = 0.2*par.s_tot; 
par.n_gn_s          = 2; 

par.I_at_ai         = 1000;
par.n_at_ai         = 2;


%% signaling network
% snf1 
par.K_s_pc            = 2.2894e+03; 
par.theta_s_pc        = 2; 
par.s_basal           = 0.0200;

% tor 
par.K_tau_aai         = 5300;
par.theta_tau_aai     = 0.5;
par.tau_basal         = 0.2; 

%% gene expression network

% transcription factors 
par.ep_r        = 0;
par.xi_r_tau    = 1;
par.xi_r_s      = 0;
par.w_r_tau     = 0.8*par.tau_tot; 
par.w_r_s       = 0.4*par.s_tot;
par.theta_r_tau = 1;
par.theta_r_s   = 1; 
par.a0_r        = 0; 

par.ep_gn       = 0;
par.w_gn_s      = 0.3*par.s_tot;
par.xi_gn_s     = 1;
par.theta_gn_s  = 3;
par.a0_gn       = 0; 

par.ep_mt       = 0;
par.xi_mt_s     = 1;
par.w_mt_s      = 0.4*par.s_tot; 
par.theta_mt_s  = 1.2000;
par.a0_mt        = 0; 

par.ep_as        = 1;
par.xi_as_tau    = 0;
par.xi_as_s      = 0;
par.w_as_tau     = 0.5*par.tau_tot;
par.w_as_s       = 0.65*par.s_tot; 
par.theta_as_tau = 2;
par.theta_as_s   = 2; 
par.a0_as         = 0; 

%_% newly added; 
par.ep_lp       = 1;
%_% par.xi_lp_s     = 1; 
par.xi_lp_tau   = 0;
%_% par.w_lp_s      = 0.2*par.s_tot;
par.w_lp_tau    = 0.5*par.tau_tot;
%_% par.theta_lp_s  = 1;
par.theta_lp_tau= 1;
par.a0_lp       = 0; 


par.ep_lo        = 0;
par.xi_lo_s      = 1;
par.xi_lo_tau    = 1;%_% 0;
par.w_lo_s       = 1*par.s_tot*0.2;
par.w_lo_tau     = 0.5*par.tau_tot; %*0.5
par.theta_lo_s   = 3; 
par.theta_lo_tau = 1;
par.a0_lo        = 0; 


%% tRNA
par.mt             = 11;
par.K_tc_aa        = 1600;
par.n_tc_aa        = 2; 


par.K_tc_ae        = 100; 
par.n_tc_ae        = 1;


par.K_rib_c       = (3/4)*20; % times 20 amino acids, I think Chen forgot to do it?
par.K_rib_u       = 15;

%% protein synthesis rate   

par.k_tj =     [    4.4157; %_%  *1.1; % 1: r
2.9834*0.5;%_% *0.8; % 2: z ; %_% newly modified
0.3229  % 3: gy
0.1139; % 4: fe
0.2256; % 5: gn
0.8526; % 6: mt
1.5841; % 7: as
0.0375; % 8: at
1.4595*0.5*(528/465);%_% *0.7; % 9: lp   %_% newly added
0.1390*0.1*(412/465);%_% *0.3; % 10: lo  %_% newly added                
]; 
% par.k_tj(num_prot.lp) =  1.4595*0.7;
% par.k_tj(num_prot.lo) =  0.1390; 
% par.k_tj(num_prot.z) =  par.k_tj(num_prot.z)*0.8;    

par.k_e          =     36720;


par.l            =     [ 12467;  % 1: r
471;   % 2:  z
421;   % 3:  gy
511;   % 4:  fe
698;   % 5:  gn
434;   % 6:  mt
550;   % 7:  as
566;   % 8:  at 
528;%465;   % 9:  lp %_% newly added 528
412;%465;   % 10: lo %_% newly added 412
];                       



% Protein_sorting_2020_10_29_update_2021_10/protein_length/protein_length.xlsx                     
%% ribosome partition 

par.I_rib_s         = 3.6* par.s_tot;
par.theta_rib_s     = 1;
par.D_rib_tau       = 0.1* par.tau_tot;
par.theta_rib_tau   = 1;

par.K_rib_ae       = 30;
par.n_rib_ae       = 1;

par.k_ro_plus          = 3.6*10^4;
par.k_ro_minus          = 9417; 
par.n_ro          = 1;


%% other fluxes 
par.d_em         = 0;

%% ATP stoichiometry 

par.q_gy        = 1.4; % =2*0.7, 0.3 goes to p1--> others 
par.q_at        = 1.00E+00;
par.q_as        = 2; 
par.q_p         = 4.00E+00;
par.q_eo        = 12480000*0.8; %_% newly modified
par.q_gn        = 1; 
par.q_mt        = 9; 
par.n_aa        = 1/20; 


%% reaction stoichiometry
par.n_gy        = 2;
par.n_gn        = 0.5; 
par.n_as        = 2; 

%% cell parameters
par.V_e         = 1;                           % L
par.V_c       = 4.20E-14;                    % Volume of 1 cell in liters

%% experiment parameters

par.D            = 0; % dilution rate  - set in run files
par.gl_in        = 0; % influx glucose - set in run files
par.eh_in        = 0;
par.aaex_in      = 0; 

par.pro_den      = 800000; 

%_% used for optimization 
par.c_3AT   = 0; % 3ATP concentration %_% 
par.K_3AT   = 30; % unit ng/mL %_% 

%
%% expand par.l and par.k_tj in order to perturb each of them in sensitivy analysis; 
% paramters below are not called in 
par.l_r	    =	12467 ; % 1
par.l_z	    =	471	;   % 2
par.l_gy	=	421	;   % 3
par.l_fe	=	511	;   % 4
par.l_gn	=	698	;   % 5
par.l_mt	=	434	;   % 6
par.l_as	=	550	;   % 7
par.l_at	=	566 ;   % 8
par.l_lp	=	465	;   % 9
par.l_lo	=	465 ;   % 10


par.K_t_r	=	par.k_tj(1)	;
par.K_t_z	=	par.k_tj(2)	;
par.K_t_gy	=	par.k_tj(3)	;
par.K_t_fe	=	par.k_tj(4)	;
par.K_t_gn	=	par.k_tj(5)	;
par.K_t_mt	=	par.k_tj(6)	;
par.K_t_as	=	par.k_tj(7)	;
par.K_t_at	=	par.k_tj(8) ;
par.K_t_lp	=	par.k_tj(9) ;
par.K_t_lo	=	par.k_tj(10) ;
%}





%% newly added 

par.k_lo          = 100; %100/(836.5/31.9);  %_%-
par.k_lp_fe       = 76.190476190476190*3; %76.190476190476190*3/(836.5/31.9); %_%-
par.n_lp          = 1/26.4; %-%
par.n_lo          = 1/26.4; %-%
par.K_lp_fe_p        = 500;
%_% par.n_lp_fe_p        = 1;

par.K_lo_lp         = 1.89;%50/(842.3/31.9); %_%-


par.I_lp_s          = 5; 
%_% par.n_lp_s          = 1; 

par.K_as_nh        = 1000; 
par.q_po            = 2.053e06; %_% 1600000*2; % = par.pro_den*4; 
par.q_no            = 0.0156;   %_%


par.nh4_in         = 0;
par.n_nh           = 1.8; %_% 1.4;

%_% par.k_no           = 0.1; %_% 
par.lp_in          = 0; 

par.n_lp_na        = 1.8/12; %-% 1.8/12;% 1.9/12;
par.q_lo           = 1.40; %-% 1.2 ; % 1.5 - 0.1; % rt 2.5-0.1
par.q_lp           = 0.85; %-% 0.9; 23oct change from 0.9 to 0.8

par.rho_cell       = 214;

%_% add
par.I_en = 1000;%_% tuned 3000;  
par.k_en = 1.2e07; %_% tuned 3e07;  1.5e07;

%
par.a_n           = 1.99e-04; %_%
par.a_c           = 1.49e-04*(842.3/31.9);%_% 




end 

