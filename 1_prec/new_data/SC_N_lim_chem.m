function [D, cell, glucose, Jgy, Jeh, nh4, gas, ... 
precursor, atp, aa,  prot_sec, ...
lipid, protein,carbo, RNA, paper_info, color, shape, Jnh4] = SC_N_lim_chem()

plt_clrs = plot_colors;
data_unit_conv; 

color.empty = [0,0,0];
shape.empty = 'o';
% colors 
% high glucose, green with different shape
color.kumar_hg   = [102, 194, 165]./255;
shape.kumar_hg   = 'o'; 

color.kumar_ln   = [102, 194, 165]./255;
shape.kumar_ln   = 'o'; 

color.Yu2021     = [102, 194, 165]./255;
% low glucose, blue with different shape
color.hackett    = plt_clrs.blue; 
color.kumar_lg   = color.hackett;
color.boer       = color.hackett;


shape.Yu2021     = '^'; 

shape.hackett    = 'o'; 
shape.kumar_lg   = 'd'; 
shape.boer       = 's';

shape.lange     = 'v'; 
color.lange     = color.hackett;
shape.ertugay   = 'v'; 
color.ertugay   = color.hackett;
shape.Yu115     = 'v'; 
color.Yu115     = color.hackett;
shape.Yu50      = 'v'; 
color.Yu50      = color.hackett;
shape.Yu30      = 'v'; 
color.Yu30      = color.hackett;
%% Boer 2010 - glucose feed = 21 g/L, (NH4)2SO4 = 0.05 g/L
% strain: FY4
paper_info.CN.boer     = 2.2* 21/0.05; %924
paper_info.gl_in.boer  = 21;
paper_info.n_in.boer   = 0.05;
paper_info.author.boer = 'boer2010'; 
paper_info.strain.boer = 'FY4';

% unit h-1 
D.boer = [0.058;
0.100;
0.172;
0.208;
0.295;];
%_% [0.058;
% 0.108;
% 0.161;
% 0.238;
% 0.311;];

% original units: *10^7 cells/mL   
% converted units: cells/L
cell.boer = [3.47;
2.84;
2.31;
1.71;
0.82;]; 

% converted units: cells/L
% cell.boer = cell.boer .* 10^7 * 1000;   
% converted units: gDW/L
cell.boer = cell.boer .* 10^7 * 1000 * gdw_per_cell;    

% orignal unit mM
% converted unit g/L
glucose.boer = [
94.482;
102.761;
108.545;
107.838;
114.234;];
glucose.boer = glucose.boer*180/1000; 


%% Hackett 2016 - glucose feed = 0.8 g/L, (NH4)2SO4 = 5 g/L
% strain: FY4
paper_info.CN.hackett     = 2.2* 21/0.05; %924
paper_info.gl_in.hackett  = 21;
paper_info.n_in.hackett   = 0.05;
paper_info.author.hackett = 'hackett2016'; 
paper_info.strain.hackett = 'FY4';

% original units: h^-1
D.hackett = [0.052297056;
0.111120512;
0.15299851;
0.221107053;
0.252175169;];
%_% [0.049863024;
%_%   0.105431457;
%_%   0.154377453;
%_%   0.212650311;
%_%   0.29384141;
%_%   ];

% original units: moles / hr / mL cells
Jgy.hackett = [0.000860796;
0.00191866;
0.002878938;
0.003691687;
0.004394276];
Jgy.hackett = Jgy.hackett *10^9;

% ethanol production
Jeh.hackett = [0.001369106;
0.003099373;
0.004805673;
0.00562562;
0.006218425;]; 
% Nlim = 0.001369106	0.003099373	0.004805673	0.00562562	0.006218425                    
Jeh.hackett = Jeh.hackett *10^9;

%  original units: volume fraction (mL cells / L culture)
%  converted unit : number of cells
cell.hackett = [0.971834666
0.744836787;
0.656165741;
0.375965236;
0.375965236;];

% converted unit cell number            
%cell.hackett =  cell.hackett *10^-3./par.V_c;  
%_%
cell.hackett  = cell.hackett.* [0.216;0.2;0.195;0.239;0.182;];

% original units: M
% converted units: uM
atp.hackett = [4.38E-04;
5.14E-04;
4.95E-04;
5.00E-04;
4.87E-04];
atp.hackett = atp.hackett * 10^6;           

% original units: M 
aa.glutamine.hackett = [4.09E-04;
3.73E-04;
7.51E-04;
2.67E-03;
3.84E-03];
aa.glutamine.hackett = aa.glutamine.hackett .* 10^6;                  



% original units: M                 
precursor.pyruvate.hackett = [6.18E-03;
7.57E-03;
7.40E-03;
8.06E-03;
7.17E-03;];                 
precursor.pyruvate.hackett = precursor.pyruvate.hackett .* 10^6;   

precursor.g6p.hackett = [5.70E-04;
1.26E-03;
9.27E-04;
7.66E-04;
6.48E-04;];

precursor.g6p.hackett = precursor.g6p.hackett .* 10^6; 


precursor.f6p.hackett = [4.4500;
9.8700;
7.2400;
5.9900;
5.0700;];
%_% [1.00E+00;
% 2.22E+00;
% 1.63E+00;
% 1.34E+00;
% 1.14E+00]; 

% precursor 2

precursor.f16bp.hackett = [6.85E-04;
2.02E-03;
1.49E-03;
1.70E-03;
1.47E-03;]; 
precursor.f16bp.hackett = precursor.f16bp.hackett .* 10^6; 

precursor.g3p.hackett = [3.6000;
4.7200;
5.0200;
6.4000;
5.7500;];
%_% [1.00E+00;
% 1.31E+00;
% 1.39E+00;
% 1.78E+00;
% 1.60E+00]; 

precursor.dhap.hackett = [3.68E-04;
5.23E-04;
5.49E-04;
7.76E-04;
7.62E-04]; 
precursor.dhap.hackett = precursor.dhap.hackett .* 10^6; 

% precursor 2

% precursor.up13bpg.hackett = [
%                        ]; 
% 
% precursor.up3pg.hackett = [
%                      ]; 

precursor.pep.hackett = [38.9000;
46.0000;
38.4000;
19.7000;
16.9000;];
%_% [1.00E+00;
% 1.18E+00;
% 9.86E-01;
% 5.05E-01;
% 4.34E-01]; 



% unit: % 
% total carbo
carbo.all.hackett = [39.3451747;
36.76202613;
33.7677008;
30.45491392;
27.4622757;];
%_% [43.66484191;
% 40.49655727;
% 38.64183387;
% 36.2278625;
% 35.50744616;];

% unit % 
% glycogen + trehalose
carbo.sc.hackett = [28.124785;
17.2885038;
14.06576688;
8.036866206;
7.275767035;];
%_% [31.21257689;
% 19.04478501;
% 16.09606263;
% 9.560312156;
% 9.407228633;];

% unit: % 
protein.hackett = [30.60313317;
34.91759152;
35.41456124;
36.90563627;
33.71787362;];
%_% [33.96302042;
% 38.46475273;
% 40.52640718;
% 43.90136579;
% 43.59564355;];

% unit: % 
RNA.hackett = [3.246989934;
4.421566898;
4.956934081;
5.894308432;
5.664835403;];
%_% [3.603473697;
% 4.870739075;
% 5.672433086;
% 7.011617104;
% 7.324368902;];


% unit: % 
lipid.hackett = [7.189618495;
5.984086347;
4.150111405;
3.402777926;
2.715152311;];
%_% [7.978959488;
% 6.591989642;
% 4.749151161;
% 4.047799023;
% 3.510565751;];

% protein fraction %_%
% prot_sec.r.hackett  = [0.18717;    0.19328;    0.19967;    0.23729;     0.22422];
% prot_sec.z.hackett  = [0.35029;    0.32722;    0.32107;    0.29906;     0.28395];
% prot_sec.gy.hackett = [0.1657;     0.15562;    0.15582;    0.16112;     0.17831];
% prot_sec.as.hackett = [0.09968;    0.12869;    0.12954;    0.14207;     0.14935];
% prot_sec.fe.hackett = [0.12243;    0.1257;     0.12583;    0.11827;     0.11938];
% prot_sec.mt.hackett = [0.06624;    0.064527;   0.0637;     0.050524;    0.047319];
% prot_sec.at.hackett = [0.00073045; 0.00088671; 0.00074145; 0.00058486;  0.00051816];
% prot_sec.sp.hackett = [0.0034451;  0.0024221;  0.0022478;  0.0014556;   0.0013416];
% prot_sec.sd.hackett = [0.00012385; 8.5012e-05; 0.00020097; 5.9956e-05;  9.5516e-05];
% prot_sec.gn.hackett = [0.0091129;  0.009945;   0.01018;    0.008285;    0.0077487];
% prot_sec.lp.hackett = [0.051729;   0.064619;   0.065823;   0.056076;    0.05037]; % lp
% prot_sec.lo.hackett = [0.0018138;  0.0021585;  0.0023665;  0.0020928;   0.0020958]; % ld
% prot_sec.e.hackett  = [0.45466;    0.47066;    0.47203;    0.46802;     0.48884]; 
prot_sec.r.hackett  = [0.21527       0.21716       0.22487        0.2634       0.24667];
prot_sec.z.hackett  = [0.33171       0.32341       0.31622       0.29467       0.28887];
prot_sec.gy.hackett = [0.15324       0.13899       0.13889       0.14112       0.15536];
prot_sec.as.hackett = [0.12025       0.14877       0.15098       0.16741       0.17256];
prot_sec.fe.hackett = [0.11134       0.11009       0.11011        0.1015       0.10172];
prot_sec.mt.hackett = [0.067695      0.064587      0.063086      0.049025      0.045956];
prot_sec.at.hackett = [0.00063796    0.00075519    0.00061788    0.00051039    0.00045365];
prot_sec.sp.hackett = [0.003121     0.0020849     0.0019313     0.0012297      0.001126];
prot_sec.sd.hackett = [0.00012911    8.3874e-05    0.00019761    5.7819e-05    9.1431e-05];
prot_sec.gn.hackett = [0.0095013      0.010038      0.010173     0.0082629     0.0077645];
prot_sec.lp.hackett = [0.079012      0.087818      0.086584      0.071795      0.070252];
prot_sec.lo.hackett = [0.001704     0.0019448     0.0020569     0.0018205     0.0017518];
prot_sec.e.hackett  = [0.45351       0.45944        0.4609       0.45596       0.47104];


%% Kumar 2019 - low nitrogen condition 

% % n - limited 
% 0.2g/L  (NH4)2SO4
paper_info.CN.kumar_ln     = 2.2* 10/0.2;%_% 110
paper_info.gl_in.kumar_ln  = 10;
paper_info.n_in.kumar_ln   = 0.2;
paper_info.author.kumar_ln = 'kumar_ln'; 
paper_info.strain.kumar_ln = 'CEN.PK113-7D';

% original units: 1/h
D.kumar_ln = [0.06;
0.12;
0.19;
0.24;
0.32;
0.34;];

% original units: dcw (g/L)         
% converted units: cells/L_culture
cell.kumar_ln = [1.34;
1.2;
0.81;
0.62;
0.52;
0.46];    

% converted to cell numeber                  
% cell.kumar_ln = cell.kumar_ln * gdw_to_cell; 

% qs                   
% original units: g/gdw/h   
% converted units: uM glucose/h
Jgy.kumar_ln = [0.37;
0.69;
1.26;
1.47;
2.21;
3.12];     

Jgy.kumar_ln = Jgy.kumar_ln * one_over_gdw_to_1_over_cell * (1/par.V_c) * gl_gPerL_to_uM; 

% qEtOH
% original units: g/gdw/h    
% converted units: uM ethanol/h
Jeh.kumar_ln = [0.12;
0.23;
0.39;
0.56;
1.03;
1.39;];   

Jeh.kumar_ln = Jeh.kumar_ln * one_over_gdw_to_1_over_cell * (1/par.V_c) * eh_gPerL_to_uM;                            

% original units: umol/gdw       
% converted units: uM 
precursor.g6p.kumar_ln = [0.953058644;
1.079742226;
1.262891325;
1.216217794;
1.274408041;
1.331936948];                      

precursor.g6p.kumar_ln= precursor.g6p.kumar_ln * one_over_gdw_to_1_over_cell * (1/par.V_c);   

precursor.f16bp.kumar_ln = [3.10243265;
8.52690506;
18.02621392;
21.09185503;
13.74994041;
18.84488469];
precursor.f16bp.kumar_ln = precursor.f16bp.kumar_ln * one_over_gdw_to_1_over_cell * (1/par.V_c); 


precursor.dhap.kumar_ln =  [3.014063415;
8.297132164;
18.08986026;
21.31823945;
14.05490819;
19.58391102];
precursor.dhap.kumar_ln  = precursor.dhap.kumar_ln  * one_over_gdw_to_1_over_cell * (1/par.V_c); 



% original units: umol/gdw  
% converted units: uM 
precursor.pyruvate.kumar_ln = [1.624063185;
2.033738755;
4.932243575;
5.770782889;
5.852754237;
8.809834234;];           

precursor.pyruvate.kumar_ln = precursor.pyruvate.kumar_ln * one_over_gdw_to_1_over_cell * (1/par.V_c);                      

%kumar_hg.pyruvate_batch_gl = 4.7495203942839;
%kumar_hg.pyruvate_batch_gl = kumar_hg.pyruvate_batch_gl * one_over_gdw_to_1_over_cell * (1/par.V_c);                      

% original units: umol/gdw      
% converted units: uM 
atp.kumar_ln = [1.83409495;
1.772683276;
3.075001425;
3.791664609;
2.973749411;
4.437429781];  

atp.kumar_ln = atp.kumar_ln * one_over_gdw_to_1_over_cell * (1/par.V_c);                 

%kumar_hg.atp_batch_gl = 3.060351241;
%kumar_hg.atp_batch_gl = kumar_hg.atp_batch_gl *  one_over_gdw_to_1_over_cell * (1/par.V_c);   

% original units: umol/gdw     
% converted units: uM 
aa.glutamate.kumar_ln = [33.80304777;
40.43235895;
57.12187163;
73.38512241;
63.84082706;
92.54973683;];     

aa.glutamate.kumar_ln = aa.glutamate.kumar_ln * one_over_gdw_to_1_over_cell * (1/par.V_c);                       

% original units: umol/gdw      
% converted units: uM 
aa.glutamine.kumar_ln = [3.631471233;
6.714386533;
8.092000892;
28.89917549;
46.95101616;
90.60211341];      

aa.glutamine.kumar_ln = aa.glutamine.kumar_ln * one_over_gdw_to_1_over_cell * (1/par.V_c);      

%_% % original units: umol/gdw      
% % converted units: uM 
% aa.glutamine.kumar_ln = [33.80304777;
%                         40.43235895;
%                         57.12187163;
%                         73.38512241;
%                         63.84082706;
%                         92.54973683;];     
% 
% aa.glutamine.kumar_ln = aa.glutamine.kumar_ln * one_over_gdw_to_1_over_cell * (1/par.V_c);   

% original units: mmol/gdw/h      
% converted units: uM/h 
gas.co2.kumar_ln = [4.77;
8.42;
15.32;
18.06;
21.97;
27.63;]; 

gas.co2.kumar_ln = gas.co2.kumar_ln * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);   


% original units: mmol/gdw/h      
% converted units: uM/h 
gas.o2.kumar_ln = [3.23;
4.14;
5.84;
6.76;
8.04;
9.52]; 

gas.o2.kumar_ln = gas.o2.kumar_ln * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c); 

%% Yu2020
% strain: CEN.PK113-7D
% ref: Nitrogen limitation reveals large reserves inmetabolic and translational capacities of yeast
% C/N = 4  ; glucsoe: 7.5g/L (NH4)2SO4 5g/L
% C/N = 33 ; glucsoe: 7.5g/L (NH4)2SO4 0.5g/L
% C/N = 57 ; glucsoe: 7.5g/L (NH4)2SO4 0.29g/L
% C/N = 132; glucsoe: 7.5g/L (NH4)2SO4 0.125g/L

paper_info.gl_in.Yu30    = 7.5;
paper_info.n_in.Yu30     = 0.5;
paper_info.CN.Yu30       = 2.2* 7.5/0.5;
paper_info.author.Yu30   = 'Yu2020_CN30'; 
paper_info.strain.Yu30   = 'CEN.PK113-7D';

D.Yu30 = 0.2;
% unit: g/L
glucose.Yu30 = [0.23];

% unit: mmol/gDW/h
% converted units: uM/h 
Jgy.Yu30 = [4.62];
Jgy.Yu30  = Jgy.Yu30  * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);

% unit: mmol/gDW/h
% converted units: uM/h 
Jeh.Yu30 = [3.52];%_%[3.53];
Jeh.Yu30  = Jeh.Yu30  * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);

% unit:g/L
%
nh4.Yu30 = [0]; 
%nh4.Yu4  = 10^6 * nh4.Yu4 /18.04;


% unit: % 
protein.Yu30 = [36.65];%_%[36.86]; 

prot_sec.r.Yu30       = [0.28288];  %_%[0.38861      ];
prot_sec.z.Yu30       = [0.2928];   %_%[0.26005      ];
prot_sec.e.Yu30       = [0.42432];  %_%[0.35135      ];
prot_sec.gy.Yu30      = [0.11175];  %_%[0.095034     ];
prot_sec.fe.Yu30      = [0.019179]; %_%[0.012743    ];
prot_sec.gn.Yu30      = [0.010656]; %_%[0.006519     ];
prot_sec.mt.Yu30      = [0.091755]; %_%[0.11015      ];
prot_sec.sp.Yu30      = [0.001558]; %_%[0.0008326    ];
prot_sec.sd.Yu30      = [0.00027923];%_%[0.0001194    ];
prot_sec.as.Yu30      = [0.20403];  %_%[0.1362       ];
prot_sec.at.Yu30      = [0.0021456];%_%[0.0012325    ];
prot_sec.lp.Yu30      = [0.076366]; %_%[0.046869     ];
prot_sec.lo.Yu30      = [0.0032625];%_%[0.0022242    ];
%% Yu2020
paper_info.gl_in.Yu50    = 7.5;
paper_info.n_in.Yu50     = 0.29;
paper_info.CN.Yu50       = 2.2* 7.5/0.29;
paper_info.author.Yu50   = 'Yu2020_CN50'; 
paper_info.strain.Yu50   = 'CEN.PK113-7D';


D.Yu50 = 0.2;
% unit: g/L
glucose.Yu50 = [1.19 ];

% unit: mmol/gDW/h
% converted units: uM/h 
Jgy.Yu50 = [5.80];%_%[5.56 ];
Jgy.Yu50  = Jgy.Yu50  * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);

% unit: mmol/gDW/h
% converted units: uM/h 
Jeh.Yu50 = [4.19];%_%[4.07 ];
Jeh.Yu50  = Jeh.Yu50  * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);

% unit:g/L
%
nh4.Yu50 = [0]; 
%nh4.Yu4  = 10^6 * nh4.Yu4 /18.04;

protein.Yu50 = [26.35];%_%[23.79]; 

prot_sec.r.Yu50       = [0.2689];   %_%[0.36777     ];
prot_sec.z.Yu50       = [0.30736];  %_%[0.27355     ];
prot_sec.e.Yu50       = [0.42374];  %_%[0.35868     ];
prot_sec.gy.Yu50      = [0.13205];  %_%[0.11316     ];
prot_sec.fe.Yu50      = [0.021253]; %_%[0.013936    ];
prot_sec.gn.Yu50      = [0.0090348];%_%[0.0055124   ];
prot_sec.mt.Yu50      = [0.089557]; %_%[0.11276     ];
prot_sec.sp.Yu50      = [0.0017143];%_%[0.00092687  ];
prot_sec.sd.Yu50      = [0.00046121];%_%[0.00019775  ];
prot_sec.as.Yu50      = [0.18316];  %_%[0.1214      ];
prot_sec.at.Yu50      = [0.0023849];%_%[0.0013767   ];
prot_sec.lp.Yu50      = [0.078098]; %_%[0.048185    ];
prot_sec.lo.Yu50      = [0.0034967];%_%[0.002361    ];

%% Yu2020
paper_info.gl_in.Yu115    = 7.5;
paper_info.n_in.Yu115     = 0.125;
paper_info.CN.Yu115       = 2.2* 7.5/0.125;
paper_info.author.Yu115   = 'Yu2020_CN115'; 
paper_info.strain.Yu115   = 'CEN.PK113-7D';

D.Yu115 = 0.2;
% unit: g/L
glucose.Yu115 = [4.01];

% unit: mmol/gDW/h
% converted units: uM/h 
Jgy.Yu115 = [5.98];%_%[5.97];
Jgy.Yu115  = Jgy.Yu115  * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);

% unit: mmol/gDW/h
% converted units: uM/h 
Jeh.Yu115 = [4.86];%_%[4.83];
Jeh.Yu115  = Jeh.Yu115  * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);

% unit:g/L
%
nh4.Yu115 = [0]; 
%nh4.Yu4  = 10^6 * nh4.Yu4 /18.04;

protein.Yu115 = [26.35];%_%[26.41]; 

prot_sec.r.Yu115       = [0.28288]; %_%[0.39347];
prot_sec.z.Yu115       = [0.30012]; %_%[0.26539];
prot_sec.e.Yu115       = [0.41701]; %_%[0.34114];
prot_sec.gy.Yu115      = [0.13438]; %_%[0.11152];
prot_sec.fe.Yu115      = [0.019359];%_%[0.012489];
prot_sec.gn.Yu115      = [0.0093413];%_%[0.0056282];
prot_sec.mt.Yu115      = [0.080282];%_%[0.099003];
prot_sec.sp.Yu115      = [0.0014993];%_%[0.00078952];
prot_sec.sd.Yu115      = [0.00039688];%_%[0.00016743];
prot_sec.as.Yu115      = [0.18458]; %_%[0.12031];
prot_sec.at.Yu115      = [0.0025211];%_%[0.0014287];
prot_sec.lp.Yu115      = [0.076968];%_%[0.044967];
prot_sec.lo.Yu115      = [0.003308];%_%[0.0022055];

%% Yu2021
paper_info.gl_in.Yu2021    = 7.5;
paper_info.n_in.Yu2021     = 0.5;
paper_info.CN.Yu2021       = 2.2* 7.5/0.5; %_%33
paper_info.author.Yu2021   = 'Yu2021'; 
paper_info.strain.Yu2021   = 'CEN.PK113-7D';


% unit: % 
D.Yu2021 = [0.05;
0.1;
0.13;
0.18;
0.3;
0.35]; 

% unit: % 
RNA.Yu2021 = [0.024168667;
0.025297333;
0.025831333;
0.038539;
0.048218667;
0.077494;]*100;


% unit: % 
protein.Yu2021 = [0.261477333;
0.225573;
0.275686667;
0.286470333;
0.414481333;
0.494313667;]*100;

% unit: (mmol/gDW h)
% converted units: uM glucose/h
Jgy.Yu2021 = [0.765;
1.503;
2.071;
3.099;
8.400;
13.088;];

Jgy.Yu2021  = Jgy.Yu2021  * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);

% original units: mmol/gdw/h      
% converted units: uM/h 
Jnh4.Yu2021 = [0.119;
0.307;
0.426;
0.636;
1.736;
2.747;];

Jnh4.Yu2021 = Jnh4.Yu2021  * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);                                                      

% original units: mmol/gdw/h      
% converted units: uM/h 
Jeh.Yu2021 = [0.000;
0.000;
0.139;
1.279;
7.314;
11.967;];

Jeh.Yu2021 = Jeh.Yu2021  * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);                                                      

% original units: mmol/gdw/h      
% converted units: uM/h 
gas.o2.Yu2021 = [2.088;
3.654;
4.956;
5.567;
8.676;
12.671;];

gas.o2.Yu2021 = gas.o2.Yu2021  * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);                                                      


% original units: mmol/gdw/h      
% converted units: uM/h 
gas.co2.Yu2021 = [2.311;
3.254;
5.395;
6.939;
15.384;
21.028;];

gas.co2.Yu2021 = gas.co2.Yu2021  * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);                                                      

%_%
% prot_sec.r.Yu2021       = [0.3243;       0.34773;       0.37614;       0.41126;       0.48692;       0.53098];
% prot_sec.z.Yu2021       = [0.2562;       0.23447;       0.22414;       0.2102;        0.17122;       0.15636];
% prot_sec.e.Yu2021       = [0.4195;       0.41779;       0.39972;       0.37855;       0.34186;       0.31265];
% prot_sec.gy.Yu2021      = [0.16462;      0.15809;       0.13533;       0.11405;       0.11388;       0.11901];
% prot_sec.fe.Yu2021      = [0.013359;     0.014047;      0.013182;      0.013679;      0.011977;      0.014359];
% prot_sec.gn.Yu2021      = [0.029678;     0.013467;      0.010774;      0.010702;      0.0082749;     0.0068022];
% prot_sec.mt.Yu2021      = [0.10985;      0.12245;       0.13225;       0.11918;       0.094652;      0.076595];
% prot_sec.sp.Yu2021      = [0.0010466;    0.00096197;    0.0010796;     0.00096255;    0.00055291;    0.00046649];
% prot_sec.sd.Yu2021      = [0.00012885;   0.00014574;    0.00012417;    0.00012057;    7.0591e-05;    5.3142e-05];
% prot_sec.as.Yu2021      = [0.10417;      0.11422;       0.11366;       0.12627;       0.11641;       0.098446];
% prot_sec.at.Yu2021      = [0.010038;     0.011354;      0.0096738;     0.0093929;     0.0054994;     0.00414];
% %prot_sec.lp.Yu2021      = [0.044515;     0.042935;      0.043953;       0.04345;      0.041485;      0.03388];
% prot_sec.lp.Yu2021      = [0.043675;      0.042127;      0.043193;      0.042683;      0.040967;      0.033517];
% prot_sec.lo.Yu2021      = [0.0022021;     0.0023802;     0.0030848;     0.0027818;     0.0019385;     0.0014351];
prot_sec.r.Yu2021    = [0.2464      0.26419       0.28117       0.29262       0.34435       0.38054];
prot_sec.z.Yu2021    = [0.27452      0.25452       0.25263       0.24341       0.21113       0.19808];
prot_sec.e.Yu2021    = [0.47908      0.48128        0.4662       0.46397       0.44451       0.42139];
prot_sec.gy.Yu2021   = [0.19232      0.18995       0.16635       0.14482        0.1493       0.15801];
prot_sec.fe.Yu2021   = [0.020251     0.021991      0.021232      0.022555      0.020846      0.026057];
prot_sec.gn.Yu2021   = [0.044384     0.021103      0.017337      0.017868      0.014621      0.012061];
prot_sec.mt.Yu2021   = [0.080699     0.091852       0.10192      0.095899      0.081704       0.07121];
prot_sec.sp.Yu2021   = [0.0018633     0.001761     0.0019878     0.0018348     0.0010833    0.00091211];
prot_sec.sd.Yu2021   = [0.00029593    0.0003417    0.00029774    0.00029724    0.00017994    0.00013716];
prot_sec.as.Yu2021   = [0.14727       0.1665       0.17135       0.19489        0.1866       0.16087];
prot_sec.at.Yu2021   = [0.017652     0.020382       0.01776       0.01773      0.010733     0.0081812];
prot_sec.lp.Yu2021   = [0.06239      0.06491      0.068329      0.069286      0.066392      0.054606];
prot_sec.lo.Yu2021   = [0.0031895    0.0034194     0.0046112     0.0042383     0.0030654     0.0023134];
end 

