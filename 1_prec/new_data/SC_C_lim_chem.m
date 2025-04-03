function [D, cell, glucose, Jgy, Jeh,... 
precursor, atp, aa,  prot_sec, ...
lipid, protein,carbo, RNA, paper_info, color, shape, ...
gas] = SC_C_lim_chem() 


plt_clrs = plot_colors;
data_unit_conv; 
color.empty = [0,0,0];
shape.empty = 'o';
% % colors 
% % high glucose, green with different shape
color.kumar_hg   = [102, 194, 165]./255;
color.xia        = color.kumar_hg;
% low glucose, blue with different shape
color.hackett    = plt_clrs.blue; 
color.kumar_lg   = color.hackett;
color.boer       = color.hackett;

shape.kumar_hg   = 'o'; 
shape.xia        = '^'; 

shape.hackett    = 'o'; 
shape.kumar_lg   = 'd'; 
shape.boer       = 's';

shape.lange     = 'v'; 
color.lange     = color.hackett;
shape.ertugay   = 'v'; 
color.ertugay   = color.hackett;
shape.Yu4       = 'v'; 
color.Yu4       = color.hackett;

shape.Suarez_Mendez       = 'v'; 
color.Suarez_Mendez       = color.kumar_hg;

%% Kumar 2019 - low glucose condition 
% 1 g L−1 glucose
paper_info.CN.kumar_lg     = 2.2* 1/5; %_% 0.44;
paper_info.gl_in.kumar_lg  = 1;
paper_info.n_in.kumar_lg   = 5;
paper_info.author.kumar_lg = 'kumar_lg'; 
paper_info.strain.kumar_lg = 'CEN.PK113-7D';

% original units: 1/h
D.kumar_lg = [0.05;
0.12;
0.18;
0.24;
0.31];

% original units: dcw (g/L)       
% converted units: cells/L_culture
cell.kumar_lg = [0.32;
0.44;
0.34;
0.33;
0.30];  

%cell.kumar_lg = cell.kumar_lg * gdw_to_cell;                   

% qs                   
% original units: g/gdw/h         
% converted units: uM glucose/h
Jgy.kumar_lg = [0.16;
0.27;
0.53;
0.73;
1.03];     

Jgy.kumar_lg = Jgy.kumar_lg * one_over_gdw_to_1_over_cell * (1/par.V_c) * gl_gPerL_to_uM;                 

% qEtOH
% original units: g/gdw/h     
% converted units: uM ethanol/h
Jeh.kumar_lg = [0;
0;
0;
0;
0];   

Jeh.kumar_lg = Jeh.kumar_lg * one_over_gdw_to_1_over_cell * (1/par.V_c) * eh_gPerL_to_uM;                                                      

% original units: umol/gdw       
% converted units: uM 
precursor.g6p.kumar_lg = [1.0645163061355;
1.8072881736487;
2.14282912259462;
2.98231109116389;
1.63578591111586];                      

precursor.g6p.kumar_lg =  precursor.g6p.kumar_lg * one_over_gdw_to_1_over_cell * (1/par.V_c);    

% original units: umol/gdw       
% converted units: uM 
precursor.pyruvate.kumar_lg = [3.325;
4.871;
8.989;
8.326;
5.893];        

precursor.pyruvate.kumar_lg = precursor.pyruvate.kumar_lg * one_over_gdw_to_1_over_cell * (1/par.V_c);                                           

% original units: umol/gdw        
% converted units: uM 
atp.kumar_lg = [5.347;
7.041;
8.852;
9.231;
6.684];  

atp.kumar_lg = atp.kumar_lg * one_over_gdw_to_1_over_cell * (1/par.V_c);                                 

% original units: umol/gdw         
% converted units: uM 
aa.glutamate.kumar_lg = [236.8;
366.4;
413.4;
426.0;
324.7];     

aa.glutamate.kumar_lg = aa.glutamate.kumar_lg * one_over_gdw_to_1_over_cell * (1/par.V_c);                                             

% original units: umol/gdw      
% converted units: uM 
aa.glutamine.kumar_lg = [43.51;
70.30;
88.44;
111.04; 
108.32;];
%_% [40.51;
%   70.30;
%   88.44;
%   111.0;
%   108.3];         

aa.glutamine.kumar_lg = aa.glutamine.kumar_lg * one_over_gdw_to_1_over_cell * (1/par.V_c);                       

% original units: mmol/gdw/h      
% converted units: uM/h 
gas.co2.kumar_lg = [2.28;
3.95;
7.52;
8.30;
10.10]; 

gas.co2.kumar_lg= gas.co2.kumar_lg * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);  


% original units: mmol/gdw/h      
% converted units: uM/h 
gas.o2.kumar_lg = [2.42;
3.47;
7.51;
8.11;
10.45;]; 

gas.o2.kumar_lg = gas.o2.kumar_lg * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c); 


%% Kumar 2019 - high glucose condition 
% 10 g L−1 glucose
paper_info.CN.kumar_hg     = 2.2* 10/5;
paper_info.gl_in.kumar_hg  = 10;
paper_info.n_in.kumar_hg   = 5;
paper_info.author.kumar_hg = 'kumar_hg'; 
paper_info.strain.kumar_hg = 'CEN.PK113-7D';


% original units: 1/h
D.kumar_hg = [0.12;
0.26;
0.36;
0.41];

% original units: dcw (g/L)         
% converted units: cells/L_culture
cell.kumar_hg = [4.25;
4.40;
1.71;
0.55];    

%cell.kumar_hg = cell.kumar_hg * gdw_to_cell; 

% qs                   
% original units: g/gdw/h   
% converted units: uM glucose/h
Jgy.kumar_hg = [0.28;
0.59;
2.05;
3.76];     

Jgy.kumar_hg = Jgy.kumar_hg * one_over_gdw_to_1_over_cell * (1/par.V_c) * gl_gPerL_to_uM; 

% qEtOH
% original units: g/gdw/h    
% converted units: uM ethanol/h
Jeh.kumar_hg = [0;
0;
0.57;
1.31];   

Jeh.kumar_hg = Jeh.kumar_hg * one_over_gdw_to_1_over_cell * (1/par.V_c) * eh_gPerL_to_uM;                            

% original units: umol/gdw       
% converted units: uM 
precursor.g6p.kumar_hg = [1.59750985308395;
1.65884964663491;
0.606088883651831;
1.65780339285714];                      

precursor.g6p.kumar_hg= precursor.g6p.kumar_hg * one_over_gdw_to_1_over_cell * (1/par.V_c);   

precursor.f16bp.kumar_hg = [1.065358019;
1.861046866;
3.788328063;
38.97298214;
];
precursor.f16bp.kumar_hg = precursor.f16bp.kumar_hg * one_over_gdw_to_1_over_cell * (1/par.V_c); 


precursor.dhap.kumar_hg =  [0.944228946;
1.666493961;
3.71221835;
43.71398438;
];
precursor.dhap.kumar_hg  = precursor.dhap.kumar_hg  * one_over_gdw_to_1_over_cell * (1/par.V_c); 



% original units: umol/gdw  
% converted units: uM 
precursor.pyruvate.kumar_hg = [3.882;
4.427;
3.68;
7.420];           

precursor.pyruvate.kumar_hg = precursor.pyruvate.kumar_hg * one_over_gdw_to_1_over_cell * (1/par.V_c);                      

%kumar_hg.pyruvate_batch_gl = 4.7495203942839;
%kumar_hg.pyruvate_batch_gl = kumar_hg.pyruvate_batch_gl * one_over_gdw_to_1_over_cell * (1/par.V_c);                      

% original units: umol/gdw      
% converted units: uM 
atp.kumar_hg = [6.291;
6.970;
5.517;
6.292];  

atp.kumar_hg = atp.kumar_hg * one_over_gdw_to_1_over_cell * (1/par.V_c);                 

%kumar_hg.atp_batch_gl = 3.060351241;
%kumar_hg.atp_batch_gl = kumar_hg.atp_batch_gl *  one_over_gdw_to_1_over_cell * (1/par.V_c);   

% original units: umol/gdw     
% converted units: uM 
aa.glutamate.kumar_hg = [237.9;
245.7;
124.8;
106.5];     

aa.glutamate.kumar_hg = aa.glutamate.kumar_hg * one_over_gdw_to_1_over_cell * (1/par.V_c);                       

% original units: umol/gdw      
% converted units: uM 
aa.glutamine.kumar_hg = [47.53;
111.8;
195.5;
293.3];      

aa.glutamine.kumar_hg = aa.glutamine.kumar_hg * one_over_gdw_to_1_over_cell * (1/par.V_c);      

% original units: mmol/gdw/h      
% converted units: uM/h 
gas.co2.kumar_hg = [3.50;
5.84;
19.26;
26.68]; 

gas.co2.kumar_hg = gas.co2.kumar_hg * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);   


% original units: mmol/gdw/h      
% converted units: uM/h 
gas.o2.kumar_hg = [3.35;
5.6;
10.14;
9.37]; 

gas.o2.kumar_hg = gas.o2.kumar_hg * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c); 


% % n - limited 
% % 0.2g/L  (NH4)2SO4
% paper_info.CN.kumar_ln     = 2.2* 10/0.2;
% paper_info.gl_in.kumar_ln  = 10;
% paper_info.n_in.kumar_ln   = 0.2;
% paper_info.author.kumar_ln = 'kumar_hg'; 
% paper_info.strain.kumar_ln = 'CEN.PK113-7D';

%% Boer 2010 - glucose feed = 0.8 g/L, (NH4)2SO4 = 5 g/L
% strain: FY4
paper_info.CN.boer     = 2.2* 0.8/5;%_%0.352
paper_info.gl_in.boer  = 0.8;
paper_info.n_in.boer   = 5;
paper_info.author.boer = 'boer2010'; 
paper_info.strain.boer = 'FY4';

% unit h-1 
D.boer = [0.058;
0.108;
0.161;
0.238;
0.311;];

% original units: *10^7 cells/mL   

cell.boer = [3.52;
4.07;
4.69;
4.46;
1.13;];
%_% [3.507846608;
%  4.063205506;
%  4.680550639;
%  4.265526057;
%  1.102261554;
%  ]; 
% converted units: cells/L
% cell.boer = cell.boer .* 10^7 * 1000;   
% converted units: gDW/L
cell.boer = cell.boer .* 10^7 * 1000 * gdw_per_cell;  

% orignal unit mM
% converted unit g/L
glucose.boer = [
0.104;
0.118;
0.149;
0.252;
1.133;];
glucose.boer = glucose.boer*180/1000; 


% original units: M/h       
% converted units: uM/h
Jgy.boer = [0.246575342;
4.52E-01;
6.44E-01;
8.49E-01;
3.287671233;
];       
Jgy.boer = Jgy.boer .* 10^6; 


% original units: M/h        
Jeh.boer = [0.00E+00;
0.00E+00;
0.00E+00;
0.01369863;
3.712328767;
];   
Jeh.boer = Jeh.boer .* 10^6; 

% fold change relative to value at lowest dilution rate
precursor.pyruvate.boer = [1;
1.439018711;
1.760144149;
3.020963998;
4.187441002;]; 



precursor.f16bp.boer = [1;
1.058587296;
1.348334043;
1.126903374;
14.32443253;]; 


aa.glutamine.boer = [1;
1.439018711;
1.760144149;
3.020963998;
4.187441002;]; 

atp.boer = [1;
1.309703238;
1.163730571;
1.022060572;
1.1615593;];

%% Hackett 2016 - glucose feed = 0.8 g/L, (NH4)2SO4 = 5 g/L
% strain: FY4
paper_info.CN.hackett = 2.2* 0.8/5; %_%0.352
paper_info.gl_in.hackett = 0.8;
paper_info.n_in.hackett = 5;
paper_info.author.hackett = 'hackett2016'; 
paper_info.strain.hackett = 'FY4';

% original units: h^-1
D.hackett = [0.049863024;
0.105431457;
0.154377453;
0.212650311;
0.29384141;
];

% original units: moles / hr / mL cells
Jgy.hackett = [0.000149923;
0.00031825;
0.0004743;
0.001085886;
0.001886523;];
% Nlim = 0.000860796	0.00191866	0.002878938	0.003691687	0.004394276            
Jgy.hackett = Jgy.hackett *10^9;

% ethanol production
Jeh.hackett = [2.5945E-06;
1.09181E-05;
2.36939E-05;
0.000976161;
0.002583142;]; 
% Nlim = 0.001369106	0.003099373	0.004805673	0.00562562	0.006218425                    
Jeh.hackett = Jeh.hackett *10^9;


% original units: M
% converted units: uM
atp.hackett = [9.86E-04;
1.03E-03;
1.06E-03;
1.02E-03;
0.0008924;
];
atp.hackett = atp.hackett * 10^6;           

% original units: M 
aa.glutamine.hackett = [2.68E-03;
2.07E-03;
2.21E-03;
3.11E-03;
0.003401716;
];
aa.glutamine.hackett = aa.glutamine.hackett .* 10^6;                  

% original units: M                 
precursor.pyruvate.hackett = [1.82E-03;
2.01E-03;
2.18E-03;
2.48E-03;
0.003439365;
];                 
precursor.pyruvate.hackett = precursor.pyruvate.hackett .* 10^6;   

precursor.g6p.hackett = [0.0006946;
0.001079212;
0.001333398;
0.0010472;
0.001732866;
];

precursor.g6p.hackett = precursor.g6p.hackett .* 10^6; 

precursor.f6p.hackett = [5.4300;
8.4300;
10.4000;
8.1800;
13.5000;];
%_% [1
% 1.553716991
% 1.919663766
% 1.507630963
% 2.494769314
% ]; 

% precursor 2

precursor.f16bp.hackett =[6.78E-05;
1.43E-04;
2.13E-04;
2.24E-04;
1.06E-03]; 
precursor.f16bp.hackett = precursor.f16bp.hackett .* 10^6; 
%_% [1
%  2.109518429
%  3.139476988
%  3.297588489
%  15.61308059
%  ]; 

precursor.g3p.hackett = [1.0600;
1.4200;
1.6900;
2.1700;
3.6100;];
%_%  [1
% 1.33858714
% 1.590362163
% 2.039832976
% 3.393954087
% ]; 
% dhap = dihydroxyacetone phosphate = glycerone phosphate
% g6p -> dhap
precursor.dhap.hackett = [9.68E-05;
1.39E-04;
1.72E-04;
2.14E-04;
4.85E-04];
%_% [1
% 1.438194796
% 1.780174682
% 2.212717054
% 5.008825635
% ]; 
precursor.dhap.hackett = precursor.dhap.hackett .* 10^6; 

% precursor 2

%_% precursor.up13bpg.hackett = [1
%                        1.273250281
%                        1.507094881
%                        1.593551817
%                        2.959605966
%                        ]; 
% 
%_% precursor.up3pg.hackett = [1
%                      1.284537473
%                      1.363993295
%                      0.765997233
%                      0.697589241
%                      ]; 

precursor.pep.hackett = [46.0000;
59.0000;
62.7000;
35.2000;
32.1000;];
%_% [1
% 1.284537473
% 1.363993295
% 0.765997233
% 0.697589241
% ]; 



% original units: volume fraction (mL cells / L culture)
cell.hackett = [1.443564631;
1.432924105;
1.390362003;
0.819320466;
0.535573119;];
% converted unit cell number            
%cell.hackett =  cell.hackett *10^-3./par.V_c;  

% Cell density = Volume_Fraction_Mean/gDCW/mL = (mL cells/L culture) * (gDCW/mL cells) = (gDCW/L culture)
% (gDCW/mL) data is from hackett SI
cell.hackett  = cell.hackett.* [0.204;0.205;0.218;0.211;0.187;];

% unit: % 
carbo.all.hackett = [34.97146709;
30.02143594;
30.70827135;
30.36643108;
28.53761839;];
%_% [42.79427713;
% 35.93818917;
% 36.30855581;
% 36.19544632;
% 34.03065389;];

% unit: % 
carbo.sc.hackett = [15.0157121;
11.96277372;
10.26023963;
9.288902439;
6.918572674;];
%_% [18.37459502;
% 14.32044842;
% 12.13140521;
% 11.07196195;
% 8.250287355;];

% unit: % 
protein.hackett = [35.22267511;
40.37811588;
39.05909607;
39.00558413;
39.79534867;];
%_% [43.10167818;
% 48.33600798;
% 46.18232506;
% 46.49293565;
% 47.45531735;];

% unit: % 
RNA.hackett = [4.217466362;
5.533685826;
5.898711223;
6.486424768;
7.379780263;];
%_% [5.160876547;
% 6.624288341;
% 6.974462457;
% 7.731532192;
% 8.800270034;];


% unit: % 
lipid.hackett = [3.081364372;
2.714001091;
3.663140322;
3.269685921;
2.825419903;];
%_% [3.77063852;
% 3.248888055;
% 4.331189253;
% 3.897321384;
% 3.369268084;];

% protein fraction
prot_sec.r.hackett = [0.21676       0.23749       0.23757       0.24183       0.29104];
%_% [0.19336;
% 0.21676;
% 0.22049;
% 0.22944;
% 0.27];


prot_sec.z.hackett = [0.32434       0.32102        0.3291        0.3479       0.31226];
%_% [0.34458;
% 0.32434;
% 0.32861;
% 0.3408;
% 0.3119];


prot_sec.gy.hackett = [0.1456       0.12535       0.11252       0.10637        0.1003];
%_% [0.16035;
% 0.1456;
% 0.13033;
% 0.12546;
% 0.11892];


prot_sec.as.hackett = [0.13355       0.15542        0.1557       0.15662       0.17488];
%_% [0.1149;
% 0.136;
% 0.13735;
% 0.13708;
% 0.15092];


prot_sec.fe.hackett = [0.051799      0.045571      0.049288      0.049646      0.062169];
%_% [0.055987;
% 0.051799;
% 0.056046;
% 0.058497;
% 0.073864];


prot_sec.mt.hackett = [0.10154       0.10147       0.11114       0.11822      0.093798];
%_% [0.097286;
% 0.10154;
% 0.11225;
% 0.12148;
% 0.097007];


prot_sec.at.hackett = [0.00032893    0.00026669    0.00025858    0.00026592    0.00032497];
%_% [0.00039301;
% 0.00032893;
% 0.00031728;
% 0.0003515;
% 0.00041931];


prot_sec.sp.hackett = [0.0031839     0.0026769     0.0017296     0.0018074     0.0016262];
%_% [0.003533;
% 0.0031839;
% 0.0020732;
% 0.0022047;
% 0.0019915];


prot_sec.sd.hackett = [7.2971e-05    7.0027e-05    5.1645e-05    5.1825e-05    6.2975e-05];
%_% [7.0763e-05;
% 7.2971e-05;
% 5.4369e-05;
% 5.5465e-05;
% 6.7762e-05];


prot_sec.gn.hackett = [0.042793      0.039622      0.036324      0.016695      0.013844];
%_% [0.043389;
% 0.042793;
% 0.040008;
% 0.018624;
% 0.015086];


prot_sec.lp.hackett = [0.070659      0.089737      0.091385      0.086333      0.090112];
%_% [0.066146;
% 0.070659;
% 0.075232;
% 0.070144;
% 0.069507]; % 

prot_sec.lo.hackett = [0.002759     0.0026167     0.0033906     0.0024383     0.0022397 ];
%_% [0.0028258;
% 0.002759;
% 0.0037707;
% 0.0026818;
% 0.0025683]; % 

prot_sec.e.hackett = [0.45889       0.44917       0.44222       0.42359        0.4242];
%_% [0.45579;
% 0.45889;
% 0.45298;
% 0.43716;
% 0.43455]; 

%% Data: ﻿lange2000: 
% ref: ﻿Statistical Reconciliation of the Elemental and Molecular Biomass Composition of Saccharomyces cerevisiae
% glcuose: 10g/L (NH4)2SO4 5g/L
% strain: CEN.PK113-7D

paper_info.CN.lange     = 2.2* 10/5;
paper_info.gl_in.lange  = 10;
paper_info.n_in.lange   = 5;
paper_info.author.lange = 'lange2000'; 
paper_info.strain.lange = 'CEN.PK113-7D';

% unit: h-1

D.lange = [0.022 ;
0.052 ;
0.087 ;
0.107 ;
0.126 ;
0.158 ;
0.211 ;];


% unit: % 
protein.lange = [38.4;
42.0 ;
42.5 ;
41.4 ;
44.5 ;
44.0 ;
46.3;] ;


% unit: % 
carbo.all.lange = [ 45.4;
42.1 ;
41.2 ;
38.2 ;
38.5 ;
31.6 ;
31.9;];

% unit: % 
lipid.lange = [10.2;
9.4;
7.6;
7.2;
10.1;
7.8;
7.7;];


% unit: % 
RNA.lange = [4.3 ;
5.2 ;
6.2 ;
7.8 ;
6.8 ;
7.0 ;
7.9 ;];


%% Xia 2021 (Xia_protein) 
% strain: CEN.PK113-7D
% ref: Proteome allocations change linearly with 
% specific growth rate of Saccharomyces cerevisiae under glucose-limitation
% glucsoe: 10g/L (NH4)2SO4 5g/L

paper_info.CN.xia       = 2.2* 10/5;
paper_info.gl_in.xia    = 10;
paper_info.n_in.xia     = 5;
paper_info.author.xia   = 'xia2021'; 
paper_info.strain.xia   = 'CEN.PK113-7D';


D.xia =   [0.027; 
0.044 ;
0.102 ;
0.152 ;
0.214 ;
0.254 ;
0.284 ;
0.334 ;
0.379 ;];

% original units: mmol/gdw/h         
% converted units: uM glucose/h
Jgy.xia =  [0.399;
0.553;
1.14;
1.66;
2.40;
2.74;
3.20;
4.88;
6.78;];

Jgy.xia  = Jgy.xia  * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);


% original units: mmol/gdw/h         
% converted units: uM glucose/h
Jeh.xia =  [0;
0;
0;
0;
0;
0;
0;
1.77;
4.795;];

Jeh.xia  = Jeh.xia  * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);

% original units: mmol/gdw/h         
% converted units: uM glucose/h
gas.o2.xia  = [0.922;
1.533;
2.597;
3.981;
5.846;
5.978;
8.002;
9.578;
9.581;];
gas.o2.xia = gas.o2.xia * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);

% original units: mmol/gdw/h         
% converted units: uM glucose/h
gas.co2.xia = [0.858;
1.533;
2.661;
3.948;
5.91;
6.78;
8.515;
12.241;
14.426;];
gas.co2.xia = gas.co2.xia * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);

% unit: % (source: table S6)
protein.xia = [30; 
31.2;
37.1;
41.9;
32;
37.8;
39.1;
41.9;
39.1;];


% protein fraction
prot_sec.r.xia = [0.21675;
0.21687;
0.23151;
0.24334;
0.25812;
0.2697;
0.27317;
0.29511;
0.30961];


prot_sec.z.xia = [0.32434;
0.31028;
0.30487;
0.3014;
0.30091;
0.29396;
0.29096;
0.2804;
0.27598];


prot_sec.gy.xia = [0.1456;
0.13568;
0.096381;
0.083352;
0.06844;
0.064459;
0.067865;
0.068689;
0.069397];


prot_sec.as.xia = [0.136;
0.15812;
0.19625;
0.20131;
0.2041;
0.20552;
0.20906;
0.21236;
0.21075];


prot_sec.fe.xia = [0.051798;
0.049615;
0.045992;
0.045843;
0.042169;
0.039195;
0.041661;
0.042323;
0.050627];


prot_sec.mt.xia = [0.10154;
0.10653;
0.11158;
0.11633;
0.12278;
0.13013;
0.1357;
0.12997;
0.11184];


prot_sec.at.xia = [0.0034551;
0.0032158;
0.0027016;
0.0024044;
0.0031752;
0.0027205;
0.0025897;
0.0025053;
0.0027454];


prot_sec.sp.xia = [0.0031839;
0.0030257;
0.0025914;
0.0022106;
0.0017963;
0.0013729;
0.0011977;
0.0012013;
0.0011006];


prot_sec.sd.xia = [ 7.297e-05;
7.5645e-05;
6.099e-05;
3.4624e-05;
2.6262e-05;
1.6389e-05;
1.7031e-05;
1.801e-05;
1.5826e-05];


prot_sec.gn.xia = [0.042792;
0.041331;
0.036239;
0.033547;
0.033024;
0.031023;
0.019144;
0.010441;
0.0083861];


prot_sec.lp.xia = [0.070658;
0.087618;
0.11079;
0.10529;
0.092367;
0.077567;
0.072494;
0.063379;
0.062047]; % lipid


prot_sec.lo.xia = [0.002759;
0.0025753;
0.0025129;
0.0025095;
0.00244;
0.0021752;
0.0019079;
0.001586;
0.0014412]; 


prot_sec.e.xia =   [0.45888;
0.47168;
0.46408;
0.45711;
0.44489;
0.44132;
0.4411;
0.43196;
0.42357]; 

%% Ertugay1997
% ref: Continuous cultivation of bakers' yeast: change in cell composition at different dilution rates and effect of heat stress on trehalose level
% strain: uknown ﻿some commercial strain
% ﻿sucrose 3 g/L, (NH4)2SO4 2g/L, ﻿yeast extract (Difco) 2g/L,
paper_info.CN.ertugay     = 2.2* 3/4;
paper_info.gl_in.ertugay  = 3;
paper_info.n_in.ertugay   = 4;
paper_info.author.ertugay = 'ertugay1997'; 
paper_info.strain.ertugay = 'unknown';

D.ertugay =[0.10;
0.15;
0.20;
0.25;
0.30;
0.35;
0.40;]; 


protein.ertugay =[41.718;
42.477;
43.731;
44.591;
44.459;
49.077;
52.409;]; 

RNA.ertugay = [ 9.142;
9.193;
9.810;
8.908;
8.953;
9.857;
11.524;];


%% Data: Yu2020 

% strain: CEN.PK113-7D
% ref: Nitrogen limitation reveals large reserves inmetabolic and translational capacities of yeast
% C/N = 4  ; glucsoe: 7.5g/L (NH4)2SO4 5g/L
% C/N = 30 ; glucsoe: 7.5g/L (NH4)2SO4 0.5g/L
% C/N = 50 ; glucsoe: 7.5g/L (NH4)2SO4 0.29g/L
% C/N = 115; glucsoe: 7.5g/L (NH4)2SO4 0.125g/L

paper_info.gl_in.Yu4    = 7.5;
paper_info.n_in.Yu4     = 5;
paper_info.CN.Yu4       = 2.2* 7.5/5;
paper_info.author.Yu4   = 'Yu2020_CN4'; 
paper_info.strain.Yu4   = 'CEN.PK113-7D';

D.Yu4 = 0.2; 

% unit: g/L
glucose.Yu4 = [0];

% unit: mmol/gDW/h
% converted units: uM/h 
Jgy.Yu4 = [2.44];%_%[2.42];
Jgy.Yu4  = Jgy.Yu4  * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);

% unit: mmol/gDW/h
% converted units: uM/h 
Jeh.Yu4 = [0];%_%[0.02];
Jeh.Yu4  = Jeh.Yu4  * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);

% unit:g/L
%
nh4.Yu4 = [0];%_%[0.016]; 
%nh4.Yu4  = 10^6 * nh4.Yu4 /18.04;

% unit: % 
protein.Yu4 = [47.7];%_%[47.84]; 


prot_sec.r.Yu4       = [0.27959];   %_%[0.3985       ];
prot_sec.z.Yu4       = [0.25351];   %_%[0.22813      ];
prot_sec.e.Yu4       = [0.46691];   %_%[0.37337      ];
prot_sec.gy.Yu4      = [0.10289];   %_%[0.08879      ];
prot_sec.fe.Yu4      = [0.018489];  %_%[0.012544     ];
prot_sec.gn.Yu4      = [0.044799];  %_%[0.029583     ];
prot_sec.mt.Yu4      = [0.078295];  %_%[0.086935     ];
prot_sec.sp.Yu4      = [0.001215];  %_%[0.00064336   ];
prot_sec.sd.Yu4      = [0.00020516];%_%[8.9573e-05   ];
prot_sec.as.Yu4      = [0.23247];   %_%[0.16322      ];
prot_sec.at.Yu4      = [0.00050738];%_%[0.00029371   ];
prot_sec.lp.Yu4      = [0.095782];  %_%[0.061796     ];
prot_sec.lo.Yu4      = [0.0027728]; %_%[0.0019844    ];


%% 
% % Suarez-Mendez, C. A., Hanemaaijer, M., ten Pierick, A., Wolters, J. C., Heijnen, J. J., & Wahl, S. A. (2016). 
% % Interaction of storage carbohydrates and other cyclic fluxes with central metabolism: 
% % A quantitative approach by non-stationary 13C metabolic flux analysis. Metabolic Engineering Communications, 3(65), 52–63. 
D.Suarez_Mendez = [0.054 0.101 0.207 0.307];

% unit umol/gDW
% convert to % 
carbo.sc.Suarez_Mendez = [435.480; 593.791; 274.755; 237.925;] +  [146.203; 162.025; 7.595; 5.696] ;
carbo.sc.Suarez_Mendez = (10^(-6))*carbo.sc.Suarez_Mendez* 180*100; 


end 
