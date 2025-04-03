

%clear; clc; close all; 
%plt_clrs = plot_colors(); 
data_unit_conv;

% data format: [dilution rate, concentration]

%% CEN.PK 
% chemostat
% Xia 2021 (Xia_protein)  --------------------------------------
% Proteome allocations change linearly with 
% specific growth rate of Saccharomyces cerevisiae under glucose-limitation

xia.dr =   [0.027; 
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
xia.Jgy =  [0.399;
0.553;
1.14;
1.66;
2.40;
2.74;
3.20;
4.88;
6.78;];

xia.Jgy = xia.Jgy * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);

% Kumar 2019 - low glucose condition --------------------------------------

% original units: 1/h
kumar_lg.dr = [0.05;
0.12;
0.18;
0.24;
0.31];

% original units: dcw (g/L)       
% converted units: cells/L_culture
kumar_lg.cells = [0.32;
0.44;
0.34;
0.33;
0.30];  

kumar_lg.cells = kumar_lg.cells * gdw_to_cell;                   

% qs                   
% original units: g/gdw/h         
% converted units: uM glucose/h
kumar_lg.Jgy = [0.16;
0.27;
0.53;
0.73;
1.03];     

kumar_lg.Jgy = kumar_lg.Jgy * one_over_gdw_to_1_over_cell * (1/par.V_c) * gl_gPerL_to_uM;                 

% qEtOH
% original units: g/gdw/h     
% converted units: uM ethanol/h
kumar_lg.Jfe_minus_Jgo = [0;
0;
0;
0;
0];   

kumar_lg.Jfe_minus_Jgo = kumar_lg.Jfe_minus_Jgo * one_over_gdw_to_1_over_cell * (1/par.V_c) * eh_gPerL_to_uM;                                                      

% original units: umol/gdw       
% converted units: uM 
kumar_lg.g6p = [1.0645163061355;
1.8072881736487;
2.14282912259462;
2.98231109116389;
1.63578591111586];                      

kumar_lg.g6p = kumar_lg.g6p * one_over_gdw_to_1_over_cell * (1/par.V_c);    

% original units: umol/gdw       
% converted units: uM 
kumar_lg.pyruvate = [3.325;
4.871;
8.989;
8.326;
5.893];        

kumar_lg.pyruvate = kumar_lg.pyruvate * one_over_gdw_to_1_over_cell * (1/par.V_c);                                           

% original units: umol/gdw        
% converted units: uM 
kumar_lg.atp = [5.347;
7.041;
8.852;
9.231;
6.684];  

kumar_lg.atp = kumar_lg.atp * one_over_gdw_to_1_over_cell * (1/par.V_c);                                 

% original units: umol/gdw         
% converted units: uM 
kumar_lg.glutamate = [236.8;
366.4;
413.4;
426.0;
324.7];     

kumar_lg.glutamate = kumar_lg.glutamate * one_over_gdw_to_1_over_cell * (1/par.V_c);                                             

% original units: umol/gdw      
% converted units: uM 
kumar_lg.glutamine = [40.51;
70.30;
88.44;
111.0;
108.3];         

kumar_lg.glutamine = kumar_lg.glutamine * one_over_gdw_to_1_over_cell * (1/par.V_c);                       

% original units: mmol/gdw/h      
% converted units: uM/h 
kumar_lg.co2_prod = [2.28;
3.95;
7.52;
8.30;
10.10]; 

kumar_lg.co2_prod = kumar_lg.co2_prod * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);  


% original units: mmol/gdw/h      
% converted units: uM/h 
kumar_lg.o2_uptake = [2.42;
3.47;
7.51;
8.11;
10.45;]; 

kumar_lg.o2_uptake = kumar_lg.o2_uptake * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c); 

% Kumar 2019 - high glucose condition -------------------------------------

% original units: 1/h
kumar_hg.dr = [0.12;
0.26;
0.36;
0.41];

% original units: dcw (g/L)         
% converted units: cells/L_culture
kumar_hg.cells = [4.25;
4.40;
1.71;
0.55];    

kumar_hg.cells = kumar_hg.cells * gdw_to_cell; 

% qs                   
% original units: g/gdw/h   
% converted units: uM glucose/h
kumar_hg.Jgy = [0.28;
0.59;
2.05;
3.76];     

kumar_hg.Jgy = kumar_hg.Jgy * one_over_gdw_to_1_over_cell * (1/par.V_c) * gl_gPerL_to_uM; 

% qEtOH
% original units: g/gdw/h    
% converted units: uM ethanol/h
kumar_hg.Jfe_minus_Jgo = [0;
0;
0.57;
1.31];   

kumar_hg.Jfe_minus_Jgo = kumar_hg.Jfe_minus_Jgo * one_over_gdw_to_1_over_cell * (1/par.V_c) * eh_gPerL_to_uM;                            

% original units: umol/gdw       
% converted units: uM 
kumar_hg.g6p = [1.59750985308395;
1.65884964663491;
0.606088883651831;
1.65780339285714];                      

kumar_hg.g6p = kumar_hg.g6p * one_over_gdw_to_1_over_cell * (1/par.V_c);   

% original units: umol/gdw  
% converted units: uM 
kumar_hg.pyruvate = [3.882;
4.427;
3.68;
7.420];           

kumar_hg.pyruvate = kumar_hg.pyruvate * one_over_gdw_to_1_over_cell * (1/par.V_c);                      

%kumar_hg.pyruvate_batch_gl = 4.7495203942839;
%kumar_hg.pyruvate_batch_gl = kumar_hg.pyruvate_batch_gl * one_over_gdw_to_1_over_cell * (1/par.V_c);                      

% original units: umol/gdw      
% converted units: uM 
kumar_hg.atp = [6.291;
6.970;
5.517;
6.292];  

kumar_hg.atp = kumar_hg.atp * one_over_gdw_to_1_over_cell * (1/par.V_c);                 

%kumar_hg.atp_batch_gl = 3.060351241;
%kumar_hg.atp_batch_gl = kumar_hg.atp_batch_gl *  one_over_gdw_to_1_over_cell * (1/par.V_c);   

% original units: umol/gdw     
% converted units: uM 
kumar_hg.glutamate = [237.9;
245.7;
124.8;
106.5];     

kumar_hg.glutamate = kumar_hg.glutamate * one_over_gdw_to_1_over_cell * (1/par.V_c);                       

% original units: umol/gdw      
% converted units: uM 
kumar_hg.glutamine = [47.53;
111.8;
195.5;
293.3];      

kumar_hg.glutamine = kumar_hg.glutamine * one_over_gdw_to_1_over_cell * (1/par.V_c);      

% original units: mmol/gdw/h      
% converted units: uM/h 
kumar_hg.co2_prod = [3.50;
5.84;
19.26;
26.68]; 

kumar_hg.co2_prod = kumar_hg.co2_prod * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);   


% original units: mmol/gdw/h      
% converted units: uM/h 
kumar_hg.o2_uptake = [3.35;
5.6;
10.14;
9.37]; 

kumar_hg.o2_uptake = kumar_hg.o2_uptake * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c); 

%kumar_hg.glutamine_batch_gl = 135.449336398504;                           
%kumar_hg.glutamine_batch_gl = kumar_hg.glutamine_batch_gl * one_over_gdw_to_1_over_cell * (1/par.V_c);  

% 1999 diderich -----------------------------------------------------------
% glucose feed = 7.5 g/L 

% original units: h^-1
diderich.dr = [0.0537209302;
0.11;
0.1637209302;
0.22;
0.2762790698;
0.3069767442;
0.33;
0.363255814;
0.3965116279;
0.4195348837;
];

% original units: gdw/g      
% converted units: cells/umol
diderich.Ysx = [0.489583333;
0.5;
0.5;
0.5;
0.479166667;
0.479166667;
0.41875;
0.270833333;
0.2125;
0.172916667;
];

diderich.Ysx = diderich.Ysx * (1/(one_over_gdw_to_1_over_cell)) * (1/gl_gPerL_to_uM);            

% original units: mmol/gdw/h
% converted units: uM/h
diderich.Jfe_minus_Jgo = [0.00E+00;
0.00E+00;
0.00E+00;
0.00E+00;
0.00E+00;
0.00E+00;
1.458333333;
6.458333333;
12.29166667;
13.125;
]; 

diderich.Jfe_minus_Jgo = diderich.Jfe_minus_Jgo * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);

% 2001 van Maris ----------------------------------------------------------
% glucose feed = 7.5 g/L

% original units: h^-1
vanmaris.dr = [0.1;
0.148917458;
0.199728269;
0.249487217;
0.279783784;
0.300406136;
0.330536158;
0.349595325;
0.360426589;
0.380645727;
];

% original units: gx/gs
% converted units: cells/umol
vanmaris.Ysx = [0.479459459;
0.5;
0.5;
0.489189189;
0.479459459;
0.449189189;
0.301081081;
0.249189189;
0.241621622;
0.160540541;
];

vanmaris.Ysx = vanmaris.Ysx * (1/(one_over_gdw_to_1_over_cell)) * (1/gl_gPerL_to_uM);                        

% original units: mmol/gdw/h
% converted units: uM/h
vanmaris.Jfe_minus_Jgo = [0;
0;
0.00E+00;
0.00E+00;
0.00E+00;
0.507246377;
4.673913043;
6.014492754;
6.847826087;
18.29710145;
];

vanmaris.Jfe_minus_Jgo = vanmaris.Jfe_minus_Jgo * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);

% 1998 vanhoek ------------------------------------------------------------
% glucose feed = 7.5 g/L

% original units: h^-1
vanhoek.dr = [0.1;
0.15;
0.2;
0.247727273;
0.277272727;
0.297727273;
0.327272727;
0.347727273;
0.356818182;
0.377272727;
];

% original units: gx/gs
% converted units: cells/umol          
vanhoek.Ysx = [0.479591837;
0.502040816;
0.502040816;
0.489795918;
0.479591837;
0.451020408;
0.3;
0.251020408;
0.240816327;
0.16122449;
];

vanhoek.Ysx = vanhoek.Ysx * (1/(one_over_gdw_to_1_over_cell)) * (1/gl_gPerL_to_uM);                        

% original units: mmol/gdw/h
% converted units: uM/h           
vanhoek.Jfe_minus_Jgo = [0.00E+00;
0.00E+00;
0.00E+00;
0;
0.00E+00;
0.612244898;
4.693877551;
6.020408163;
6.734693878;
18.36734694;
];

vanhoek.Jfe_minus_Jgo = vanhoek.Jfe_minus_Jgo * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);

% 2000 bakker -------------------------------------------------------------
% glucose feed = 7.5 g/L

% original units: h^-1
bakker.dr_ysx = [0.048639554;
0.098143038;
0.148299977;
0.197835278;
0.24798977;
0.278393774;
0.298391449;
0.328580082;
0.358484817;
0.379042946;
];

% original units: gx/gs
% converted units: cells/umol               
bakker.Ysx = [0.490215922;
0.500590434;
0.499899045;
0.499996941;
0.50009606;
0.479602788;
0.420353771;
0.26942529;
0.210195853;
0.169920276;
];

bakker.Ysx = bakker.Ysx * (1/(one_over_gdw_to_1_over_cell)) * (1/gl_gPerL_to_uM);                        

% original units: h^-1          
bakker.dr_Jfe_minus_Jgo = [0.099690402;
0.149845201;
0.19876161;
0.249535604;
0.279256966;
0.299071207;
0.329411765;
0.359133127;
0.379566563;
];

% original units: mmol/gdw/h
% converted units: uM/h                           
bakker.Jfe_minus_Jgo = [0
0
0
0
0
1.487210674
6.512826183
12.44316674
13.3253354
];

bakker.Jfe_minus_Jgo = bakker.Jfe_minus_Jgo * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);

%% FY4 

% Boer 2010 - glucose feed = 0.8 g/L

% original units: h^-1
boer.dr = [0.057692308
0.108333333
0.160897436
0.226923077
0.310897436
];

% original units: M/h       
% converted units: uM/h
boer.Jgy = [0.246575342
4.52E-01
6.44E-01
8.49E-01
3.287671233
];       
boer.Jgy = boer.Jgy .* 10^6;         

% original units: M/h        
boer.Jfe_minus_Jgo = [0.00E+00
0.00E+00
0.00E+00
0.01369863
3.712328767
];   
boer.Jfe_minus_Jgo = boer.Jfe_minus_Jgo .* 10^6;                                       

% original units: *10^7 cells/mL   
% converted units: cells/L
boer.cells = [3.507846608
4.063205506
4.680550639
4.265526057
1.102261554
];
boer.cells = boer.cells .* 10^7 * 1000;   

% fold change relative to value at lowest dilution rate
boer.pyr = [1
1.439018711
1.760144149
3.020963998
4.187441002]; 

boer.f16bp = [1
1.058587296
1.348334043
1.126903374
14.32443253]; 

boer.glutamine = [1
1.439018711
1.760144149
3.020963998
4.187441002]; 

boer.atp = [1
1.309703238
1.163730571
1.022060572
1.1615593];             

% Hackett 2016 - glucose feed = 0.8 g/L

% original units: h^-1
hackett.dr = [0.049863024
0.105431457
0.154377453
0.212650311
0.29384141
];

% original units: moles / hr / mL cells
hackett.Jgy = [0.000149923
0.00031825
0.0004743
0.001085886
0.001886523];

hackett.Jgy = hackett.Jgy *10^9;

hackett.Jfe_minus_Jgo = [2.5945E-06
1.09181E-05
2.36939E-05
0.000976161
0.002583142]; 

hackett.Jfe_minus_Jgo = hackett.Jfe_minus_Jgo *10^9;

%  original units: volume fraction (mL cells / L culture)
hackett.cells = [1.443564631
1.432924105
1.390362003
0.819320466
0.535573119];

hackett.cells =  hackett.cells *10^-3./par.V_c;    
% original units:(fL)
hackett.cells_v =  [19.7
18.4
21
21
23.5];

%hackett.cells =  hackett.cells *10^-3./(hackett.cells_v* 10^-15);    

% original units: M
% converted units: uM
hackett.atp = [9.86E-04
1.03E-03
1.06E-03
1.02E-03
0.0008924
];
hackett.atp = hackett.atp * 10^6;           

% original units: M 
hackett.glutamine = [2.68E-03
2.07E-03
2.21E-03
3.11E-03
0.003401716
];
hackett.glutamine = hackett.glutamine .* 10^6;                  

% original units: M                 
hackett.pyruvate = [1.82E-03
2.01E-03
2.18E-03
2.48E-03
0.003439365
];                 
hackett.pyruvate = hackett.pyruvate .* 10^6;   

hackett.g6p = [0.0006946
0.001079212
0.001333398
0.0010472
0.001732866
];

hackett.g6p = hackett.g6p .* 10^6; 

%
hackett.gdw_per_mL = [0.204
0.205
0.218
0.211
0.187
]; 
%}

%hackett.gdw_per_mL = 2.38;                   

% original units: g/gdw    
% converted units: uM
hackett.glycogen = [0.053261229
0.048197953
0.054575099
0.056255205
0.057700495
];

hackett.glycogen = hackett.glycogen * gl_gPerL_to_uM .* hackett.gdw_per_mL * 1000;

% original units: g/gdw    
% converted units: uM                  
hackett.trehalose = [0.096895892
0.071429784
0.048027297
0.03663382
0.011485232
];        

hackett.trehalose = hackett.trehalose * gl_gPerL_to_uM .* hackett.gdw_per_mL * 1000;                 
hackett.sc = hackett.glycogen + hackett.trehalose; 

% yifei 2021/10/05 added 
% pyruvate data --------------------------------------------
% strain: ?CBS 8066
postma.dr         =  [0.1
0.149
0.2
0.249
0.299
0.322
0.344
0.359
0.379
0.402
0.392
0.398
0.424
0.474
0.486
0.433
0.414
0.45
0.452];
% unit mM 
postma.pyruvate_ex =  [0.014
0.01
0.006
0.029
0.014
0.06
0.1
0.163
0.276
0.328
0.497
0.48
1.515
1.533
1.714
1.803
1.814
2.001
2.199];

% Sonia Cortassa, Miguel A. Aon, Distributed control of the glycolytic flux 
% in wild-type cells and catabolite repression mutants of Saccharomyces cerevisiae 
% growing in carbon-limited chemostat cultures, Enzyme and Microbial Technology, 1997
% strain: CEN.Pk122                   
cortassa.dr =   [0.1 
0.225 
0.3];

% unit mM
cortassa.pyruvate =   [0.1 
0.1 
5.98];  

% strain: CEN.PK 113-5D with pYX212-EFE (ethylene forming enzyme)                

johansson.dr         = [0.03
0.1
0.15
0.2
0.25
0.3
0.35];
% unit: g/L 
johansson.pyruvate_ex = [0
1.06
0.11
0.64
3.19
17.77
24.04];         
%% plot by glucose feed concentration 
%{
x_label = 'dilution rate (h^{-1})';
y_lim_jgy           = [0 10^7];
y_lim_Jfe_minus_Jgo = [0 1.2*10^7];
y_lim_cells         = [0 2.8*10^11];
y_lim_aa            = [0 1.2*10^5];
y_lim_atp           = [0 4.0*10^3];
y_lim_pc1           = [0 2.0*10^3];
y_lim_pc2           = [0 4.0*10^3];
y_lim_sc            = [0 1.8*10^5];
a = 8;
b = 4; 
figure; 

% 0.8 g/L 
subplot(a, b, 1)
plot(boer.dr, boer.Jgy, 'o', 'Color', plt_clrs.blue)
ylabel('J_{ht} (\mu M/h)')
%xlabel(x_label)
xlim([0 0.45])
ylim(y_lim_jgy)
%legend({'boer'})
box on; 
axis square; 
title('[G_{l}^{f}] = 0.8 g/L')

% ethanol production rate
subplot(a, b, 1 + b)
plot(boer.dr, boer.Jfe_minus_Jgo, 'o', 'Color', plt_clrs.blue)
ylabel({'J_{fe} - J_{gn}'; '(\mu M/h)'})
%xlabel(x_label)
xlim([0 0.45])
ylim(y_lim_Jfe_minus_Jgo)
%legend({'boer'})
box on; 
axis square; 

% biomass
subplot(a, b, 1 + 2*b)
plot(boer.dr, boer.cells, 'o', 'Color', plt_clrs.blue)
ylabel('cells')
%xlabel(x_label)
xlim([0 0.45])
ylim(y_lim_cells)
%legend({'boer'})
box on; 
axis square; 

% pyruvate
subplot(a, b, 1 + 3*b)
plot(hackett.dr, hackett.g6p, 'o', 'Color', plt_clrs.blue)
ylabel({'glucose_{in}'; '(\mu M)'; '(g6p)'})
%xlabel(x_label)
xlim([0 0.45])
ylim(y_lim_pc1)
%legend({'hackett'})
box on; 
axis square; 

% pyruvate
subplot(a, b, 1 + 4*b)
plot(hackett.dr, hackett.pyruvate, 'o', 'Color', plt_clrs.blue)
ylabel({'precursor'; '(\mu M)'; '(pyruvate)'})
%xlabel(x_label)
xlim([0 0.45])
ylim(y_lim_pc2)
%legend({'hackett'})
box on; 
axis square; 

% glutamine
subplot(a, b, 1 + 5*b)
plot(hackett.dr, hackett.glutamine, 'o', 'Color', plt_clrs.blue)
ylabel({'amino acids'; '(\mu M)'; '(glutamine)'})
%xlabel(x_label)
xlim([0 0.45])
ylim(y_lim_aa)
%legend({'hackett'})
box on; 
axis square; 

% atp
subplot(a, b, 1 + 6*b)
plot(hackett.dr, hackett.atp, 'o', 'Color', plt_clrs.blue)
ylabel('atp (\mu M)')
%xlabel(x_label)
xlim([0 0.45])
ylim(y_lim_atp)
%legend({'hackett'})
box on; 
axis square; 

% storage carbohydrates
subplot(a, b, 1 + 7*b)
plot(hackett.dr, hackett.glycogen + hackett.trehalose, 'o', 'Color', plt_clrs.blue)
ylabel({'storage'; 'carbohydrates'; '(\mu M)'})
%ylabel('glycogen + trehalose')
xlabel(x_label)
xlim([0 0.45])
ylim(y_lim_sc)
%legend({'hackett'})
box on; 
axis square;

% 1 g/L -------------------------------------------------------------------

% glucose uptake rate
subplot(a, b, 2)
plot(kumar_lg.dr, kumar_lg.Jgy, 'o', 'Color', plt_clrs.green)
ylabel('J_{ht} (\mu M/h)')
%xlabel(x_label)
xlim([0 0.45])
ylim(y_lim_jgy)
%legend({'kumar'})
box on; 
axis square; 
title('[G_{l}]^{f} = 1 g/L')

% ethanol production rate
subplot(a, b, 2 + b)
plot(kumar_lg.dr, kumar_lg.Jfe_minus_Jgo, 'o', 'Color', plt_clrs.green)
ylabel({'J_{fe} - J_{gn}'; '(\mu M/h)'})
%xlabel(x_label)
xlim([0 0.45])
ylim(y_lim_Jfe_minus_Jgo)
%legend({'kumar'})
box on; 
axis square; 

% biomass
subplot(a, b, 2 + 2*b)
plot(kumar_lg.dr, kumar_lg.cells, 'o', 'Color', plt_clrs.green)
ylabel('cells')
%xlabel(x_label)
xlim([0 0.45])
ylim(y_lim_cells)
%legend({'kumar'})
box on; 
axis square; 

% g6p
subplot(a, b, 2 + 3*b)
plot(kumar_lg.dr, kumar_lg.g6p, 'o', 'Color', plt_clrs.green)
ylabel({'glucose_{in}'; '(\mu M)'; '(g6p)'})
%xlabel(x_label)
xlim([0 0.45])
ylim(y_lim_pc1)
%legend({'kumar'})
box on; 
axis square; 

% pyruvate
subplot(a, b, 2 + 4*b)
plot(kumar_lg.dr, kumar_lg.pyruvate, 'o', 'Color', plt_clrs.green)
ylabel({'precursor'; '(\mu M)'; '(pyruvate)'})
%xlabel(x_label)
xlim([0 0.45])
ylim(y_lim_pc2)
%legend({'kumar'})
box on; 
axis square; 

% glutamine
subplot(a, b, 2 + 5*b)
plot(kumar_lg.dr, kumar_lg.glutamine, 'o', 'Color', plt_clrs.green)
ylabel({'amino acids'; '(\mu M)'; '(glutamine)'})
%xlabel(x_label)
xlim([0 0.45])
ylim(y_lim_aa)
%legend({'kumar'})
box on; 
axis square; 

% atp
subplot(a, b, 2 + 6*b)
plot(kumar_lg.dr, kumar_lg.atp, 'o', 'Color', plt_clrs.green)
ylabel('atp (\mu M)')
xlabel(x_label)
xlim([0 0.45])
ylim(y_lim_atp)
%legend({'kumar'})
box on; 
axis square; 

% 10 g/L ------------------------------------------------------------------

% glucose uptake rate
subplot(a, b, 4)
plot(kumar_hg.dr, kumar_hg.Jgy, 'o', 'Color', plt_clrs.green)
ylabel('J_{ht} (\mu M/h)')
%xlabel(x_label)
xlim([0 0.45])
ylim(y_lim_jgy)
%legend({'kumar'})
box on; 
axis square; 
title('[G_{l}]^{f} = 10 g/L')

% ethanol production rate
subplot(a, b, 4 + b)
plot(kumar_hg.dr, kumar_hg.Jfe_minus_Jgo, 'o', 'Color', plt_clrs.green)
ylabel({'J_{fe} - J_{gn}'; '(\mu M/h)'})
%xlabel(x_label)
xlim([0 0.45])
ylim(y_lim_Jfe_minus_Jgo)
%legend({'kumar'})
box on; 
axis square; 

% biomass
subplot(a, b, 4 + 2*b)
plot(kumar_hg.dr, kumar_hg.cells, 'o', 'Color', plt_clrs.green)
ylabel('cells')
%xlabel(x_label)
xlim([0 0.45])
ylim(y_lim_cells)
%legend({'kumar'})
box on; 
axis square; 

% g6p
subplot(a, b, 4 + 3*b)
plot(kumar_hg.dr, kumar_hg.g6p, 'o', 'Color', plt_clrs.green)
ylabel({'glucose_{in}'; '(\mu M)'; '(g6p)'})
%xlabel(x_label)
xlim([0 0.45])
ylim(y_lim_pc1)
%legend({'kumar'})
box on; 
axis square; 

% pyruvate
subplot(a, b, 4 + 4*b)
plot(kumar_hg.dr, kumar_hg.pyruvate, 'o', 'Color', plt_clrs.green)
ylabel({'precursor'; '(\mu M)'; '(pyruvate)'})
%xlabel(x_label)
xlim([0 0.45])
ylim(y_lim_pc2)
%legend({'kumar'})
box on; 
axis square; 

% glutamine
subplot(a, b, 4 + 5*b)
plot(kumar_hg.dr, kumar_hg.glutamine, 'o', 'Color', plt_clrs.green)
ylabel({'amino acids'; '(\mu M)'; '(glutamine)'})
%xlabel(x_label)
xlim([0 0.45])
ylim(y_lim_aa)
%legend({'kumar'})
box on; 
axis square; 

% atp
subplot(a, b, 4 + 6*b)
plot(kumar_hg.dr, kumar_hg.atp, 'o', 'Color', plt_clrs.green)
ylabel('atp (\mu M)')
xlabel(x_label)
xlim([0 0.45])
ylim(y_lim_atp)
%legend({'kumar'})
box on; 
axis square; 

% 7.5 g/L
% -------------------------------------------------------------------------

subplot(a, b, 3)
axis off;
title('[G_{l}]^{f} = 7.5 g/L')

% biomass yield 
subplot(a, b, 3 + 2*b)
hold on; 
plot(diderich.dr,     diderich.Ysx, 'o', 'Color', plt_clrs.green)
plot(vanmaris.dr,     vanmaris.Ysx, 's', 'Color', plt_clrs.green)
plot(vanhoek.dr,      vanhoek.Ysx,  '^', 'Color', plt_clrs.green)
plot(bakker.dr_ysx,   bakker.Ysx,   'v', 'Color', plt_clrs.green)
ylabel({'biomass yield'; '(cells/ glucose'; 'consumed)'})
xlabel(x_label)
xlim([0 0.45])
%legend({'diderich'; 'van Maris'; 'van Hoek'; 'bakker'})
hold off; 
box on; 
axis square; 

% Jfe minus 0.5 Jgn
subplot(a, b, 3 + b)
hold on; 
plot(diderich.dr,             diderich.Jfe_minus_Jgo, 'o', 'Color', plt_clrs.green)
plot(vanmaris.dr,             vanmaris.Jfe_minus_Jgo, 's', 'Color', plt_clrs.green)
plot(vanhoek.dr,              vanhoek.Jfe_minus_Jgo,  '^', 'Color', plt_clrs.green)
plot(bakker.dr_Jfe_minus_Jgo, bakker.Jfe_minus_Jgo,   'v', 'Color', plt_clrs.green)
ylabel({'J_{fe} - J_{gn}'; '(\mu M/h)'})
%xlabel(x_label)
xlim([0 0.45])
%legend({'diderich'; 'van Maris'; 'van Hoek'; 'bakker'})
hold off; 
box on; 
axis square; 

subplot(a, b, 3+4*b)
hold on; 
plot(0, 0, '^', 'Color', plt_clrs.green) % van Hoek
plot(0, 0, 'x', 'Color', plt_clrs.green) % diderich
plot(0, 0, 'v', 'Color', plt_clrs.green) % Bakker
plot(0, 0, 's', 'Color', plt_clrs.green) % van Maris
plot(0, 0, '*', 'Color', plt_clrs.blue) % Boer
plot(0, 0, 'o', 'Color', plt_clrs.blue) % Hackett


legend({'1998 van Hoek, CEN.PK113-7D';
'1999 Diderich, CEN.PK113-7D';
'2000 Bakker, CEN.PK113-7D';
'2001 van Maris, CEN.PK113-7D';
'2010 Boer, FY4';
'2016 Hackett, FY4'})
legend boxoff; 

hold off; 
axis off; 

%sgtitle('Chemostat metabolite data - organized by glucose feed concentration - colored by strain')
%}
