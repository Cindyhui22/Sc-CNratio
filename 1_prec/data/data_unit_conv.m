
%{
% yifei
par.V_c                  =  4.2*10^-14; 
vol_per_cell                =  4.2*10^-14;       %(4.2*10^-14 L/cell)
gdw_per_cell                =  1.9E-11;        % used this number to convert  
gdw_to_cell                 =  1/1.9E-11; %1/(1.65*10^(-11));
one_over_gdw_to_1_over_cell =  1.9E-11;  % 1.65*10^(-11);  % 1.62*10^-11 g/cell
gl_gPerL_to_uM              =  5550.621669627;
eh_gPerL_to_uM              =  21706.0994139353;  
mmol_to_umol                =  1000; 
OD_to_gdw                   =  0.62;  % gDW/L
gdw_per_L                   =  gdw_per_cell./vol_per_cell;
par.V_c                  = 4.2E-14; 
%}

%
% V
fr_gPerL_to_uM              = 5550.621669627;
su_gPerL_to_uM              = 2921.44400666931;
gdw_to_cell                 = 1/(7.170731694E-12);
one_over_gdw_to_1_over_cell = 7.170731694E-12; 
gdw_per_cell                = one_over_gdw_to_1_over_cell; 
gdw_per_L                   = 170.731707; 
gl_gPerL_to_uM              = 5550.621669627;
eh_gPerL_to_uM              = 21706.0994139353; 
nh4_gPerL_to_uM             = 7567.7;     % (nh4)2so4  132.14 g/mol
mmol_to_umol                = 1000; 

lp_uM_to_gPerL              = 842.3/(10^6);%836.5/(10^6); %_%- 31.9/(10^6);
pt_uM_to_gPerL              = 110/(10^6); 
gl_uM_to_gPerL              = 180/(10^6); 


cell_volume  = 4.2E-14; % L
vol_per_cell = 4.2E-14; % L
par.V_c   = 4.2E-14; 
OD_to_gdw = 0.486;       % 0.486 gdw/OD

hour_to_min = 60; 

% https://bionumbers.hms.harvard.edu/bionumber.aspx?id=100987&ver=5
OD660_to_cells =  1.85e+7 * 1000;  % cell/ml * 1000 = cell/L

lipid_sc_c = 2.4; 
lipid_rt_c = 8; 
%}
