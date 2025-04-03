function paper_plot_sc_atp_lerp(D, D_crit_simu, ...
WT_y_steady_clim, WT_y_steady_nlim, WT_y_steady_cnlim1, WT_y_steady_cnlim2, ...
WT_met_reac_steady_clim, WT_met_reac_steady_nlim, WT_met_reac_steady_cnlim1, WT_met_reac_steady_cnlim2,...
WT_sig_steady_clim, WT_sig_steady_nlim, WT_sig_steady_cnlim1, WT_sig_steady_cnlim2, ...
par, num_y, num_flux, num_prot ,figure_output_format, figure_position)
yeast_type = 'sc'; 
%% load data
data_unit_conv; 


[D_clim, cell_clim, glucose_clim, Jgy_clim, Jeh_clim,... 
precursor_clim, atp_clim, aa_clim,  prot_sec_clim, ...
lipid_clim, protein_clim,carbo_clim, RNA_clim, paperinfo_clim,...
color_clim, shape_clim, gas_clim] = SC_C_lim_chem();   

[D_nlim, cell_nlim, glucose_nlim, Jgy_nlim, Jeh_nlim, nh4_nlim, gas_nlim, ... 
precursor_nlim, atp_nlim, aa_nlim,  prot_sec_nlim, ...
lipid_nlim, protein_nlim,carbo_nlim, RNA_nlim, paperinfo_nlim,...                    
color_nlim, shape_nlim, Jnh4_nlim]  = SC_N_lim_chem();  

lipid_constant_part = lipid_sc_c; 


%% plotting stuff
%biomass_log_climale = 1;
paper_figure_color; 
paper_figure_label; 


figure_output = figure_output_format;

[b, a]   = size(figure_output); 
position = figure_position;

x_label  = 'dilution rate (h^{-1})'; 

y_lim_prot = [-3 1.3]; % log2 fold change range 
x_lim      = [0  0.42]; 
ylim_RNA   = [0 15]; 
ylim_prot  = [0 57];   % biomass range 
ylim_lipid = [0 50]; 
ylim_carbo = [0 45];

% Patch area 
x_left   = 0.21; % h-1 experiment critical dilution lower bound
x_middle = 0.33;  % h-1 experiment critical dilution upper bound

color_clim.default =  [102, 194, 165]./255; % green = color_clim.kumar_hg;
color_nlim.default  = [31,  120, 180]./255 ; % blue = color_clim.hackett; 

if strcmp(yeast_type,'sc')

shape_clim.kumar_hg       = 'o'; 
shape_clim.xia            = '^';    
shape_clim.hackett        = 'v'; 
shape_clim.kumar_lg       = 'd'; 
shape_clim.boer           = 's';
shape_clim.Suarez_Mendez  = 's';
shape_clim.lange          = 'v'; 
shape_clim.ertugay        = 'v'; 
shape_clim.Yu4            = 'o';     

plt_clrs.darkgreen   = '#419945';
color_cnlim1.default = '#51c5b8';%'#368E6A';
color_cnlim2.default = '#0095dd';%'#2B838F';
%_% plt_clrs.green   = '#61c95e';'#1a9c53';
%_% plt_clrs.blue   = '#0e80be';'#1f8cbe';
%_% color_cnlim1.default = '#38b07b';'#4cbfb0';'#56afb2';
%_% color_cnlim2.default = '#3eb0c8';'#36a4c5';'#3aa9d1';

color_clim.default    = plt_clrs.green; %_% darkgreen;
color_clim.kumar_hg   = plt_clrs.green;
color_clim.xia        = plt_clrs.green;
color_clim.kumar_lg   = plt_clrs.green;
color_clim.Suarez_Mendez  = plt_clrs.green;
color_clim.boer       = plt_clrs.green;
color_clim.boer       = plt_clrs.green;

color_clim.lange     = plt_clrs.green;    
color_clim.ertugay   = plt_clrs.green;
color_clim.Yu4       = plt_clrs.green;





% nlim
shape_nlim.hackett    = shape_clim.hackett; 
shape_nlim.kumar_ln   = shape_clim.kumar_lg; 
shape_nlim.Yu30       = shape_clim.Yu4;  
shape_nlim.Yu50       = shape_clim.Yu4;     
shape_nlim.Yu115      = shape_clim.Yu4;       
shape_nlim.Yu2021     = shape_clim.xia;   
shape_nlim.boer       = shape_clim.boer;     

color_nlim.default    = plt_clrs.blue;


color_nlim.Yu30       = color_cnlim1.default;    
color_nlim.Yu50       = color_cnlim1.default;     
color_nlim.Yu2021     = color_cnlim1.default;    

color_nlim.Yu115      = color_cnlim2.default;   


color_nlim.kumar_ln   = color_cnlim2.default;     
color_nlim.hackett    = plt_clrs.blue;     
color_nlim.boer       = plt_clrs.blue;   

end


if strcmp(yeast_type,'rt')
plt_clrs.darkgreen    = '#419945';
color_clim.default    = plt_clrs.darkgreen;   
color_clim.default    = plt_clrs.green;
color_clim.shen       = plt_clrs.green;
color_clim.shen2017   = plt_clrs.green;    
color_clim.Dinh       = plt_clrs.green; 
color_clim.Zhu        = plt_clrs.green; 


shape_clim.shen       = 'o'; 
shape_clim.shen2017   = '^';  
shape_clim.Dinh       = 's';     
shape_clim.Zhu        = 'v';  

color_nlim.default    = plt_clrs.blue;
color_nlim.shen       = plt_clrs.blue;
color_nlim.shen2017   = plt_clrs.blue;    
color_nlim.Dinh       = plt_clrs.blue; 
color_nlim.Zhu       = plt_clrs.blue; 

shape_nlim.shen       = 'o'; 
shape_nlim.shen2017   = '^';  
shape_nlim.Dinh       = 's'; 
shape_nlim.Zhu        = 'v';    
end 

%% proteins frac (sim)
total_protein_con_clim = sum(WT_y_steady_clim(:,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);
total_protein_con_nlim = sum(WT_y_steady_nlim(:,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);
total_protein_con_cnlim1 = sum(WT_y_steady_cnlim1(:,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);
total_protein_con_cnlim2 = sum(WT_y_steady_cnlim2(:,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);

%_% WT_y_steady_clim(:,[num_y.r: num_y.sd]) = (WT_y_steady_clim(:,[num_y.r: num_y.sd]).*par.l')./ total_protein_con_clim;
%_% WT_y_steady_nlim(:,[num_y.r: num_y.sd]) = (WT_y_steady_nlim(:,[num_y.r: num_y.sd]).*par.l')./ total_protein_con_nlim;
%_% WT_y_steady_cnlim1(:,[num_y.r: num_y.sd]) = (WT_y_steady_cnlim1(:,[num_y.r: num_y.sd]).*par.l')./ total_protein_con_cnlim1;
%_% WT_y_steady_cnlim2(:,[num_y.r: num_y.sd]) = (WT_y_steady_cnlim2(:,[num_y.r: num_y.sd]).*par.l')./ total_protein_con_cnlim2;
WT_y_steady_clim(:,[num_y.r: num_y.lo]) = (WT_y_steady_clim(:,[num_y.r: num_y.lo]).*par.l')./ total_protein_con_clim;
WT_y_steady_nlim(:,[num_y.r: num_y.lo]) = (WT_y_steady_nlim(:,[num_y.r: num_y.lo]).*par.l')./ total_protein_con_nlim;
WT_y_steady_cnlim1(:,[num_y.r: num_y.lo]) = (WT_y_steady_cnlim1(:,[num_y.r: num_y.lo]).*par.l')./ total_protein_con_cnlim1;
WT_y_steady_cnlim2(:,[num_y.r: num_y.lo]) = (WT_y_steady_cnlim2(:,[num_y.r: num_y.lo]).*par.l')./ total_protein_con_cnlim2;

WT_y_steady_clim_frac = WT_y_steady_clim;
WT_y_steady_nlim_frac = WT_y_steady_nlim;
WT_y_steady_cnlim1_frac = WT_y_steady_cnlim1;
WT_y_steady_cnlim2_frac = WT_y_steady_cnlim2;



%% REZ fraction calculation before normalizing
% total_protein_con_clim = sum(WT_y_steady_clim(:,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);
% total_protein_con_nlim = sum(WT_y_steady_nlim(:,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);
% total_protein_con_cnlim1 = sum(WT_y_steady_cnlim1(:,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);
% total_protein_con_cnlim2 = sum(WT_y_steady_cnlim2(:,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);
% 
% WT_y_steady_clim(:,[num_y.r: num_y.sd]) = (WT_y_steady_clim(:,[num_y.r: num_y.sd]).*par.l')./ total_protein_con_clim;
% WT_y_steady_nlim(:,[num_y.r: num_y.sd]) = (WT_y_steady_nlim(:,[num_y.r: num_y.sd]).*par.l')./ total_protein_con_nlim;
% WT_y_steady_cnlim1(:,[num_y.r: num_y.sd]) = (WT_y_steady_cnlim1(:,[num_y.r: num_y.sd]).*par.l')./ total_protein_con_cnlim1;
% WT_y_steady_cnlim2(:,[num_y.r: num_y.sd]) = (WT_y_steady_cnlim2(:,[num_y.r: num_y.sd]).*par.l')./ total_protein_con_cnlim2;
% 
% WT_y_steady_clim_frac = WT_y_steady_clim;
% WT_y_steady_nlim_frac = WT_y_steady_nlim;
% WT_y_steady_cnlim1_frac = WT_y_steady_cnlim1;
% WT_y_steady_cnlim2_frac = WT_y_steady_cnlim2;

% fraction
R_chemo_sim_clim = WT_y_steady_clim(:,num_y.r).*par.l(num_prot.r)' ./ total_protein_con_clim;
Z_chemo_sim_clim = WT_y_steady_clim(:,num_y.z).*par.l(num_prot.z)' ./ total_protein_con_clim;
E_chemo_sim_clim = sum(WT_y_steady_clim(:,num_y.gy:num_y.lo).*par.l(num_prot.gy:num_prot.lo)',2) ./ total_protein_con_clim;
% 
% R_chemo_sim_nlim = WT_y_steady_nlim(:,num_y.r).*par.l(num_prot.r)' ./ total_protein_con_nlim;
% Z_chemo_sim_nlim = WT_y_steady_nlim(:,num_y.z).*par.l(num_prot.z)' ./ total_protein_con_nlim;
% E_chemo_sim_nlim = sum(WT_y_steady_nlim(:,num_y.gy:num_y.lo).*par.l(num_prot.gy:num_prot.lo)',2) ./ total_protein_con_nlim;


% R_chemo_sim_clim = WT_y_steady_clim(:,num_y.r);
% Z_chemo_sim_clim = WT_y_steady_clim(:,num_y.z);
% E_chemo_sim_clim = sum(WT_y_steady_clim(:,num_y.gy:num_y.lo),2);

R_chemo_sim_nlim = WT_y_steady_nlim(:,num_y.r).*par.l(num_prot.r)' ;
Z_chemo_sim_nlim = WT_y_steady_nlim(:,num_y.z).*par.l(num_prot.z)' ;
E_chemo_sim_nlim = sum(WT_y_steady_nlim(:,num_y.gy:num_y.lo), 2);

% relative fold change
R_chemo_sim_fc = R_chemo_sim_clim./R_chemo_sim_clim(1); 
Z_chemo_sim_fc = Z_chemo_sim_clim./Z_chemo_sim_clim(1);
E_chemo_sim_fc = E_chemo_sim_clim./E_chemo_sim_clim(1);

R_chemo_sim_fc_nlim = R_chemo_sim_nlim./R_chemo_sim_nlim(1); 
Z_chemo_sim_fc_nlim = Z_chemo_sim_nlim./Z_chemo_sim_nlim(1);
E_chemo_sim_fc_nlim = E_chemo_sim_nlim./E_chemo_sim_nlim(1);




%% normalize precursor, amino acid, atp, protein with respect to first data point 
D_exp = 0.05; % first D in exp
[norm_D, norm_indx] = min(abs(D - D_exp)); 
% aa atp (sim)
WT_y_steady_clim(:,[num_y.aa_in, num_y.ae]) = WT_y_steady_clim(:,[num_y.aa_in, num_y.ae]) ./ WT_y_steady_clim(norm_indx,[num_y.aa_in, num_y.ae]);
WT_y_steady_nlim(:,[num_y.aa_in, num_y.ae]) = WT_y_steady_nlim(:,[num_y.aa_in, num_y.ae]) ./ WT_y_steady_nlim(norm_indx,[num_y.aa_in, num_y.ae]);
WT_y_steady_cnlim1(:,[num_y.aa_in, num_y.ae]) = WT_y_steady_cnlim1(:,[num_y.aa_in, num_y.ae]) ./ WT_y_steady_cnlim1(norm_indx,[num_y.aa_in, num_y.ae]);
WT_y_steady_cnlim2(:,[num_y.aa_in, num_y.ae]) = WT_y_steady_cnlim2(:,[num_y.aa_in, num_y.ae]) ./ WT_y_steady_cnlim2(norm_indx,[num_y.aa_in, num_y.ae]);

% precursor (sim)
%norm_indx = numel(WT_y_steady_clim(:,[num_y.pc])); 
WT_y_steady_clim(:,[num_y.pc]) = WT_y_steady_clim(:,[num_y.pc]) ./ WT_y_steady_clim(norm_indx,[num_y.pc]);
WT_y_steady_nlim(:,[num_y.pc]) = WT_y_steady_nlim(:,[num_y.pc]) ./ WT_y_steady_nlim(norm_indx,[num_y.pc]);
WT_y_steady_cnlim1(:,[num_y.pc]) = WT_y_steady_cnlim1(:,[num_y.pc]) ./ WT_y_steady_cnlim1(norm_indx,[num_y.pc]);
WT_y_steady_cnlim2(:,[num_y.pc]) = WT_y_steady_cnlim2(:,[num_y.pc]) ./ WT_y_steady_cnlim2(norm_indx,[num_y.pc]);



% protein log2fc 
%_% WT_y_steady_clim(:,[num_y.r: num_y.sd]) = log2(WT_y_steady_clim(:,[num_y.r: num_y.sd])./WT_y_steady_clim(norm_indx,[num_y.r: num_y.sd])); 
%_% WT_y_steady_nlim(:,[num_y.r: num_y.sd]) = log2(WT_y_steady_nlim(:,[num_y.r: num_y.sd])./WT_y_steady_nlim(norm_indx,[num_y.r: num_y.sd])); 
%_% WT_y_steady_cnlim1(:,[num_y.r: num_y.sd]) = log2(WT_y_steady_cnlim1(:,[num_y.r: num_y.sd])./WT_y_steady_cnlim1(norm_indx,[num_y.r: num_y.sd])); 
%_% WT_y_steady_cnlim2(:,[num_y.r: num_y.sd]) = log2(WT_y_steady_cnlim2(:,[num_y.r: num_y.sd])./WT_y_steady_cnlim2(norm_indx,[num_y.r: num_y.sd])); 
WT_y_steady_clim(:,[num_y.r: num_y.lo]) = log2(WT_y_steady_clim(:,[num_y.r: num_y.lo])./WT_y_steady_clim(norm_indx,[num_y.r: num_y.lo])); 
WT_y_steady_nlim(:,[num_y.r: num_y.lo]) = log2(WT_y_steady_nlim(:,[num_y.r: num_y.lo])./WT_y_steady_nlim(norm_indx,[num_y.r: num_y.lo])); 
WT_y_steady_cnlim1(:,[num_y.r: num_y.lo]) = log2(WT_y_steady_cnlim1(:,[num_y.r: num_y.lo])./WT_y_steady_cnlim1(norm_indx,[num_y.r: num_y.lo])); 
WT_y_steady_cnlim2(:,[num_y.r: num_y.lo]) = log2(WT_y_steady_cnlim2(:,[num_y.r: num_y.lo])./WT_y_steady_cnlim2(norm_indx,[num_y.r: num_y.lo])); 

%% plotting SC -------------------------------------------------------------

%% ATP flux  --------------------------------------------------------------
x = 0:0.1:0.5; 
Jatp_cal_clim.kumar_hg = 1.4*Jgy_clim.kumar_hg + 9*gas_clim.o2.kumar_hg/3; 
Jatp_cal_clim.kumar_lg = 1.4*Jgy_clim.kumar_lg + 9*gas_clim.o2.kumar_lg/3; 
Jatp_cal_nlim.kumar_ln = 1.4*Jgy_nlim.kumar_ln + 9*gas_nlim.o2.kumar_ln/3; 


Jatp_gy_cal_clim.kumar_hg = 1.4*Jgy_clim.kumar_hg./Jatp_cal_clim.kumar_hg; 
Jatp_gy_cal_clim.kumar_lg = 1.4*Jgy_clim.kumar_lg./Jatp_cal_clim.kumar_lg;
Jatp_gy_cal_nlim.kumar_ln = 1.4*Jgy_nlim.kumar_ln./Jatp_cal_nlim.kumar_ln; 

Jatp_mt_cal_clim.kumar_hg  = 1 - Jatp_gy_cal_clim.kumar_hg; 
Jatp_mt_cal_clim.kumar_lg  = 1 - Jatp_gy_cal_clim.kumar_lg; 
Jatp_mt_cal_nlim.kumar_ln  = 1 - Jatp_gy_cal_nlim.kumar_ln; 

JATP_fit_clim.kumar_lc        = 0.9*1.98E+07.*x;               % carbon-limited fitting result %_% manual linear fitting, need change
JATP_fit_nlim.kumar_ln        = 0.9*2.4E+07.*x + 459484*1.5;     % nitrogen-limited fitting result %_% manual linear fitting, need change

%_% Jatp_cal_clim.xia      = 1.4*Jgy_clim.xia + 9*gas_clim.o2.xia/3; 
%_% Jatp_cal_nlim.Yu2021   = 1.4*Jgy_nlim.Yu2021 + 9*gas_nlim.o2.Yu2021/3; 

%_% JATP_fit_clim.xia      = 0.9*1.87E+07.*x;               % carbon-limited fitting result
%_% JATP_fit_nlim.Yu2021   = 2.42E+07.*x + 459484*1;     % nitrogen-limited fitting result

JATP_clim   =  par.q_gy*  WT_met_reac_steady_clim.flux(:, num_flux.gy)   + par.q_mt*  WT_met_reac_steady_clim.flux(:, num_flux.mt);  
JATP_nlim   =  par.q_gy*  WT_met_reac_steady_nlim.flux(:, num_flux.gy)   + par.q_mt*  WT_met_reac_steady_nlim.flux(:, num_flux.mt);  
JATP_cnlim1 =  par.q_gy*  WT_met_reac_steady_cnlim1.flux(:, num_flux.gy) + par.q_mt*  WT_met_reac_steady_cnlim1.flux(:, num_flux.mt);  
JATP_cnlim2 =  par.q_gy*  WT_met_reac_steady_cnlim2.flux(:, num_flux.gy) + par.q_mt*  WT_met_reac_steady_cnlim2.flux(:, num_flux.mt);  

JATP_gy_clim = par.q_gy*  WT_met_reac_steady_clim.flux(:, num_flux.gy)./JATP_clim;
JATP_gy_nlim = par.q_gy*  WT_met_reac_steady_nlim.flux(:, num_flux.gy)./JATP_nlim;
JATP_gy_cnlim1 = par.q_gy*  WT_met_reac_steady_cnlim1.flux(:, num_flux.gy)./JATP_cnlim1;
JATP_gy_cnlim2 = par.q_gy*  WT_met_reac_steady_cnlim2.flux(:, num_flux.gy)./JATP_cnlim2;

JATP_as_clim   = (par.q_as + par.q_p) *  WT_met_reac_steady_clim.flux(:, num_flux.as); 
JATP_as_cnlim1 = (par.q_as + par.q_p) *  WT_met_reac_steady_cnlim1.flux(:, num_flux.as); 
JATP_as_nlim   = (par.q_as + par.q_p) *  WT_met_reac_steady_nlim.flux(:, num_flux.as); 

JATP_lp_clim   = par.q_lp *  WT_met_reac_steady_clim.flux(:, num_flux.lp_fe); 
JATP_lp_cnlim1 = par.q_lp *  WT_met_reac_steady_cnlim1.flux(:, num_flux.lp_fe);
JATP_lp_nlim   = par.q_lp *  WT_met_reac_steady_nlim.flux(:, num_flux.lp_fe);

JATP_mt_clim = 1 - JATP_gy_clim    ;
JATP_mt_nlim = 1 - JATP_gy_nlim    ;
JATP_mt_cnlim1 = 1- JATP_gy_cnlim1 ;
JATP_mt_cnlim2 = 1 - JATP_gy_cnlim2;


figure; 
a = 2; b = 4; 
x_lim = [0, 0.41]; 

subplot(a,b,1)
hold on 
plot(D, JATP_clim, '-', 'Color', color_clim.default, 'LineWidth',2);
plot(D, JATP_nlim, '-', 'Color', color_nlim.default, 'LineWidth',2);
hold off

ylabel('J_{ATP} (\muM)', 'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
ylim([0 10^7])
yticks([0 5 10]*10^6)%_%
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'NH4')), 'fontweight', 'bold', 'FontSize', fontsize);
axis square; 
box on;
hold off;



% subplot(a,b,4)
% hold on 
% plot(D, JATP_clim, '-', 'Color', color_clim.default);
% plot(D, JATP_cnlim1, '-', 'Color', color_nlim.default);
% hold off
% 
% ylabel('J_{ATP} (\muM)', 'Color',xy_label_color)
% xlabel(x_label,'Color',xy_label_color);
% xlim(x_lim);
% ylim([0 1.1*max(ylim)])
% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
% %text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'NH4')), 'fontweight', 'bold', 'FontSize', fontsize);
% axis square; 
% box on;
% hold off;

% subplot(a,b,4+b)
% hold on 
% plot(D, JATP_clim, '-', 'Color', color_clim.default);
% plot(D, JATP_cnlim2, '-', 'Color', color_nlim.default);
% hold off
% 
% ylabel('J_{ATP} (\muM)', 'Color',xy_label_color)
% xlabel(x_label,'Color',xy_label_color);
% xlim(x_lim);
% ylim([0 1.1*max(ylim)])
% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
% %text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'NH4')), 'fontweight', 'bold', 'FontSize', fontsize);
% axis square; 
% box on;
% hold off;


subplot(a,b,1+b)
hold on 
%plot(D_clim.kumar_lg, Jatp_cal_clim.kumar_lg , 'o','MarkerFaceColor',plt_clrs.blue)
%_% plot(D_clim.kumar_hg, Jatp_cal_clim.kumar_hg ,  'o','Color',color_clim.default,  'MarkerFaceColor',color_clim.default)
%_% plot(D_nlim.kumar_ln, Jatp_cal_nlim.kumar_ln,   'o','Color',color_nlim.default,  'MarkerFaceColor',color_nlim.default)
plot(D_clim.kumar_hg, Jatp_cal_clim.kumar_hg ,  'o-','Color',color_clim.default,  'MarkerFaceColor',color_clim.default, 'LineWidth',2)
plot(D_nlim.kumar_ln, Jatp_cal_nlim.kumar_ln,   'o-','Color',color_nlim.default,  'MarkerFaceColor',color_nlim.default, 'LineWidth',2)

%_% plot(x, JATP_fit_clim.kumar_lc , '-','Color',color_clim.default)
%_% plot(x, JATP_fit_nlim.kumar_ln , '-','Color',color_nlim.default)
hold off
xlabel('D (h^{-1})')
ylabel('ATP prod. (uM/h)')
ylim([0 10^7])
yticks([0 5 10]*10^6)%_%
xlim(x_lim)
axis square 
box on 
legend('Clim', 'Nlim')
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
sgtitle('data from kumar')
set(gcf,'Position',position)


subplot(a,b,2)
hold on 
plot(D, JATP_gy_clim, '-', 'Color', color_clim.default, 'LineWidth',2);
plot(D, JATP_mt_clim, '--', 'Color', color_clim.default, 'LineWidth',2);
hold off

ylabel('J_{ATP} (\muM)', 'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
ylim([0 1])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'NH4')), 'fontweight', 'bold', 'FontSize', fontsize);
axis square; 
box on;
hold off;



subplot(a,b,3)
hold on 
%_% plot(D, JATP_gy_cnlim1, '-', 'Color', color_nlim.default);
%_% plot(D, JATP_mt_cnlim1, '--', 'Color', color_nlim.default);
plot(D, JATP_gy_cnlim2, '-', 'Color', color_nlim.default, 'LineWidth',2);
plot(D, JATP_mt_cnlim2, '--', 'Color', color_nlim.default, 'LineWidth',2);
hold off

ylabel('J_{ATP} (\muM)', 'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
ylim([0 1])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'NH4')), 'fontweight', 'bold', 'FontSize', fontsize);
axis square; 
box on;
hold off;


subplot(a,b,4)
hold on 
plot(D, JATP_as_clim, '-', 'Color', color_nlim.default, 'LineWidth',2);
plot(D, JATP_lp_clim, '--', 'Color', color_nlim.default, 'LineWidth',2);
hold off

ylabel('J_{ATP} (\muM)', 'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
ylim([0 9*10^6])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'NH4')), 'fontweight', 'bold', 'FontSize', fontsize);
axis square; 
box on;
hold off;

subplot(a,b,4+b)
hold on 
plot(D, JATP_as_nlim, '-', 'Color', color_nlim.default, 'LineWidth',2);
plot(D, JATP_lp_nlim, '--', 'Color', color_nlim.default, 'LineWidth',2);
hold off

ylabel('J_{ATP} (\muM)', 'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
ylim([0 9*10^6])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'NH4')), 'fontweight', 'bold', 'FontSize', fontsize);
axis square; 
box on;
hold off;


subplot(a,b,2+b)
hold on 
%plot(D_clim.kumar_lg, Jatp_cal_clim.kumar_lg , 'o','MarkerFaceColor',plt_clrs.blue)
plot(D_clim.kumar_hg, Jatp_gy_cal_clim.kumar_hg ,  '--^','Color',color_clim.default,  'MarkerFaceColor',color_clim.default, 'LineWidth',2)
plot(D_clim.kumar_hg, Jatp_mt_cal_clim.kumar_hg,   '-o','Color',color_clim.default,  'MarkerFaceColor',color_clim.default, 'LineWidth',2)
hold off
xlabel('D (h^{-1})')
ylabel('ATP prod. (uM/h)')
ylim([0 1])
xlim(x_lim)
axis square 
box on 
legend('Clim', 'Nlim')
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
sgtitle('data from kumar')
set(gcf,'Position',position)


subplot(a,b,3+b)
hold on 
%plot(D_clim.kumar_lg, Jatp_cal_clim.kumar_lg , 'o','MarkerFaceColor',plt_clrs.blue)
plot(D_nlim.kumar_ln, Jatp_gy_cal_nlim.kumar_ln ,  '--^','Color',color_nlim.default,  'MarkerFaceColor',color_nlim.default, 'LineWidth',2)
plot(D_nlim.kumar_ln, Jatp_mt_cal_nlim.kumar_ln,   '-o','Color',color_nlim.default,  'MarkerFaceColor',color_nlim.default, 'LineWidth',2)
hold off
xlabel('D (h^{-1})')
ylabel('ATP prod. (uM/h)')
ylim([0 1])
xlim(x_lim)
axis square 
box on 
legend('Clim', 'Nlim')
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)

set(gcf,'Position',position)

set(gcf, 'Position', [101 492 753 337])




% if exist('prot_sec_nlim','var')
%     papers = fieldnames(prot_sec_nlim.r); 
%     for i = 1: length(papers)
%         plot(D_nlim.(papers{i}),  log2(prot_sec_nlim.r.(papers{i})./prot_sec_nlim.r.(papers{i})(1)) ,...
%             shape_nlim.(papers{i}),...
%             'Color', color_nlim.(papers{i}), ...
%             'MarkerFaceColor', color_nlim.(papers{i}))
%     end 
% end 
