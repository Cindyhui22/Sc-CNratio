function paper_plot_sc_met_lerp(D, D_crit_simu, ...
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

color_clim.default    = plt_clrs.green; %_%darkgreen;
color_clim.kumar_hg   = plt_clrs.green;
color_clim.xia        = plt_clrs.green;
color_clim.kumar_lg   = plt_clrs.green;
color_clim.Suarez_Mendez  = plt_clrs.green;
color_clim.boer       = plt_clrs.green;
color_clim.hackett    = plt_clrs.green;

color_clim.lange     = plt_clrs.green;    
color_clim.ertugay   = plt_clrs.green;
color_clim.Yu4       = plt_clrs.green;





% nlim
shape_nlim.hackett    = shape_clim.hackett; %'v'
shape_nlim.kumar_ln   = shape_clim.kumar_lg; %'d'
shape_nlim.Yu30       = shape_clim.Yu4;  %'o'
shape_nlim.Yu50       = shape_clim.Yu4;     
shape_nlim.Yu115      = shape_clim.Yu4;       
shape_nlim.Yu2021     = shape_clim.xia;  %'^'   
shape_nlim.boer       = shape_clim.boer; %'s'    

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


%% real metabolites
figure;
% glucose uptake rate
subplot(a, b, find(strcmp(figure_output,'glucose')))
hold on;
if exist('Jgy_clim','var') 
papers = fieldnames(Jgy_clim);
for i = 1: length(papers)
plot(D_clim.(papers{i}),  Jgy_clim.(papers{i}) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end    
end
if exist('Jgy_nlim','var') 
papers = fieldnames(Jgy_nlim); 
for i = 1: length(papers)
plot(D_nlim.(papers{i}),  Jgy_nlim.(papers{i}) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end
plot(D, WT_met_reac_steady_clim.flux(:,num_flux.gy), '-', 'Color', color_clim.default, 'LineWidth',2);
plot(D, WT_met_reac_steady_nlim.flux(:,num_flux.gy), '-', 'Color', color_nlim.default, 'LineWidth',2);
plot(D, WT_met_reac_steady_cnlim1.flux(:,num_flux.gy), '-', 'Color', color_cnlim1.default, 'LineWidth',2);
plot(D, WT_met_reac_steady_cnlim2.flux(:,num_flux.gy), '-', 'Color', color_cnlim2.default, 'LineWidth',2);

ylabel(figcabbi.y_label.glu,'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
%max_y = max([max(ylim), max(kumar_lg.Jgy), max(kumar_hg.Jgy), max(boer.Jgy), max(hackett.Jgy)]); 
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'glucose')), 'fontweight', 'bold', 'FontSize', fontsize);
yticks([0 2 4]*10^6)
ylim([0 5*10^6])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
axis square; 
box on;
hold off;

% ethanol production rate
%
subplot(a, b, find(strcmp(figure_output,'ethanol')))
hold on;
if exist('Jeh_clim','var') 
papers = fieldnames(Jeh_clim);
for i = 1: length(papers)
plot(D_clim.(papers{i}),  Jeh_clim.(papers{i}) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end 
if exist('Jeh_nlim','var') 
papers = fieldnames(Jeh_nlim);
for i = 1: length(papers)
plot(D_nlim.(papers{i}),  Jeh_nlim.(papers{i}) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end 
%plot(D, WT_met_reac_steady_1gL.flux(:,num_flux.fe)  - WT_met_reac_steady_1gL.flux(:,num_flux.gn), '-', 'Color',  color_clim.kumar_lg);
plot(D, WT_met_reac_steady_clim.flux(:,num_flux.fe) - WT_met_reac_steady_clim.flux(:,num_flux.gn), '-', 'Color',  color_clim.default, 'LineWidth',2);
plot(D, WT_met_reac_steady_nlim.flux(:,num_flux.fe) - WT_met_reac_steady_nlim.flux(:,num_flux.gn), '-', 'Color',  color_nlim.default, 'LineWidth',2);
plot(D, WT_met_reac_steady_cnlim1.flux(:,num_flux.fe) - WT_met_reac_steady_cnlim1.flux(:,num_flux.gn), '-', 'Color',  color_cnlim1.default, 'LineWidth',2);
plot(D, WT_met_reac_steady_cnlim2.flux(:,num_flux.fe) - WT_met_reac_steady_cnlim2.flux(:,num_flux.gn), '-', 'Color',  color_cnlim2.default, 'LineWidth',2);
%max_y = max([max(ylim), max(kumar_lg.Jfe_minus_Jgo), max(kumar_hg.Jfe_minus_Jgo), max(boer.Jfe_minus_Jgo), max(hackett.Jfe_minus_Jgo)]); 

ylabel(figcabbi.y_label.eth,'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
yticks([0 3*10^6 6*10^6])
ylim([0 6*10^6])

%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'ethanol')), 'fontweight', 'bold', 'FontSize', fontsize);
axis square; 
box on;
hold off;
%}

%{
% nh4 uptake rate
subplot(a, b, find(strcmp(figure_output,'JNH4')))
hold on;
if exist('Jnh4_clim','var') 
papers = fieldnames(Jnh4_clim);
for i = 1: length(papers)
plot(D_clim.(papers{i}),  Jnh4_clim.(papers{i}) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end 


plot(D, par.n_nh *WT_met_reac_steady_clim.flux(:,num_flux.as), '-', 'Color', color_clim.default);
plot(D, par.n_nh *WT_met_reac_steady_nlim.flux(:,num_flux.as), '-', 'Color', color_nlim.default);
plot(D, par.n_nh *WT_met_reac_steady_cnlim1.flux(:,num_flux.as), '-', 'Color', color_cnlim1.default);
plot(D, par.n_nh *WT_met_reac_steady_cnlim2.flux(:,num_flux.as), '-', 'Color', color_cnlim2.default);

ylabel(figcabbi.y_label.jnh4,'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
%max_y = max([max(ylim), max(kumar_lg.Jgy), max(kumar_hg.Jgy), max(boer.Jgy), max(hackett.Jgy)]); 
patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'glucose')), 'fontweight', 'bold', 'FontSize', fontsize);

ylim([0 max(ylim)])
axis square; 
box on;
hold off;
%}

% biomass
subplot(a, b, find(strcmp(figure_output,'cells')))
hold on;
if exist('cell_clim','var') 
papers = fieldnames(cell_clim);
for i = 1: length(papers)
plot(D_clim.(papers{i}),  cell_clim.(papers{i}) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end 
if exist('cell_nlim','var') 
papers = fieldnames(cell_nlim); 
for i = 1: length(papers)
plot(D_nlim.(papers{i}),  cell_nlim.(papers{i}) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end 
%plot(D, WT_y_steady_1gL(:,end), '-', 'Color', color_clim.kumar_lg);
plot(D, WT_y_steady_clim(:,end), '-', 'Color', color_clim.default, 'LineWidth',2);
plot(D, WT_y_steady_nlim(:,end), '-', 'Color', color_nlim.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim1(:,end), '-', 'Color', color_cnlim1.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim2(:,end), '-', 'Color', color_cnlim2.default, 'LineWidth',2);
ylabel(figcabbi.y_label.cell,'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
ylim([0 6.6])
yticks([0 3 6])
% if biomass_log_climale == 0
%     yticks([0 3*10^11 6*10^11])
%     ylim([0 7.4*10^11])
% end 
% if biomass_log_climale == 1
%     ylim([1*10^-1 1*10^1])
%     set(gca,'Yscale','log')
%     set(gca,'Yminortick','off')
% end 
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'cells')), 'fontweight', 'bold', 'FontSize', fontsize);
axis square; 
box on;
hold off;

%% virtual metabolites 

% precursor 
%
subplot(a, b, find(strcmp(figure_output,'precursor')))
hold on;
if exist('precursor_clim','var') 
precs = fieldnames(precursor_clim);
for k = 1: length(precs)
prec = precs{k};
papers = fieldnames(precursor_clim.(prec)); 
for i = 1: length(papers)
plot(D_clim.(papers{i}),  precursor_clim.(prec).(papers{i})./precursor_clim.(prec).(papers{i})(1) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end 
end 
if exist('precursor_nlim','var') 
precs = fieldnames(precursor_nlim);
for k = 1: length(precs)
prec = precs{k};
papers = fieldnames(precursor_nlim.(prec)); 
for i = 1: length(papers)
plot(D_nlim.(papers{i}),  precursor_nlim.(prec).(papers{i})./precursor_nlim.(prec).(papers{i})(1) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end 
end 
% plot(postma.dr,   postma.pyruvate_ex,       shape_postma,    'MarkerFaceColor', color_postma,    'MarkerEdgeColor', 'k')
% plot(cortassa.dr, cortassa.pyruvate,        shape_cortassa,  'MarkerFaceColor', color_cortassa,    'MarkerEdgeColor', 'k')
% plot(johansson.dr, johansson.pyruvate_ex,   shape_johansson, 'MarkerFaceColor', color_johansson,    'MarkerEdgeColor', 'k')

%plot(D, WT_y_steady_1gL(:,num_y.pc),  '-', 'Color', color_clim.kumar_lg);
plot(D, WT_y_steady_clim(:,num_y.pc), '-', 'Color', color_clim.default, 'LineWidth',2);
plot(D, WT_y_steady_nlim(:,num_y.pc), '-', 'Color', color_nlim.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim1(:,num_y.pc), '-', 'Color', color_cnlim1.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim2(:,num_y.pc), '-', 'Color', color_cnlim2.default, 'LineWidth',2);

ylabel(figcabbi.y_label.pc,'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
ylim([0 max(ylim)*1.05])
yticks([0 20 40])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'precursor')),'fontweight','bold','FontSize',fontsize);
axis square; 
box on;
hold off;
%}

% intracellular amino acids 
%
subplot(a, b, find(strcmp(figure_output,'aa_in')))
hold on;
if exist('aa_clim','var') 
papers = fieldnames(aa_clim.glutamine); 
for i = 1: length(papers)
plot(D_clim.(papers{i}),  aa_clim.glutamine.(papers{i})./aa_clim.glutamine.(papers{i})(1) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end 

if exist('aa_nlim','var') 
papers = fieldnames(aa_nlim.glutamine); 
for i = 1: length(papers)
plot(D_nlim.(papers{i}),  aa_nlim.glutamine.(papers{i})./aa_nlim.glutamine.(papers{i})(1) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end 


%plot(D, WT_y_steady_1gL(:,num_y.aa_in),  '-', 'Color', color_clim.kumar_lg);
plot(D, WT_y_steady_clim(:,num_y.aa_in), '-', 'Color', color_clim.default, 'LineWidth',2);
plot(D, WT_y_steady_nlim(:,num_y.aa_in), '-', 'Color', color_nlim.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim1(:,num_y.aa_in), '-', 'Color', color_cnlim1.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim2(:,num_y.aa_in), '-', 'Color', color_cnlim2.default, 'LineWidth',2);


ylabel(figcabbi.y_label.aa,'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
yticks([0 5 10])%_%
%_% ylim([0 1.05*max(ylim)])
ylim([0 10.5])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'aa_in')),'fontweight','bold','FontSize',fontsize);
axis square; 
box on;
hold off;
%}



% lipid 
%
% subplot(a, b, find(strcmp(figure_output,'lipid')))
% hold on;
% 
% %plot(D, WT_y_steady_1gL(:,num_y.ae),  '-', 'Color', color_clim.kumar_lg);
% plot(D, WT_y_steady_clim(:,num_y.lp_e), '-', 'Color', color_clim.default);
% plot(D, WT_y_steady_nlim(:,num_y.lp_e), '-', 'Color', color_nlim.default);
% 
% ylabel(figcabbi.y_label.lipid,'Color',xy_label_color)
% xlabel(x_label,'Color',xy_label_color);
% xlim(x_lim);
% ylim([0 1.1*max(ylim)])
% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
% %text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'lipid')), 'fontweight', 'bold', 'FontSize', fontsize);
% axis square; 
% box on;
% hold off;



% NH4 
%
subplot(a, b, find(strcmp(figure_output,'NH4')))
hold on;
if exist('nh4_clim','var') 
papers = fieldnames(nh4_clim);
for i = 1: length(papers)
plot(D_clim.(papers{i}),  nh4_clim.(papers{i}) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end 
if exist('nh4_nlim','var') 
papers = fieldnames(nh4_nlim); 
for i = 1: length(papers)
plot(D_nlim.(papers{i}),  nh4_nlim.(papers{i}) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end 
%plot(D, WT_y_steady_1gL(:,num_y.ae),  '-', 'Color', color_clim.kumar_lg);
plot(D, WT_y_steady_clim(:,num_y.nh4), '-', 'Color', color_clim.default, 'LineWidth',2);
plot(D, WT_y_steady_nlim(:,num_y.nh4), '-', 'Color', color_nlim.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim1(:,num_y.nh4), '-', 'Color', color_cnlim1.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim2(:,num_y.nh4), '-', 'Color', color_cnlim2.default, 'LineWidth',2);

ylabel(figcabbi.y_label.nh4,'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
yticks([0 2 4]*10^5);%_%
ylim([0 1.1*max(ylim)])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'NH4')), 'fontweight', 'bold', 'FontSize', fontsize);
axis square; 
box on;
hold off;

%}

%% signaling 
%{\
% SNF1
subplot(a, b, find(strcmp(figure_output,'snf1')))
hold on;
%plot(D, WT_sig_steady_1gL.snf1(:,1),  '-', 'Color', color_clim.kumar_lg);
%_% plot(D, WT_sig_steady_clim.snf1(:,1), '-', 'Color', color_clim.default);
%_% plot(D, WT_sig_steady_nlim.snf1(:,1), '-', 'Color', color_nlim.default);
%_% plot(D, WT_sig_steady_cnlim1.snf1(:,1), '-', 'Color', color_cnlim1.default);
%_% plot(D,  WT_sig_steady_cnlim2.snf1(:,1), '-', 'Color', color_cnlim2.default);
plot(D, WT_sig_steady_clim.snf1(:,1)./par.s_tot, '-', 'Color', color_clim.default, 'LineWidth',2);
plot(D, WT_sig_steady_nlim.snf1(:,1)./par.s_tot, '-', 'Color', color_nlim.default, 'LineWidth',2);
plot(D, WT_sig_steady_cnlim1.snf1(:,1)./par.s_tot, '-', 'Color', color_cnlim1.default, 'LineWidth',2);
plot(D,  WT_sig_steady_cnlim2.snf1(:,1)./par.s_tot, '-', 'Color', color_cnlim2.default, 'LineWidth',2);

ylabel(figcabbi.y_label.snf1,'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
yticks([0 0.5 1])
ylim([0 1.08])

%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'snf1')), 'fontweight', 'bold', 'FontSize', fontsize);
axis square; 
box on;
hold off; 


% TORC1
subplot(a, b, find(strcmp(figure_output,'tor')))
hold on;
%plot(D, WT_sig_steady_1gL.tor(:,1),  '-', 'Color', color_clim.kumar_lg);
plot(D, WT_sig_steady_clim.tor(:,1)./par.tau_tot, '-', 'Color', color_clim.default, 'LineWidth',2);
plot(D, WT_sig_steady_nlim.tor(:,1)./par.tau_tot, '-', 'Color', color_nlim.default, 'LineWidth',2);
plot(D, WT_sig_steady_cnlim1.tor(:,1)./par.tau_tot, '-', 'Color', color_cnlim1.default, 'LineWidth',2);
plot(D,  WT_sig_steady_cnlim2.tor(:,1)./par.tau_tot, '-', 'Color', color_cnlim2.default, 'LineWidth',2);
hold off; 

ylabel(figcabbi.y_label.tor,'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
ylim([0.35 0.6])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'tor')), 'fontweight', 'bold', 'FontSize', fontsize);
axis square; 
box on;

%}

%% Jeh vs Jgy

subplot(a, b, find(strcmp(figure_output,'JehJgy')))
hold on 
if exist('Jeh_clim','var') 
papers = fieldnames(Jeh_clim);
for i = 1: length(papers)
plot(Jgy_clim.(papers{i}),  Jeh_clim.(papers{i}) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end 
if exist('Jeh_nlim','var') 
papers = fieldnames(Jeh_nlim);
for i = 1: length(papers)
plot(Jgy_nlim.(papers{i}),  Jeh_nlim.(papers{i}) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end 
plot(WT_met_reac_steady_clim.flux(:,num_flux.gy), WT_met_reac_steady_clim.flux(:,num_flux.fe)- WT_met_reac_steady_clim.flux(:,num_flux.gn), '-', 'Color', color_clim.default, 'LineWidth',2);
plot(WT_met_reac_steady_nlim.flux(:,num_flux.gy), WT_met_reac_steady_nlim.flux(:,num_flux.fe)- WT_met_reac_steady_nlim.flux(:,num_flux.gn), '-', 'Color', color_nlim.default, 'LineWidth',2);
plot(WT_met_reac_steady_cnlim1.flux(:,num_flux.gy),  WT_met_reac_steady_cnlim1.flux(:,num_flux.fe) - WT_met_reac_steady_cnlim1.flux(:,num_flux.gn), '-', 'Color', color_cnlim1.default, 'LineWidth',2);
plot(WT_met_reac_steady_cnlim2.flux(:,num_flux.gy),  WT_met_reac_steady_cnlim2.flux(:,num_flux.fe) - WT_met_reac_steady_cnlim2.flux(:,num_flux.gn), '-', 'Color', color_cnlim2.default, 'LineWidth',2);

hold off
xticks([0 2 4]*10^6)%_%
yticks([0 2 4 6]*10^6)%_%
ylim([0 max(ylim)*1.05])

%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'precursor')),'fontweight','bold','FontSize',fontsize);
axis square; 
box on;
hold off;


%% Biomass composition
% protein (%)
subplot(a, b, find(strcmp(figure_output,'protein')))
hold on;
if exist('protein_clim','var')
papers = fieldnames(protein_clim); 
for i = 1: length(papers)
plot(D_clim.(papers{i}),  protein_clim.(papers{i}) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end 
if exist('protein_nlim','var')
papers = fieldnames(protein_nlim); 
for i = 1: length(papers)
plot(D_nlim.(papers{i}),  protein_nlim.(papers{i}) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end
%_% plot(D, 100*(total_protein_con_clim*110/10^6)./par.rho_cell, '-', 'Color', color_clim.default);
%_% plot(D, 100*(total_protein_con_nlim*110/10^6)./par.rho_cell, '-', 'Color', color_nlim.default);  
%_% plot(D, 100*(total_protein_con_cnlim1*110/10^6)./par.rho_cell, '-', 'Color', color_cnlim1.default);
%_% plot(D, 100*(total_protein_con_cnlim2*110/10^6)./par.rho_cell, '-', 'Color', color_cnlim2.default);  
plot(D, 100.*(total_protein_con_clim.*pt_uM_to_gPerL)./par.rho_cell, '-', 'Color', color_clim.default, 'LineWidth',2);
plot(D, 100.*(total_protein_con_nlim.*pt_uM_to_gPerL)./par.rho_cell, '-', 'Color', color_nlim.default, 'LineWidth',2); 
plot(D, 100.*(total_protein_con_cnlim1.*pt_uM_to_gPerL)./par.rho_cell, '-', 'Color', color_cnlim1.default, 'LineWidth',2);
plot(D, 100.*(total_protein_con_cnlim2.*pt_uM_to_gPerL)./par.rho_cell, '-', 'Color', color_cnlim2.default, 'LineWidth',2); 

hold off 
ylabel('protein (%)')
xlabel(x_label);
xlim(x_lim);
ylim(ylim_prot)
yticks([0 20 40])%_%
%_% patch_background(x_left, x_middle, max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), min(ylim)+labely*(max(ylim)-min(ylim)), figure_id(strcmp(figure_output,'Eat')), 'fontweight', 'bold', 'FontSize', fontsize);
axis square; 
box on;
hold off;

% lipid (%) 
subplot(a, b, find(strcmp(figure_output,'lipid')))
hold on;
if exist('lipid_clim','var')
papers = fieldnames(lipid_clim); 
for i = 1: length(papers)
plot(D_clim.(papers{i}),  lipid_clim.(papers{i}) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end
if exist('lipid_nlim','var')
papers = fieldnames(lipid_nlim); 
for i = 1: length(papers)
plot(D_nlim.(papers{i}),  lipid_nlim.(papers{i})  ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end     
end 
% add small number to avoid dividing 0 
%_% plot(D,   100.*WT_y_steady_clim(:,num_y.lp_e).*lp_uM_to_gPerL./(WT_y_steady_clim(:,num_y.cells)+10^-10)+lipid_constant_part , 'Color', color_clim.default);
%_% plot(D,   100.*WT_y_steady_nlim(:,num_y.lp_e).*lp_uM_to_gPerL./(WT_y_steady_nlim(:,num_y.cells)+10^-10)+lipid_constant_part , 'Color', color_nlim.default); 
%_% plot(D,   100.*WT_y_steady_cnlim1(:,num_y.lp_e).*lp_uM_to_gPerL./(WT_y_steady_cnlim1(:,num_y.cells)+10^-10)+lipid_constant_part , 'Color', color_cnlim1.default);
%_% plot(D,   100.*WT_y_steady_cnlim2(:,num_y.lp_e).*lp_uM_to_gPerL./(WT_y_steady_cnlim2(:,num_y.cells)+10^-10)+lipid_constant_part , 'Color', color_cnlim2.default);
plot(D,   100.*WT_y_steady_clim(:,num_y.lp_e).*lp_uM_to_gPerL./par.rho_cell +lipid_constant_part , 'Color', color_clim.default, 'LineWidth',2);
plot(D,   100.*WT_y_steady_nlim(:,num_y.lp_e).*lp_uM_to_gPerL./par.rho_cell +lipid_constant_part , 'Color', color_nlim.default, 'LineWidth',2); 
plot(D,   100.*WT_y_steady_cnlim1(:,num_y.lp_e).*lp_uM_to_gPerL./par.rho_cell +lipid_constant_part , 'Color', color_cnlim1.default, 'LineWidth',2);
plot(D,   100.*WT_y_steady_cnlim2(:,num_y.lp_e).*lp_uM_to_gPerL./par.rho_cell +lipid_constant_part , 'Color', color_cnlim2.default, 'LineWidth',2);

hold off 
ylabel('lipid (%)')
xlabel(x_label);
xlim(x_lim);
ylim(ylim_lipid)
yticks([0 20 40])%_%
%_% patch_background(x_left, x_middle, max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), min(ylim)+labely*(max(ylim)-min(ylim)), figure_id(strcmp(figure_output,'Eat')), 'fontweight', 'bold', 'FontSize', fontsize);
axis square; 
box on;
hold off;
set(gcf,'Position',position) % tight    



%% legends

figure
subplot(2,1,1)
papers = fieldnames(paperinfo_clim.CN); 
hold on 
for i = 1:length(papers)
plot(0,  0, shape_clim.(papers{i}) ,...
'Color', color_clim.(papers{i}), 'MarkerFaceColor', color_clim.(papers{i}))
end
hold off

%papers = fieldnames(D_Nlim);
% ASCII character, 32 = space, 44 = comma
for i = 1:length(papers)
legs{i} =  strcat(paperinfo_clim.author.(papers{i}), ', ', 32, ...
paperinfo_clim.strain.(papers{i}), ...
', C/N=' , ...
num2str(paperinfo_clim.CN.(papers{i}),'%.1f'),...
', Glu=',...
num2str(paperinfo_clim.gl_in.(papers{i}),'%.1f'),...
'g/L', ...
', (NH4)2SO4=',...
num2str(paperinfo_clim.n_in.(papers{i}),'%.1f'),...
'g/L');
end

legend(legs, Interpreter="none")
legend box off


legs = {}; 
subplot(2,1,2)
papers = fieldnames(paperinfo_nlim.CN); 
hold on 
for i = 1:length(papers)
plot(0,  0, shape_nlim.(papers{i}) ,...
'Color', color_nlim.(papers{i}), 'MarkerFaceColor', color_nlim.(papers{i}))
end
hold off

%papers = fieldnames(D_Nlim);
% ASCII character, 32 = space, 44 = comma
for i = 1:length(papers)
legs{i} =  strcat(paperinfo_nlim.author.(papers{i}), ', ', 32, ...
paperinfo_nlim.strain.(papers{i}), ...
', C/N=' , ...
num2str(paperinfo_nlim.CN.(papers{i}),'%.1f'),...
', Glu=',...
num2str(paperinfo_nlim.gl_in.(papers{i}),'%.1f'),...
'g/L', ...
', (NH4)2SO4=',...
num2str(paperinfo_nlim.n_in.(papers{i}),'%.1f'),...
'g/L');
end

legend(legs, Interpreter="none")
legend box off

end
