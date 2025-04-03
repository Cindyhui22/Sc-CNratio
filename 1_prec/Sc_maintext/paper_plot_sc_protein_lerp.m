function paper_plot_sc_protein_lerp(D, D_crit_simu, ...
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
color_clim.hackett    = plt_clrs.green;

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

%% proteins 

figure; 
% gy
subplot(a, b, find(strcmp(figure_output,'Egy')))
hold on;
if exist('prot_sec_clim','var') 
papers = fieldnames(prot_sec_clim.gy); 
for i = 1: length(papers)
plot(D_clim.(papers{i}),  log2(prot_sec_clim.gy.(papers{i})./prot_sec_clim.gy.(papers{i})(1)) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end
if exist('prot_sec_nlim','var') 
papers = fieldnames(prot_sec_nlim.gy); 
for i = 1: length(papers)
plot(D_nlim.(papers{i}),  log2(prot_sec_nlim.gy.(papers{i})./prot_sec_nlim.gy.(papers{i})(1)) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end  
end
%plot(D, WT_y_steady_1gL(:,num_y.gy),  '-', 'Color', color_clim.kumar_lg);
plot(D, WT_y_steady_clim(:,num_y.gy), '-', 'Color', color_clim.default, 'LineWidth',2);
plot(D, WT_y_steady_nlim(:,num_y.gy), '-', 'Color', color_nlim.default, 'LineWidth',2);

ylabel(figcabbi.y_label.gy,'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
ylim(y_lim_prot);
yticks([-3 -1 1])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), min(ylim)+labely*(max(ylim)-min(ylim)), figure_id(strcmp(figure_output,'Egy')), 'fontweight', 'bold', 'FontSize', fontsize);
axis square; 
box on;
hold off;

% fe 
subplot(a, b, find(strcmp(figure_output,'Efe')))
hold on;
if exist('prot_sec_clim','var')
papers = fieldnames(prot_sec_clim.fe); 
for i = 1: length(papers)
plot(D_clim.(papers{i}),  log2(prot_sec_clim.fe.(papers{i})./prot_sec_clim.fe.(papers{i})(1)) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end
if exist('prot_sec_nlim','var')
papers = fieldnames(prot_sec_nlim.fe); 
for i = 1: length(papers)
plot(D_nlim.(papers{i}),  log2(prot_sec_nlim.fe.(papers{i})./prot_sec_nlim.fe.(papers{i})(1)) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end  
end 
%plot(D, WT_y_steady_1gL(:,num_y.fe),  '-', 'Color', color_clim.kumar_lg);
plot(D, WT_y_steady_clim(:,num_y.fe), '-', 'Color', color_clim.default, 'LineWidth',2);
plot(D, WT_y_steady_nlim(:,num_y.fe), '-', 'Color', color_nlim.default, 'LineWidth',2);

ylabel(figcabbi.y_label.fe,'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
ylim(y_lim_prot);
yticks([-3 -1 1])
%_% patch_background(x_left, x_middle, max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), min(ylim)+labely*(max(ylim)-min(ylim)), figure_id(strcmp(figure_output,'Efe')), 'fontweight', 'bold', 'FontSize', fontsize);
axis square; 
box on;
hold off;

% gn 
subplot(a, b, find(strcmp(figure_output,'Egn')))
hold on;
if exist('prot_sec_clim','var')
papers = fieldnames(prot_sec_clim.gn); 
for i = 1: length(papers)
plot(D_clim.(papers{i}),  log2(prot_sec_clim.gn.(papers{i})./prot_sec_clim.gn.(papers{i})(1)) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end 
if exist('prot_sec_nlim','var')
papers = fieldnames(prot_sec_nlim.gn); 
for i = 1: length(papers)
plot(D_nlim.(papers{i}),  log2(prot_sec_nlim.gn.(papers{i})./prot_sec_nlim.gn.(papers{i})(1)) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end   
end 
%plot(D, WT_y_steady_1gL(:,num_y.gn),  '-', 'Color', color_clim.kumar_lg);
plot(D, WT_y_steady_clim(:,num_y.gn), '-', 'Color', color_clim.default, 'LineWidth',2);
plot(D, WT_y_steady_nlim(:,num_y.gn), '-', 'Color', color_nlim.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim1(:,num_y.gn), '-', 'Color', color_cnlim1.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim2(:,num_y.gn), '-', 'Color', color_cnlim2.default, 'LineWidth',2);

ylabel(figcabbi.y_label.gn,'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
ylim(y_lim_prot);
yticks([-3 -1 1])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), min(ylim)+labely*(max(ylim)-min(ylim)), figure_id(strcmp(figure_output,'Egn')),'fontweight','bold','FontSize',fontsize);
axis square; 
box on;
hold off;

% mt 
subplot(a, b, find(strcmp(figure_output,'Emt')))
hold on;
if exist('prot_sec_clim','var')
papers = fieldnames(prot_sec_clim.mt); 
for i = 1: length(papers)
plot(D_clim.(papers{i}),  log2(prot_sec_clim.mt.(papers{i})./prot_sec_clim.mt.(papers{i})(1)) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end 
if exist('prot_sec_nlim','var')
papers = fieldnames(prot_sec_nlim.mt); 
for i = 1: length(papers)
plot(D_nlim.(papers{i}),  log2(prot_sec_nlim.mt.(papers{i})./prot_sec_nlim.mt.(papers{i})(1))  ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end
%plot(D, WT_y_steady_1gL(:,num_y.mt),  '-', 'Color', color_clim.kumar_lg);
plot(D, WT_y_steady_clim(:,num_y.mt), '-', 'Color', color_clim.default, 'LineWidth',2);
plot(D, WT_y_steady_nlim(:,num_y.mt), '-', 'Color', color_nlim.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim1(:,num_y.mt), '-', 'Color', color_cnlim1.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim2(:,num_y.mt), '-', 'Color', color_cnlim2.default, 'LineWidth',2);


ylabel(figcabbi.y_label.mt,'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
ylim(y_lim_prot);
yticks([-3 -1 1])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), min(ylim)+labely*(max(ylim)-min(ylim)), figure_id(strcmp(figure_output,'Emt')),'fontweight','bold','FontSize',fontsize);
axis square; 
box on;
hold off;

% as 
subplot(a, b, find(strcmp(figure_output,'Eas')))
hold on;
if  exist('prot_sec_clim','var')
papers = fieldnames(prot_sec_clim.as); 
for i = 1: length(papers)
plot(D_clim.(papers{i}),  log2(prot_sec_clim.as.(papers{i})./prot_sec_clim.as.(papers{i})(1)) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end 
if exist('prot_sec_nlim','var')
papers = fieldnames(prot_sec_nlim.as); 
for i = 1: length(papers)
plot(D_nlim.(papers{i}),  log2(prot_sec_nlim.as.(papers{i})./prot_sec_nlim.as.(papers{i})(1)) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end 
%plot(D, WT_y_steady_1gL(:,num_y.as),  '-', 'Color', color_clim.kumar_lg);
plot(D, WT_y_steady_clim(:,num_y.as), '-', 'Color', color_clim.default, 'LineWidth',2);
plot(D, WT_y_steady_nlim(:,num_y.as), '-', 'Color', color_nlim.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim1(:,num_y.as), '-', 'Color', color_cnlim1.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim2(:,num_y.as), '-', 'Color', color_cnlim2.default, 'LineWidth',2);

ylabel(figcabbi.y_label.as,'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
ylim(y_lim_prot);
yticks([-3 -1 1])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), min(ylim)+labely*(max(ylim)-min(ylim)), figure_id(strcmp(figure_output,'Eas')), 'fontweight', 'bold', 'FontSize', fontsize);
axis square; 
box on;
hold off;

% at 
%
subplot(a, b, find(strcmp(figure_output,'Eat')))
hold on;
if exist('prot_sec_clim','var')
papers = fieldnames(prot_sec_clim.at); 
for i = 1: length(papers)
plot(D_clim.(papers{i}),  log2(prot_sec_clim.at.(papers{i})./prot_sec_clim.at.(papers{i})(1)) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end
if exist('prot_sec_nlim','var')
papers = fieldnames(prot_sec_nlim.at); 
for i = 1: length(papers)
plot(D_nlim.(papers{i}),  log2(prot_sec_nlim.at.(papers{i})./prot_sec_nlim.at.(papers{i})(1)) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end
%plot(D, WT_y_steady_1gL(:,num_y.at),  '-', 'Color', color_clim.kumar_lg);
plot(D, WT_y_steady_clim(:,num_y.at), '-', 'Color', color_clim.default, 'LineWidth',2);
plot(D, WT_y_steady_nlim(:,num_y.at), '-', 'Color', color_nlim.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim1(:,num_y.at), '-', 'Color', color_cnlim1.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim2(:,num_y.at), '-', 'Color', color_cnlim2.default, 'LineWidth',2);

ylabel(figcabbi.y_label.at,'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
ylim(y_lim_prot);
yticks([-3 -1 1])
%_% patch_background(x_left, x_middle, max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), min(ylim)+labely*(max(ylim)-min(ylim)), figure_id(strcmp(figure_output,'Eat')), 'fontweight', 'bold', 'FontSize', fontsize);
axis square; 
box on;
hold off;
%}

% lp 
subplot(a, b, find(strcmp(figure_output,'Elp')))
hold on;
if exist('prot_sec_clim','var')
papers = fieldnames(prot_sec_clim.lp); 
for i = 1: length(papers)
plot(D_clim.(papers{i}),  log2(prot_sec_clim.lp.(papers{i})./prot_sec_clim.lp.(papers{i})(1)),...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end
if exist('prot_sec_nlim','var')
papers = fieldnames(prot_sec_nlim.lp); 
for i = 1: length(papers)
plot(D_nlim.(papers{i}),  log2(prot_sec_nlim.lp.(papers{i})./prot_sec_nlim.lp.(papers{i})(1)) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end    
end 
%plot(D, WT_y_steady_1gL(:,num_y.as),  '-', 'Color', color_clim.kumar_lg);
plot(D, WT_y_steady_clim(:,num_y.lp), '-', 'Color', color_clim.default, 'LineWidth',2);
plot(D, WT_y_steady_nlim(:,num_y.lp), '-', 'Color', color_nlim.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim1(:,num_y.lp), '-', 'Color', color_cnlim1.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim2(:,num_y.lp), '-', 'Color', color_cnlim2.default, 'LineWidth',2);


ylabel(figcabbi.y_label.lp,'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
ylim(y_lim_prot);
yticks([-3 -1 1])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), min(ylim)+labely*(max(ylim)-min(ylim)), figure_id(strcmp(figure_output,'Eas')), 'fontweight', 'bold', 'FontSize', fontsize);
axis square; 
box on;
hold off;

% Elo
%
subplot(a, b, find(strcmp(figure_output,'Elo')))
hold on;
if exist('prot_sec_clim','var')
papers = fieldnames(prot_sec_clim.lo); 
for i = 1: length(papers)
plot(D_clim.(papers{i}),  log2(prot_sec_clim.lo.(papers{i})./prot_sec_clim.lo.(papers{i})(1)),...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end
if exist('prot_sec_nlim','var')
papers = fieldnames(prot_sec_nlim.lo); 
for i = 1: length(papers)
plot(D_nlim.(papers{i}),  log2(prot_sec_nlim.lo.(papers{i})./prot_sec_nlim.lo.(papers{i})(1)) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end    
end 
plot(D, WT_y_steady_clim(:,num_y.lo), '-', 'Color', color_clim.default, 'LineWidth',2);
plot(D, WT_y_steady_nlim(:,num_y.lo), '-', 'Color', color_nlim.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim1(:,num_y.lo), '-', 'Color', color_cnlim1.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim2(:,num_y.lo), '-', 'Color', color_cnlim2.default, 'LineWidth',2);

ylabel(figcabbi.y_label.lo,'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
ylim(y_lim_prot);
yticks([-3 -1 1])
%_% patch_background(x_left, x_middle, max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), min(ylim)+labely*(max(ylim)-min(ylim)), figure_id(strcmp(figure_output,'Eat')), 'fontweight', 'bold', 'FontSize', fontsize);
axis square; 
box on;
hold off;

% sp
%{
subplot(a, b, find(strcmp(figure_output,'Esp')))
hold on;
if exist('prot_sec_clim','var')
papers = fieldnames(prot_sec_clim.sp); 
for i = 1: length(papers)
plot(D_clim.(papers{i}),  log2(prot_sec_clim.sp.(papers{i})./prot_sec_clim.sp.(papers{i})(1)) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end
if exist('prot_sec_nlim','var')
papers = fieldnames(prot_sec_nlim.sp); 
for i = 1: length(papers)
plot(D_nlim.(papers{i}),  log2(prot_sec_nlim.sp.(papers{i})./prot_sec_nlim.sp.(papers{i})(1)) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end 
plot(D, WT_y_steady_clim(:,num_y.sp), '-', 'Color', color_clim.default);
plot(D, WT_y_steady_nlim(:,num_y.sp), '-', 'Color', color_nlim.default);
plot(D, WT_y_steady_cnlim1(:,num_y.sp), '-', 'Color', color_cnlim1.default);
plot(D, WT_y_steady_cnlim2(:,num_y.sp), '-', 'Color', color_cnlim2.default);

ylabel(figcabbi.y_label.sp,'Color',xy_label_color)
xlabel(x_label);
xlim(x_lim);
ylim(y_lim_prot);
yticks([-3 -1 1])
patch_background(x_left, x_middle, max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), min(ylim)+labely*(max(ylim)-min(ylim)), figure_id(strcmp(figure_output,'Eat')), 'fontweight', 'bold', 'FontSize', fontsize);
axis square; 
box on;
hold off;
%}

% sd
%{
subplot(a, b, find(strcmp(figure_output,'Esd')))
hold on;
if exist('prot_sec_clim','var')
papers = fieldnames(prot_sec_clim.sd); 
for i = 1: length(papers)
plot(D_clim.(papers{i}),  log2(prot_sec_clim.sd.(papers{i})./prot_sec_clim.sd.(papers{i})(1)) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end 
if exist('prot_sec_nlim','var')
papers = fieldnames(prot_sec_nlim.sd); 
for i = 1: length(papers)
plot(D_nlim.(papers{i}),  log2(prot_sec_nlim.sd.(papers{i})./prot_sec_nlim.sd.(papers{i})(1)) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end 
plot(D, WT_y_steady_clim(:,num_y.sd), '-', 'Color', color_clim.default);
plot(D, WT_y_steady_nlim(:,num_y.sd), '-', 'Color', color_nlim.default);
plot(D, WT_y_steady_cnlim1(:,num_y.sd), '-', 'Color', color_cnlim1.default);
plot(D, WT_y_steady_cnlim2(:,num_y.sd), '-', 'Color', color_cnlim2.default);

ylabel(figcabbi.y_label.sd,'Color',xy_label_color)
xlabel(x_label);
xlim(x_lim);
ylim(y_lim_prot);
yticks([-3 -1 1])
patch_background(x_left, x_middle, max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), min(ylim)+labely*(max(ylim)-min(ylim)), figure_id(strcmp(figure_output,'Eat')), 'fontweight', 'bold', 'FontSize', fontsize);
axis square; 
box on;
hold off;
%}


% R
%
subplot(a, b, find(strcmp(figure_output,'R')))
hold on;
if exist('prot_sec_clim','var')
papers = fieldnames(prot_sec_clim.r); 
for i = 1: length(papers)
plot(D_clim.(papers{i}),  log2(prot_sec_clim.r.(papers{i})./prot_sec_clim.r.(papers{i})(1)) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end
if exist('prot_sec_nlim','var')
papers = fieldnames(prot_sec_nlim.r); 
for i = 1: length(papers)
plot(D_nlim.(papers{i}),  log2(prot_sec_nlim.r.(papers{i})./prot_sec_nlim.r.(papers{i})(1)) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end 
plot(D, WT_y_steady_clim(:,num_y.r), '-', 'Color', color_clim.default, 'LineWidth',2);
plot(D, WT_y_steady_nlim(:,num_y.r), '-', 'Color', color_nlim.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim1(:,num_y.r), '-', 'Color', color_cnlim1.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim2(:,num_y.r), '-', 'Color', color_cnlim2.default, 'LineWidth',2);

ylabel(figcabbi.y_label.r,'Color',xy_label_color)
xlabel(x_label);
xlim(x_lim);
ylim(y_lim_prot);
yticks([-3 -1 1])
%_% patch_background(x_left, x_middle, max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), min(ylim)+labely*(max(ylim)-min(ylim)), figure_id(strcmp(figure_output,'Eat')), 'fontweight', 'bold', 'FontSize', fontsize);
axis square; 
box on;
hold off;


% Z
%
subplot(a, b, find(strcmp(figure_output,'Z')))
hold on;
if exist('prot_sec_clim','var')
papers = fieldnames(prot_sec_clim.z); 
for i = 1: length(papers)
plot(D_clim.(papers{i}),  log2(prot_sec_clim.z.(papers{i})./prot_sec_clim.z.(papers{i})(1)) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end
if exist('prot_sec_nlim','var')
papers = fieldnames(prot_sec_nlim.z); 
for i = 1: length(papers)
plot(D_nlim.(papers{i}),  log2(prot_sec_nlim.z.(papers{i})./prot_sec_nlim.z.(papers{i})(1)) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end 
plot(D, WT_y_steady_clim(:,num_y.z), '-', 'Color', color_clim.default, 'LineWidth',2);
plot(D, WT_y_steady_nlim(:,num_y.z), '-', 'Color', color_nlim.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim1(:,num_y.z), '-', 'Color', color_cnlim1.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim2(:,num_y.z), '-', 'Color', color_cnlim2.default, 'LineWidth',2);

ylabel(figcabbi.y_label.z,'Color',xy_label_color)
xlabel(x_label);
xlim(x_lim);
ylim(y_lim_prot);
yticks([-3 -1 1])
%_% patch_background(x_left, x_middle, max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), min(ylim)+labely*(max(ylim)-min(ylim)), figure_id(strcmp(figure_output,'Eat')), 'fontweight', 'bold', 'FontSize', fontsize);
axis square; 
box on;
hold off;

set(gcf,'Position',position)



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

protein_colors = {'#ABD037';
'#0FA64A';
'#0778A9';
'#3D57A7';
'#6E4EA0';
'#773A96';
'#A92179';
'#EE3324';
'#F47C20';
'#FBA919';
'#FED304';};

index = 12;
figure; 
subplot(1, 3, 1)
hold on;
%plot(hackett.dr, hackett.lo, 'o', 'MarkerFaceColor', color_sc.hackett)
plot(100*WT_y_steady_clim_frac(index,num_y.r), 100*WT_y_steady_nlim_frac(index,num_y.r) , 'o', 'Color', protein_colors{11}, 'MarkerSize', 8, 'LineWidth', 2);
plot(100*WT_y_steady_clim_frac(index,num_y.z), 100*WT_y_steady_nlim_frac(index,num_y.z) , 'o', 'Color', protein_colors{2}, 'MarkerSize', 8, 'LineWidth', 2);
plot(100*WT_y_steady_clim_frac(index,num_y.gy), 100*WT_y_steady_nlim_frac(index,num_y.gy) , 'o', 'Color', protein_colors{3}, 'MarkerSize', 8, 'LineWidth', 2);
plot(100*WT_y_steady_clim_frac(index,num_y.fe), 100*WT_y_steady_nlim_frac(index,num_y.fe) , 'o', 'Color', protein_colors{4}, 'MarkerSize', 8, 'LineWidth', 2);

plot(100*WT_y_steady_clim_frac(index,num_y.gn), 100*WT_y_steady_nlim_frac(index,num_y.gn) , 'o', 'Color', protein_colors{5}, 'MarkerSize', 8, 'LineWidth', 2);
plot(100*WT_y_steady_clim_frac(index,num_y.mt), 100*WT_y_steady_nlim_frac(index,num_y.mt) , 'o', 'Color', protein_colors{6}, 'MarkerSize', 8, 'LineWidth', 2);
plot(100*WT_y_steady_clim_frac(index,num_y.as), 100*WT_y_steady_nlim_frac(index,num_y.as) , 'o', 'Color', protein_colors{7}, 'MarkerSize', 8, 'LineWidth', 2);
plot(100*WT_y_steady_clim_frac(index,num_y.at), 100*WT_y_steady_nlim_frac(index,num_y.at) , 'o', 'Color', protein_colors{8}, 'MarkerSize', 8, 'LineWidth', 2);
plot(100*WT_y_steady_clim_frac(index,num_y.lp), 100*WT_y_steady_nlim_frac(index,num_y.lp) , 'o', 'Color', protein_colors{9}, 'MarkerSize', 8, 'LineWidth', 2);
plot(100*WT_y_steady_clim_frac(index,num_y.lo), 100*WT_y_steady_nlim_frac(index,num_y.lo) , 'o', 'Color', protein_colors{10}, 'MarkerSize', 8, 'LineWidth', 2);

plot([0.01 100], [0.01 100], 'k-', 'LineWidth', 0.5)
plot([0.01 50], [0.02 100], 'k:', 'LineWidth', 0.5)
plot([0.01 200], [0.005 100], 'k:', 'LineWidth', 0.5)
xlabel(x_label);
xlim([0.01 100])
ylim([0.01 100])
xticks([0.01 1 100])
yticks([0.01 1 100])
set(gca, 'XScale','log')
set(gca, 'YScale','log')
set(gca, 'XMinorTick', 'off')
set(gca, 'YMinorTick', 'off')
xlabel('fraction clim')
ylabel('fraction nlim')
axis square; 
box on;
hold off;
%legend('R', 'Z', 'E_{gy}', 'E_{fe}', 'E_{gn}', 'E_{mt}', 'E_{as}', 'E_{lp}', 'E_{lo}')

if strcmp(yeast_type,'sc')
subplot(1, 3, 2)
hold on;
%plot(hackett.dr, hackett.lo, 'o', 'MarkerFaceColor', color_sc.hackett)
plot(100*prot_sec_clim.r.Yu4, 100*prot_sec_nlim.r.Yu50 , 'o', 'Color', protein_colors{11}, 'MarkerSize', 8, 'LineWidth', 2);
plot(100*prot_sec_clim.z.Yu4, 100*prot_sec_nlim.z.Yu50, 'o', 'Color', protein_colors{2}, 'MarkerSize', 8, 'LineWidth', 2);
plot(100*prot_sec_clim.gy.Yu4, 100*prot_sec_nlim.gy.Yu50 , 'o', 'Color', protein_colors{3}, 'MarkerSize', 8, 'LineWidth', 2);
plot(100*prot_sec_clim.fe.Yu4, 100*prot_sec_nlim.fe.Yu50 , 'o', 'Color', protein_colors{4}, 'MarkerSize', 8, 'LineWidth', 2);

plot(100*prot_sec_clim.gn.Yu4, 100*prot_sec_nlim.gn.Yu50 , 'o', 'Color', protein_colors{5}, 'MarkerSize', 8, 'LineWidth', 2);
plot(100*prot_sec_clim.mt.Yu4, 100*prot_sec_nlim.mt.Yu50 , 'o', 'Color', protein_colors{6}, 'MarkerSize', 8, 'LineWidth', 2);
plot(100*prot_sec_clim.as.Yu4, 100*prot_sec_nlim.as.Yu50 , 'o', 'Color', protein_colors{7}, 'MarkerSize', 8, 'LineWidth', 2);
plot(100*prot_sec_clim.at.Yu4, 100*prot_sec_nlim.at.Yu50 , 'o', 'Color', protein_colors{8}, 'MarkerSize', 8, 'LineWidth', 2);
plot(100*prot_sec_clim.lp.Yu4, 100*prot_sec_nlim.lp.Yu50 , 'o', 'Color', protein_colors{9}, 'MarkerSize', 8, 'LineWidth', 2);
plot(100*prot_sec_clim.lo.Yu4, 100*prot_sec_nlim.lo.Yu50 , 'o', 'Color', protein_colors{10}, 'MarkerSize', 8, 'LineWidth', 2);

plot([0.01 100], [0.01 100], 'k-', 'LineWidth', 0.5)
plot([0.01 50], [0.02 100], 'k:', 'LineWidth', 0.5)
plot([0.01 200], [0.005 100], 'k:', 'LineWidth', 0.5)
xlabel(x_label);
xlim([0.01 100])
ylim([0.01 100])
xticks([0.01 1 100])
yticks([0.01 1 100])
set(gca, 'XScale','log')
set(gca, 'YScale','log')
set(gca, 'XMinorTick', 'off')
set(gca, 'YMinorTick', 'off')
xlabel('fraction clim')
ylabel('fraction nlim')
axis square; 
box on;
hold off;
legend('R', 'Z', 'E_{gy}', 'E_{fe}', 'E_{gn}', 'E_{mt}', 'E_{as}', 'E_{lp}', 'E_{lo}')
%_% legend('R', 'Z', 'E_{gy}', 'E_{fe}', 'E_{gn}', 'E_{mt}', 'E_{as}', 'E_{at}', 'E_{lp}', 'E_{lo}')%_%
end 

if strcmp(yeast_type,'rt')
subplot(1, 3, 2)
hold on;
%plot(hackett.dr, hackett.lo, 'o', 'MarkerFaceColor', color_sc.hackett)
plot(100*mRNA_sec_clim.r.Zhu, 100*mRNA_sec_nlim.r.Zhu , 'o', 'Color', protein_colors{11}, 'MarkerSize', 8, 'LineWidth', 2);
plot(100*mRNA_sec_clim.z.Zhu, 100*mRNA_sec_nlim.z.Zhu, 'o', 'Color', protein_colors{2}, 'MarkerSize', 8, 'LineWidth', 2);
plot(100*mRNA_sec_clim.gy.Zhu, 100*mRNA_sec_nlim.gy.Zhu , 'o', 'Color', protein_colors{3}, 'MarkerSize', 8, 'LineWidth', 2);
plot(100*mRNA_sec_clim.fe.Zhu, 100*mRNA_sec_nlim.fe.Zhu , 'o', 'Color', protein_colors{4}, 'MarkerSize', 8, 'LineWidth', 2);

plot(100*mRNA_sec_clim.gn.Zhu, 100*mRNA_sec_nlim.gn.Zhu , 'o', 'Color', protein_colors{5}, 'MarkerSize', 8, 'LineWidth', 2);
plot(100*mRNA_sec_clim.mt.Zhu, 100*mRNA_sec_nlim.mt.Zhu , 'o', 'Color', protein_colors{6}, 'MarkerSize', 8, 'LineWidth', 2);
plot(100*mRNA_sec_clim.as.Zhu, 100*mRNA_sec_nlim.as.Zhu , 'o', 'Color', protein_colors{7}, 'MarkerSize', 8, 'LineWidth', 2);
%plot(100*prot_sec_clim.at.Yu4, 100*prot_sec_nlim.at.Yu50 , 'o', 'Color', protein_colors{8}, 'MarkerSize', 8, 'LineWidth', 2);
plot(100*mRNA_sec_clim.lp.Zhu, 100*mRNA_sec_nlim.lp.Zhu , 'o', 'Color', protein_colors{9}, 'MarkerSize', 8, 'LineWidth', 2);
plot(100*mRNA_sec_clim.lo.Zhu, 100*mRNA_sec_nlim.lo.Zhu , 'o', 'Color', protein_colors{10}, 'MarkerSize', 8, 'LineWidth', 2);

plot([0.01 100], [0.01 100], 'k-', 'LineWidth', 0.5)
plot([0.01 50], [0.02 100], 'k:', 'LineWidth', 0.5)
plot([0.01 200], [0.005 100], 'k:', 'LineWidth', 0.5)
xlabel(x_label);
xlim([0.01 100])
ylim([0.01 100])
xticks([0.01 1 100])
yticks([0.01 1 100])
set(gca, 'XScale','log')
set(gca, 'YScale','log')
set(gca, 'XMinorTick', 'off')
set(gca, 'YMinorTick', 'off')
xlabel('fraction clim')
ylabel('fraction nlim')
axis square; 
box on;
hold off;
legend('R', 'Z', 'E_{gy}', 'E_{fe}', 'E_{gn}', 'E_{mt}', 'E_{as}', 'E_{lp}', 'E_{lo}')
%_% legend('R', 'Z', 'E_{gy}', 'E_{fe}', 'E_{gn}', 'E_{mt}', 'E_{as}', 'E_{at}', 'E_{lp}', 'E_{lo}')%_%

end 


set(gcf,'Position',[440 684 596 114])








% if exist('prot_sec_nlim','var')
%     papers = fieldnames(prot_sec_nlim.r); 
%     for i = 1: length(papers)
%         plot(D_nlim.(papers{i}),  log2(prot_sec_nlim.r.(papers{i})./prot_sec_nlim.r.(papers{i})(1)) ,...
%             shape_nlim.(papers{i}),...
%             'Color', color_nlim.(papers{i}), ...
%             'MarkerFaceColor', color_nlim.(papers{i}))
%     end 
% end 
