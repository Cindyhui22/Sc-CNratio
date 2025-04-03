function paper_plot_sc_SI_lerp(D, D_crit_simu, ...
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



%% plotting stuff
%biomass_log_climale = 1;
paper_figure_color; 
paper_figure_label; 


figure_output = figure_output_format;

[n, m]   = size(figure_output); 
position = figure_position;

x_label  = 'dilution rate (h^{-1})'; 


x_lim      = [0  0.42]; 


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

%% plot before normalizing  ----------------------------------------------
% aa, atp, pc, proteins 

figure 
subplot(m, n, find(strcmp(figure_output,'precursor')))
hold on
plot(D, WT_y_steady_clim(:,num_y.pc), '-', 'Color', color_clim.default, 'LineWidth',2);
plot(D, WT_y_steady_nlim(:,num_y.pc), '-', 'Color', color_nlim.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim1(:,num_y.pc), '-', 'Color', color_cnlim1.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim2(:,num_y.pc), '-', 'Color', color_cnlim2.default, 'LineWidth',2);

ylabel(figcabbi.y_label.pc,'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
ylim([0 max(ylim)*1.05])

%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'precursor')),'fontweight','bold','FontSize',fontsize);
axis square; 
box on;
hold off;


% subplot(m, n, find(strcmp(figure_output,'atp')))
% 
% hold on
% plot(D, WT_y_steady_clim(:,num_y.ae), '-', 'Color', color_clim.default);
% plot(D, WT_y_steady_nlim(:,num_y.ae), '-', 'Color', color_nlim.default);
% plot(D, WT_y_steady_cnlim1(:,num_y.ae), '-', 'Color', color_cnlim1.default);
% plot(D, WT_y_steady_cnlim2(:,num_y.ae), '-', 'Color', color_cnlim2.default);
% 
% ylabel(figcabbi.y_label.ae,'Color',xy_label_color)
% xlabel(x_label,'Color',xy_label_color);
% xlim(x_lim);
% ylim([0 max(ylim)*1.05])
% 
% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
% %text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'precursor')),'fontweight','bold','FontSize',fontsize);
% axis square; 
% box on;
% hold off;



subplot(m, n, find(strcmp(figure_output,'aa_in')))

hold on
plot(D, WT_y_steady_clim(:,num_y.aa_in), '-', 'Color', color_clim.default, 'LineWidth',2);
plot(D, WT_y_steady_nlim(:,num_y.aa_in), '-', 'Color', color_nlim.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim1(:,num_y.aa_in), '-', 'Color', color_cnlim1.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim2(:,num_y.aa_in), '-', 'Color', color_cnlim2.default, 'LineWidth',2);

ylabel(figcabbi.y_label.aa,'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
ylim([0 max(ylim)*1.05])

%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'precursor')),'fontweight','bold','FontSize',fontsize);
axis square; 
box on;
hold off;


subplot(m, n, find(strcmp(figure_output,'Egn')))
hold on
plot(D, WT_y_steady_clim_frac(:,num_y.gn), '-', 'Color', color_clim.default, 'LineWidth',2);
plot(D, WT_y_steady_nlim_frac(:,num_y.gn), '-', 'Color', color_nlim.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim1_frac(:,num_y.gn), '-', 'Color', color_cnlim1.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim2_frac(:,num_y.gn), '-', 'Color', color_cnlim2.default, 'LineWidth',2);

ylabel(figcabbi.y_label.gn,'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
ylim([0 max(ylim)*1.05])

%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'precursor')),'fontweight','bold','FontSize',fontsize);
axis square; 
box on;
hold off;



subplot(m, n, find(strcmp(figure_output,'R')))
hold on
plot(D, WT_y_steady_clim_frac(:,num_y.r), '-', 'Color', color_clim.default, 'LineWidth',2);
plot(D, WT_y_steady_nlim_frac(:,num_y.r), '-', 'Color', color_nlim.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim1_frac(:,num_y.r), '-', 'Color', color_cnlim1.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim2_frac(:,num_y.r), '-', 'Color', color_cnlim2.default, 'LineWidth',2);

ylabel(figcabbi.y_label.r,'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
ylim([0 max(ylim)*1.05])

%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'precursor')),'fontweight','bold','FontSize',fontsize);
axis square; 
box on;
hold off;


% ATP 
%
subplot(m, n, find(strcmp(figure_output,'atp')))
hold on;
if exist('atp_clim','var') 
papers = fieldnames(atp_clim); 
for i = 1: length(papers)
plot(D_clim.(papers{i}),  atp_clim.(papers{i})./atp_clim.(papers{i})(1) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end  
end 

if exist('atp_nlim','var') 
papers = fieldnames(atp_nlim); 
for i = 1: length(papers)
plot(D_nlim.(papers{i}),  atp_nlim.(papers{i})./atp_nlim.(papers{i})(1) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end  
end 

%plot(D, WT_y_steady_1gL(:,num_y.ae),  '-', 'Color', color_clim.kumar_lg);
plot(D, WT_y_steady_clim(:,num_y.ae)./WT_y_steady_clim(1,num_y.ae), '-', 'Color', color_clim.default, 'LineWidth',2);
plot(D, WT_y_steady_nlim(:,num_y.ae)./WT_y_steady_nlim(1,num_y.ae), '-', 'Color', color_nlim.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim1(:,num_y.ae)./WT_y_steady_cnlim1(1,num_y.ae), '-', 'Color', color_cnlim1.default, 'LineWidth',2);
plot(D, WT_y_steady_cnlim2(:,num_y.ae)./WT_y_steady_cnlim2(1,num_y.ae), '-', 'Color', color_cnlim2.default, 'LineWidth',2);


ylabel(figcabbi.y_label.ae,'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
ylim([0 2.2])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'atp')), 'fontweight', 'bold', 'FontSize', fontsize);
axis square; 
box on;
hold off;



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





% [x_1, x_2] = find(strcmp(figure_output,'REZ'));
% x1 = x_1+0.3;
% x2 = x_1+1-0.3;
%subplot(a, b, [x1+1*b x1+2*b x2+1*b x2+2*b])
subplot(m, n, find(strcmp(figure_output,'REZ')))
hold on 
if exist('prot_sec_clim','var') 
papers = fieldnames(prot_sec_clim.r); 
for i = 1: length(papers)
plot(D_clim.(papers{i}), ...
log2(prot_sec_clim.r.(papers{i})./prot_sec_clim.r.(papers{i})(1)) ,...
'o', 'Color', plt_clrs.yellow,...
'MarkerFaceColor', plt_clrs.yellow)
end     
for i = 1: length(papers)
plot(D_clim.(papers{i}), ...
log2(prot_sec_clim.e.(papers{i})./prot_sec_clim.e.(papers{i})(1)) ,...
'o', 'Color', plt_clrs.green,...
'MarkerFaceColor', plt_clrs.green)
end  
for i = 1: length(papers)
plot(D_clim.(papers{i}), ...
log2(prot_sec_clim.z.(papers{i})./prot_sec_clim.z.(papers{i})(1)) ,...
'o', 'Color', plt_clrs.gray,...
'MarkerFaceColor', plt_clrs.gray)
end  
end
plot(D, log2(R_chemo_sim_fc), '-', 'Color', plt_clrs.yellow, 'LineWidth',2);
plot(D, log2(Z_chemo_sim_fc), '-', 'Color', plt_clrs.gray, 'LineWidth',2);
plot(D, log2(E_chemo_sim_fc), '-', 'Color', plt_clrs.green, 'LineWidth',2);
ylabel(figcabbi.y_label.rez,'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
%ylim(y_lim_rez);
ylim([-1 1.45]);
yticks([-1 0 1])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
hold off
%legend([rsim, zsim, esim],'R (sim)','Z(sim)','E(sim)')
%legend([rexp, eexp, zexp],'R ','E ','Z ','Location','southeast')
%legend boxoff 
axis square; 
box on;
%text(labelx*max(xlim)/2, ((labely-1)/2+1)*max(ylim), figure_id(strcmp(figure_output,'REZ')), 'fontweight', 'bold', 'FontSize', fontsize);

subplot(m, n, find(strcmp(figure_output,'REZ_nlim')))
hold on 
if exist('prot_sec_clim','var') 
papers = fieldnames(prot_sec_nlim.r); 
for i = 1: length(papers)
plot(D_nlim.(papers{i}), ...
log2(prot_sec_nlim.r.(papers{i})./prot_sec_nlim.r.(papers{i})(1)) ,...
'o', 'Color', plt_clrs.yellow,...
'MarkerFaceColor', plt_clrs.yellow)
end     
for i = 1: length(papers)
plot(D_nlim.(papers{i}), ...
log2(prot_sec_nlim.e.(papers{i})./prot_sec_nlim.e.(papers{i})(1)) ,...
'o', 'Color', plt_clrs.green,...
'MarkerFaceColor', plt_clrs.green)
end  
for i = 1: length(papers)
plot(D_nlim.(papers{i}), ...
log2(prot_sec_nlim.z.(papers{i})./prot_sec_nlim.z.(papers{i})(1)) ,...
'o', 'Color', plt_clrs.gray,...
'MarkerFaceColor', plt_clrs.gray)
end  
end
plot(D, log2(R_chemo_sim_fc_nlim), '-', 'Color', plt_clrs.yellow, 'LineWidth',2);
plot(D, log2(Z_chemo_sim_fc_nlim), '-', 'Color', plt_clrs.gray, 'LineWidth',2);
plot(D, log2(E_chemo_sim_fc_nlim), '-', 'Color', plt_clrs.green, 'LineWidth',2);
ylabel(figcabbi.y_label.rez,'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
%ylim(y_lim_rez);
ylim([-1 1.45]);
yticks([-1 0 1])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
hold off
%legend([rsim, zsim, esim],'R (sim)','Z(sim)','E(sim)')
%legend([rexp, eexp, zexp],'R ','E ','Z ','Location','southeast')
%legend boxoff 
axis square; 
box on;
%text(labelx*max(xlim)/2, ((labely-1)/2+1)*max(ylim), figure_id(strcmp(figure_output,'REZ')), 'fontweight', 'bold', 'FontSize', fontsize);


set(gcf,'Position',position)

