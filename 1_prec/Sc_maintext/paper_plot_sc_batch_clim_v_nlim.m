function paper_plot_sc_batch_clim_v_nlim(clim_WT_batch_min_t, nlim_WT_batch_min_t,clim_WT_batch_YPD_t,   ...
clim_WT_batch_min_y, nlim_WT_batch_min_y,clim_WT_batch_YPD_y, ...
clim_WT_batch_min_sig, nlim_WT_batch_min_sig,clim_WT_batch_YPD_sig, ...
clim_WT_batch_min_met_reac, nlim_WT_batch_min_met_reac,clim_WT_batch_YPD_met_reac, ...
par, num_y, num_prot, num_flux, figure_output_format, figure_position )

%% load data
yeast_type = 'sc';
data_unit_conv;

[glucose_clim, ethanol_clim,  cell_clim, ...
precursor_clim, atp_clim, aa_clim, prot_sec_clim,  ...
tag_clim, pl_clim, paperinfo_clim,...
color_clim, shape_clim] =  SC_C_lim_batch(); 

lipid_c = lipid_sc_c; 

%% plotting stuff

plt_clrs = plot_colors; 

paper_figure_color; 
paper_figure_label; 

x_label      = 'time (h)';
x_lim        = [0 max(clim_WT_batch_min_t)]; 
x_lim_long   = [0 35]; 

ylim_prot    = [0 60]; 
ylim_RNA     = [0 20];
ylim_lipid   = [0 70];
ylim_clim      = [0 20]; 
figure_output = figure_output_format;
%figure_id     = fig_batch.figure_id;

[b, a]     = size(figure_output);
position   = figure_position; 

%[206 315 530 749]
y_lim_prot = [-1 3.5];
% Patch area 


color_clim.default    = [102, 194, 165]./255; % color_clim.climolomeo;
color_nlim.default    = [31,  120, 180]./255; % color_clim.nlimar;
plt_clrs.darkgreen   = '#419945';
color_cnlim1.default = '#51c5b8';%'#368E6A';
color_cnlim2.default = '#0095dd';%'#2B838F';
%_% plt_clrs.green   = '#61c95e';'#1a9c53';
%_% plt_clrs.blue   = '#0e80be';'#1f8cbe';
%_% color_cnlim1.default = '#38b07b';'#4cbfb0';'#56afb2';
%_% color_cnlim2.default = '#3eb0c8';'#36a4c5';'#3aa9d1';

if  strcmp(yeast_type,'sc')
color_clim.default    = plt_clrs.green;%_% '#419945'; 
color_clim.bartolomeo = plt_clrs.green; 
color_clim.murphy     = plt_clrs.green;   
color_clim.zampar     = plt_clrs.green;  
color_clim.solis      = plt_clrs.green;
color_clim.casanovas  = plt_clrs.green;


shape_clim.casanovas  ='s';
shape_clim.zampar     ='o';
shape_clim.murphy     ='d';
shape_clim.bartolomeo ='^';
shape_clim.solis      ='v';    
x_left   = 7; 
x_middle = 14;
end


if  strcmp(yeast_type,'rt')

color_clim.default    = '#419945'; 
color_clim.jagtap     = plt_clrs.green; 
color_clim.kim        = plt_clrs.green; 
color_clim.shi        = plt_clrs.green; 

shape_clim.default    = 'o'; 
shape_clim.jagtap     = 's'; 
shape_clim.kim        = 'd'; 
shape_clim.shi        = '^'; 


color_nlim.shi        = plt_clrs.blue; 
color_nlim.evans_nh4  = plt_clrs.blue; 
color_nlim.evans_glu  = plt_clrs.blue; 
color_nlim.evansb_nh4 = plt_clrs.blue; 
color_nlim.evansb_glu = plt_clrs.blue; 
color_nlim.tiukova    = plt_clrs.blue; 

shape_nlim.shi        ='o';
shape_nlim.evansb_glu ='^';
shape_nlim.evansb_nh4 ='^';   
shape_nlim.evans_glu  ='^';
shape_nlim.evans_nh4  ='^';     
shape_nlim.tiukova    ='s';  

x_left   = 60; 
x_middle = 60;
end


%% REZ fraction calculation before normalizing
% clim_WT_batch_min_t, nlim_WT_batch_min_t, clim_WT_batch_YPD_t,
clim_YPD_total_protein_con = sum(clim_WT_batch_YPD_y(:,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);
nlim_total_protein_con   = sum(nlim_WT_batch_min_y(:,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);
clim_total_protein_con   = sum(clim_WT_batch_min_y(:,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);

% fraction
R_clim_YPD_sim = clim_WT_batch_YPD_y(:,num_y.r).*par.l(num_prot.r)'./clim_YPD_total_protein_con;
Z_clim_YPD_sim = clim_WT_batch_YPD_y(:,num_y.z).*par.l(num_prot.z)'./clim_YPD_total_protein_con;
E_clim_YPD_sim = sum(clim_WT_batch_YPD_y(:,num_y.gy:num_y.lo).*par.l(num_prot.gy:num_prot.lo)',2)./clim_YPD_total_protein_con;     

R_clim_sim = clim_WT_batch_min_y(:,num_y.r).*par.l(num_prot.r)'./clim_total_protein_con;
Z_clim_sim = clim_WT_batch_min_y(:,num_y.z).*par.l(num_prot.z)'./clim_total_protein_con;
E_clim_sim = sum(clim_WT_batch_min_y(:,num_y.gy:num_y.lo).*par.l(num_prot.gy:num_prot.lo)',2)./clim_total_protein_con;

R_nlim_sim = nlim_WT_batch_min_y(:,num_y.r).*par.l(num_prot.r)'./nlim_total_protein_con;
Z_nlim_sim = nlim_WT_batch_min_y(:,num_y.z).*par.l(num_prot.z)'./nlim_total_protein_con;
E_nlim_sim = sum(nlim_WT_batch_min_y(:,num_y.gy:num_y.lo).*par.l(num_prot.gy:num_prot.lo)',2)./nlim_total_protein_con;



%% normalizing  
clim_WT_batch_min_y(:,num_y.pc)   = clim_WT_batch_min_y(:,num_y.pc)   ./ clim_WT_batch_min_y(1,num_y.pc);
nlim_WT_batch_min_y(:,num_y.pc)   = nlim_WT_batch_min_y(:,num_y.pc)   ./ nlim_WT_batch_min_y(1,num_y.pc);
%clim_WT_batch_YPD_y(:,num_y.pc) = clim_WT_batch_YPD_y(:,num_y.pc) ./ clim_WT_batch_YPD_y(1,num_y.pc);

clim_WT_batch_min_y(:,num_y.aa_in)   = clim_WT_batch_min_y(:,num_y.aa_in)   ./ clim_WT_batch_min_y(1,num_y.aa_in);
nlim_WT_batch_min_y(:,num_y.aa_in)   = nlim_WT_batch_min_y(:,num_y.aa_in)   ./ nlim_WT_batch_min_y(1,num_y.aa_in);
%clim_WT_batch_YPD_y(:,num_y.aa_in) = clim_WT_batch_YPD_y(:,num_y.aa_in) ./ clim_WT_batch_YPD_y(1,num_y.aa_in);

clim_WT_batch_min_y(:,num_y.ae)   = clim_WT_batch_min_y(:,num_y.ae)   ./ clim_WT_batch_min_y(1,num_y.ae);
nlim_WT_batch_min_y(:,num_y.ae)   = nlim_WT_batch_min_y(:,num_y.ae)   ./ nlim_WT_batch_min_y(1,num_y.ae);
%clim_WT_batch_YPD_y(:,num_y.ae) = clim_WT_batch_YPD_y(:,num_y.ae) ./ clim_WT_batch_YPD_y(1,num_y.ae);

% proteins frac (sim)
%_% clim_WT_batch_min_y(:,[num_y.r: num_y.sd]) = (clim_WT_batch_min_y(:,[num_y.r: num_y.sd]).*par.l')./ clim_total_protein_con;
%_% nlim_WT_batch_min_y(:,[num_y.r: num_y.sd]) = (nlim_WT_batch_min_y(:,[num_y.r: num_y.sd]).*par.l')./ nlim_total_protein_con;
%_% %clim_WT_batch_YPD_y(:,[num_y.r: num_y.sd]) = (clim_WT_batch_YPD_y(:,[num_y.r: num_y.sd]).*par.l')./ clim_total_protein_con;
clim_WT_batch_min_y(:,[num_y.r: num_y.lo]) = (clim_WT_batch_min_y(:,[num_y.r: num_y.lo]).*par.l')./ clim_total_protein_con;
nlim_WT_batch_min_y(:,[num_y.r: num_y.lo]) = (nlim_WT_batch_min_y(:,[num_y.r: num_y.lo]).*par.l')./ nlim_total_protein_con;
%clim_WT_batch_YPD_y(:,[num_y.r: num_y.lo]) = (clim_WT_batch_YPD_y(:,[num_y.r: num_y.lo]).*par.l')./ clim_total_protein_con;

% protein log2fc 
%_% clim_WT_batch_min_y(:,[num_y.r: num_y.sd])   = log2(clim_WT_batch_min_y(:,[num_y.r: num_y.sd])   ./ clim_WT_batch_min_y(1,[num_y.r: num_y.sd]));
%_% nlim_WT_batch_min_y(:,[num_y.r: num_y.sd])   = log2(nlim_WT_batch_min_y(:,[num_y.r: num_y.sd])   ./ nlim_WT_batch_min_y(1,[num_y.r: num_y.sd]));
%_% %clim_WT_batch_YPD_y(:,[num_y.r: num_y.sd]) = log2(clim_WT_batch_YPD_y(:,[num_y.r: num_y.sd]) ./ clim_WT_batch_YPD_y(1,[num_y.r: num_y.sd]));
clim_WT_batch_min_y(:,[num_y.r: num_y.lo])   = log2(clim_WT_batch_min_y(:,[num_y.r: num_y.lo])   ./ clim_WT_batch_min_y(1,[num_y.r: num_y.lo]));
nlim_WT_batch_min_y(:,[num_y.r: num_y.lo])   = log2(nlim_WT_batch_min_y(:,[num_y.r: num_y.lo])   ./ nlim_WT_batch_min_y(1,[num_y.r: num_y.lo]));
%clim_WT_batch_YPD_y(:,[num_y.r: num_y.lo]) = log2(clim_WT_batch_YPD_y(:,[num_y.r: num_y.lo]) ./ clim_WT_batch_YPD_y(1,[num_y.r: num_y.lo]));






% ------------------------------- plotting --------------------------------
% ------------------------------------------------------------------------- 
figure; 
%% real metabolites
% glucose
subplot(a, b, find(strcmp(figure_output,'glucose')))
hold on;
if exist('glucose_clim','var')    
papers = fieldnames(glucose_clim);
for i = 1: length(papers)
plot(glucose_clim.(papers{i})(:,1),  glucose_clim.(papers{i})(:,2) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end   
end 
if exist('glucose_nlim','var')    
papers = fieldnames(glucose_nlim);
for i = 1: length(papers)
plot(glucose_nlim.(papers{i})(:,1),  glucose_nlim.(papers{i})(:,2) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end 
plot(clim_WT_batch_min_t,   clim_WT_batch_min_y(:,num_y.gl),   'Color', color_clim.default, 'LineWidth',2);
plot(nlim_WT_batch_min_t,   nlim_WT_batch_min_y(:,num_y.gl),   'Color', color_nlim.default, 'LineWidth',2);
%plot(clim_WT_batch_YPD_t, clim_WT_batch_YPD_y(:,num_y.gl), 'Color', color_clim.murphy);

ylabel(fig_batch.y_label.glu,'Color', xy_label_color)
xlabel(x_label,'Color', xy_label_color);
xlim(x_lim);
ylim([0 1.6*10^5]) %_%
yticks([0 0.8*10^5 1.6*10^5]) %_%
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2) 
axis square; 
box on;
hold off;
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'glucose')), 'fontweight', 'bold', 'FontSize', fontsize);

% ethanol 
%
subplot(a, b, find(strcmp(figure_output,'ethanol')))
hold on;
if exist('ethanol_clim','var')    
papers = fieldnames(ethanol_clim);
for i = 1: length(papers)
plot(ethanol_clim.(papers{i})(:,1),  ethanol_clim.(papers{i})(:,2) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end  
end 
plot(clim_WT_batch_min_t,   clim_WT_batch_min_y(:,num_y.eh),   'Color', color_clim.default, 'LineWidth',2);
plot(nlim_WT_batch_min_t,   nlim_WT_batch_min_y(:,num_y.eh),   'Color', color_nlim.default, 'LineWidth',2);
%plot(clim_WT_batch_YPD_t,   clim_WT_batch_YPD_y(:,num_y.eh), 'Color', color_clim.murphy);

ylabel(fig_batch.y_label.eth,'Color', xy_label_color)
xlabel(x_label,'Color', xy_label_color);
xlim(x_lim);
ylim([0 2.5*10^5])
% yticks([0 0.9*10^5 1.8*10^5])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
axis square; 
box on;
hold off;
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'ethanol')), 'fontweight', 'bold', 'FontSize', fontsize);
%}

% nh4 uptake rate
%{
subplot(a, b, find(strcmp(figure_output,'JNH4')))
hold on;
plot(clim_WT_batch_min_t,   par.n_nh *clim_WT_batch_min_met_reac.flux(:,num_flux.as),   'Color', color_clim.default);
plot(nlim_WT_batch_min_t,   par.n_nh *nlim_WT_batch_min_met_reac.flux(:,num_flux.as),   'Color', color_nlim.default);
%plot(clim_WT_batch_YPD_t, par.n_nh *clim_WT_batch_YPD_met_reac.flux(:,num_flux.as), 'Color', color_clim.murphy);

ylabel(fig2.y_label.jnh4,'Color',xy_label_color)
xlabel(x_label,'Color',xy_label_color);
xlim(x_lim);
ylim([0 1.1*max(ylim)])
%max_y = max([max(ylim), max(kumar_lg.Jgy), max(kumar_hg.Jgy), max(boer.Jgy), max(hackett.Jgy)]); 
patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'glucose')), 'fontweight', 'bold', 'FontSize', fontsize);
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
plot(cell_clim.(papers{i})(:,1),  cell_clim.(papers{i})(:,2) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end  
end 
if exist('cell_nlim','var')    
papers = fieldnames(cell_nlim);
for i = 1: length(papers)
plot(cell_nlim.(papers{i})(:,1),  cell_nlim.(papers{i})(:,2) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end  
end 
plot(clim_WT_batch_min_t,   clim_WT_batch_min_y(:,end),   'Color', color_clim.default, 'LineWidth',2);
plot(nlim_WT_batch_min_t,   nlim_WT_batch_min_y(:,end),   'Color', color_nlim.default, 'LineWidth',2);
%plot(clim_WT_batch_YPD_t, clim_WT_batch_YPD_y(:,end), 'Color', color_clim.murphy);

ylabel(fig_batch.y_label.cell,'Color', xy_label_color)
xlabel(x_label,'Color', xy_label_color);
xlim(x_lim);
%ylim([0 6.0*10^11])
yticks([0 5 10])%_%
ylim([0 10])%_%
%_% ylim([0 1.2*max(ylim)])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
axis square; 
box on;
hold off;
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'cells')), 'fontweight', 'bold', 'FontSize', fontsize);

%% virtual metabolites 
% 
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
plot(precursor_clim.(prec).(papers{i})(:,1),...
precursor_clim.(prec).(papers{i})(:,2)./precursor_clim.(prec).(papers{i})(1,2) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end 
end 

plot(clim_WT_batch_min_t,   clim_WT_batch_min_y(:,num_y.pc),   'Color', color_clim.default, 'LineWidth',2);
plot(nlim_WT_batch_min_t,   nlim_WT_batch_min_y(:,num_y.pc),   'Color', color_nlim.default, 'LineWidth',2);
%plot(clim_WT_batch_YPD_t, clim_WT_batch_YPD_y(:,num_y.pc), 'Color', color_clim.murphy);

ylabel(fig_batch.y_label.pc,'Color', xy_label_color)
xlabel(x_label,'Color', xy_label_color);
xlim(x_lim);
ylim([0 20])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
axis square; 
box on;
hold off;
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'precursor')), 'fontweight', 'bold', 'FontSize', fontsize);

% 
% % precursor 
% %
% subplot(a, b, find(strcmp(figure_output,'precursor')))
% hold on;
% 
% plot(clim_WT_batch_min_t,   clim_WT_batch_min_y(:,num_y.pc),   'Color', color_clim.default);
% plot(nlim_WT_batch_min_t,   nlim_WT_batch_min_y(:,num_y.pc),   'Color', color_nlim.default);
% %plot(clim_WT_batch_YPD_t, clim_WT_batch_YPD_y(:,num_y.pc), 'Color', color_clim.murphy);
% 
% ylabel(fig_batch.y_label.pc,'Color', xy_label_color)
% xlabel(x_label,'Color', xy_label_color);
% xlim(x_lim);
% ylim([0 1.05*max(ylim)])
% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
% axis square; 
% box on;
% hold off;
% %text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'precursor')), 'fontweight', 'bold', 'FontSize', fontsize);


%}

% intracellular amino acids 
%
subplot(a, b, find(strcmp(figure_output,'aa_in')))
hold on;
if exist('aa_clim','var')    
papers = fieldnames(aa_clim.glutamine);
for i = 1: length(papers)
plot(aa_clim.glutamine.(papers{i})(:,1),...
aa_clim.glutamine.(papers{i})(:,2)./ aa_clim.glutamine.(papers{i})(1,2) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end 
plot(clim_WT_batch_min_t,   clim_WT_batch_min_y(:,num_y.aa_in),   'Color', color_clim.default, 'LineWidth',2);
plot(nlim_WT_batch_min_t,   nlim_WT_batch_min_y(:,num_y.aa_in),   'Color', color_nlim.default, 'LineWidth',2);
%plot(clim_WT_batch_YPD_t, clim_WT_batch_YPD_y(:,num_y.aa_in), 'Color', color_clim.murphy);

ylabel(fig_batch.y_label.aa,'Color', xy_label_color)
xlabel(x_label,'Color', xy_label_color);
xlim(x_lim);
ylim([0 1.5])
%yticks([0 1 2 ])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
axis square; 
box on;
hold off;
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'aa_in')), 'fontweight', 'bold', 'FontSize', fontsize);


%}

% ATP 
%
subplot(a, b, find(strcmp(figure_output,'atp')))
hold on;
if exist('atp_clim','var')    
papers = fieldnames(atp_clim);
for i = 1: length(papers)
plot(atp_clim.(papers{i})(:,1),...
atp_clim.(papers{i})(:,2)./ atp_clim.(papers{i})(1,2) ,...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end
%plot(atp_clim.nlimar(:,1),         atp_clim.nlimar(:,2),    shape_clim.nlimar, 'Color', color_nlim.default, 'MarkerFaceColor', color_nlim.default)

plot(clim_WT_batch_min_t,   clim_WT_batch_min_y(:,num_y.ae),   'Color', color_clim.default, 'LineWidth',2);
plot(nlim_WT_batch_min_t,   nlim_WT_batch_min_y(:,num_y.ae),   'Color', color_nlim.default, 'LineWidth',2);
%plot(clim_WT_batch_YPD_t, clim_WT_batch_YPD_y(:,num_y.ae), 'Color', color_clim.murphy);

ylabel(fig_batch.y_label.ae,'Color', xy_label_color)
xlabel(x_label,'Color', xy_label_color);
xlim(x_lim);
ylim([0 4])
yticks([0 2 4])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
axis square; 
box on;
hold off;
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'atp')), 'fontweight', 'bold', 'FontSize', fontsize);


% NH4 
% 
subplot(a, b, find(strcmp(figure_output,'NH4')))
hold on;
if exist('nh4_nlim','var')    
papers = fieldnames(nh4_nlim);
for i = 1: length(papers)
plot(nh4_nlim.(papers{i})(:,1),  nh4_nlim.(papers{i})(:,2) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end     

plot(clim_WT_batch_min_t,   clim_WT_batch_min_y(:,num_y.nh4),   'Color', color_clim.default, 'LineWidth',2);
plot(nlim_WT_batch_min_t,   nlim_WT_batch_min_y(:,num_y.nh4),   'Color', color_nlim.default, 'LineWidth',2);
%plot(clim_WT_batch_YPD_t, clim_WT_batch_YPD_y(:,num_y.nh4), 'Color', color_clim.murphy);

ylabel(fig_batch.y_label.nh4,'Color', xy_label_color)
xlabel(x_label,'Color', xy_label_color);
xlim(x_lim);
ylim([0, 1.1*max(ylim)])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
axis square; 
box on;
hold off;



%}

%% signaling 

% SNF1 
subplot(a, b, find(strcmp(figure_output,'snf1')))
hold on;

plot(clim_WT_batch_min_t,   clim_WT_batch_min_sig.snf1(:,1)   ./par.s_tot, 'Color', color_clim.default, 'LineWidth',2);
plot(nlim_WT_batch_min_t,   nlim_WT_batch_min_sig.snf1(:,1)   ./par.s_tot, 'Color', color_nlim.default, 'LineWidth',2);
%plot(clim_WT_batch_YPD_t, clim_WT_batch_YPD_sig.snf1(:,1) ./par.s_tot, 'Color', color_clim.murphy);
ylabel(fig_batch.y_label.snf1,'Color', xy_label_color)
xlabel(x_label,'Color', xy_label_color);
xlim(x_lim);
ylim([0 1.1])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%ylim([0 5])
axis square; 
box on;
hold off; 
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'snf1')), 'fontweight', 'bold', 'FontSize', fontsize);


% TORC1
subplot(a, b, find(strcmp(figure_output,'tor')))
hold on;

plot(clim_WT_batch_min_t,   clim_WT_batch_min_sig.tor(:,1)  ./ par.tau_tot, 'Color', color_clim.default, 'LineWidth',2);
plot(nlim_WT_batch_min_t,   nlim_WT_batch_min_sig.tor(:,1)  ./ par.tau_tot, 'Color', color_nlim.default, 'LineWidth',2);
%plot(clim_WT_batch_YPD_t, clim_WT_batch_YPD_sig.tor(:,1)./ par.tau_tot, 'Color', color_clim.murphy);
ylabel(fig_batch.y_label.tor,'Color', xy_label_color)
xlabel(x_label,'Color', xy_label_color);
xlim(x_lim);
ylim([0 1.1*max(ylim)])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim),  color1, color2, alpha1, alpha2)
axis square; 
box on;
hold off; 
%text(labelx*max(xlim), labely*max(ylim), figure_id(strcmp(figure_output,'tor')), 'fontweight', 'bold', 'FontSize', fontsize);

%% proteins 

%subplot(a, b, [x1+b x1+2*b x2+b x2+2*b])
% REZ
subplot(a, b, find(strcmp(figure_output,'REZ')))
hold on 
if exist('prot_sec_clim','var')    
papers = fieldnames(prot_sec_clim.r); 
for i = 1: length(papers)
plot(prot_sec_clim.time_prot.(papers{i}), ...
prot_sec_clim.r.(papers{i}),...
'o', 'Color', plt_clrs.yellow,...
'MarkerFaceColor', plt_clrs.yellow)
end     
for i = 1: length(papers)
plot(prot_sec_clim.time_prot.(papers{i}), ...
prot_sec_clim.e.(papers{i}),...
'o', 'Color', plt_clrs.green,...
'MarkerFaceColor', plt_clrs.green)
end  
for i = 1: length(papers)
plot(prot_sec_clim.time_prot.(papers{i}), ...
prot_sec_clim.z.(papers{i}),...
'o', 'Color', plt_clrs.gray,...
'MarkerFaceColor', plt_clrs.gray)
end         
end 
rsim = plot(clim_WT_batch_YPD_t, R_clim_YPD_sim, '-', 'Color', plt_clrs.yellow, 'LineWidth',2);
esim = plot(clim_WT_batch_YPD_t, E_clim_YPD_sim, '-', 'Color', plt_clrs.green, 'LineWidth',2);
zsim = plot(clim_WT_batch_YPD_t, Z_clim_YPD_sim, '-', 'Color', plt_clrs.gray, 'LineWidth',2);
ylabel(fig_batch.y_label.rez,'Color', xy_label_color)
xlabel(x_label,'Color', xy_label_color);
xlim(x_lim);
ylim([0 0.6]);
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
hold off
axis square; 
box on;
hold off;
%text(labelx*max(xlim)/2, ((labely-1)/2+1)*max(ylim), figure_id(strcmp(figure_output,'REZ')), 'fontweight', 'bold', 'FontSize', fontsize);


% gy
subplot(a, b, find(strcmp(figure_output,'Egy')))
hold on;
if exist('prot_sec_clim','var')    
papers = fieldnames(prot_sec_clim.gy); 
for i = 1: length(papers)
plot(prot_sec_clim.time_prot.(papers{i}), ...
log2(prot_sec_clim.gy.(papers{i})./prot_sec_clim.gy.(papers{i})(1)), ...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end
if exist('prot_sec_nlim','var')    
papers = fieldnames(prot_sec_nlim.gy); 
for i = 1: length(papers)
plot(prot_sec_nlim.time_prot.(papers{i}), ...
log2(prot_sec_nlim.gy.(papers{i})./prot_sec_nlim.gy.(papers{i})(1)), ...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end 
plot(clim_WT_batch_min_t,   clim_WT_batch_min_y(:,num_y.gy),   'Color', color_clim.default, 'LineWidth',2);
plot(nlim_WT_batch_min_t,   nlim_WT_batch_min_y(:,num_y.gy),   'Color', color_nlim.default, 'LineWidth',2);
%plot(clim_WT_batch_YPD_t, clim_WT_batch_YPD_y(:,num_y.gy), 'Color', color_clim.murphy);
ylabel(fig_batch.y_label.gy,'Color', xy_label_color)
xlabel(x_label,'Color', xy_label_color);
xlim(x_lim);
ylim(y_lim_prot);
yticks([-1 1 3])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
axis square; 
box on;
hold off;
%text(labelx*max(xlim), min(ylim)+labely*(max(ylim)-min(ylim)), figure_id(strcmp(figure_output,'Egy')), 'fontweight', 'bold', 'FontSize', fontsize);

% fe 
subplot(a, b,  find(strcmp(figure_output,'Efe')))
hold on;
if exist('prot_sec_clim','var')    
papers = fieldnames(prot_sec_clim.fe); 
for i = 1: length(papers)
plot(prot_sec_clim.time_prot.(papers{i}), ...
log2(prot_sec_clim.fe.(papers{i})./prot_sec_clim.fe.(papers{i})(1)), ...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end 
if exist('prot_sec_nlim','var')    
papers = fieldnames(prot_sec_nlim.fe); 
for i = 1: length(papers)
plot(prot_sec_nlim.time_prot.(papers{i}), ...
log2(prot_sec_nlim.fe.(papers{i})./prot_sec_nlim.fe.(papers{i})(1)), ...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end 
plot(clim_WT_batch_min_t,   clim_WT_batch_min_y(:,num_y.fe),   'Color', color_clim.default, 'LineWidth',2);
plot(nlim_WT_batch_min_t,   nlim_WT_batch_min_y(:,num_y.fe),   'Color', color_nlim.default, 'LineWidth',2);
%plot(clim_WT_batch_YPD_t, clim_WT_batch_YPD_y(:,num_y.fe), 'Color', color_clim.murphy);
ylabel(fig_batch.y_label.fe,'Color', xy_label_color)
xlabel(x_label,'Color', xy_label_color);
xlim(x_lim);
%ylim([-1 1])
ylim(y_lim_prot);
yticks([-1 1 3])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
axis square; 
box on;
hold off;
%text(labelx*max(xlim),  min(ylim)+labely*(max(ylim)-min(ylim)), figure_id(strcmp(figure_output,'Efe')), 'fontweight', 'bold', 'FontSize', fontsize);

% gn 
subplot(a, b, find(strcmp(figure_output,'Egn')))
hold on;
if exist('prot_sec_clim','var')    
papers = fieldnames(prot_sec_clim.gn); 
for i = 1: length(papers)
plot(prot_sec_clim.time_prot.(papers{i}), ...
log2(prot_sec_clim.gn.(papers{i})./prot_sec_clim.gn.(papers{i})(1)), ...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end 
if exist('prot_sec_nlim','var')    
papers = fieldnames(prot_sec_nlim.gn); 
for i = 1: length(papers)
plot(prot_sec_nlim.time_prot.(papers{i}), ...
log2(prot_sec_nlim.gn.(papers{i})./prot_sec_nlim.gn.(papers{i})(1)), ...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end

plot(clim_WT_batch_min_t,   clim_WT_batch_min_y(:,num_y.gn),   'Color', color_clim.default, 'LineWidth',2);
plot(nlim_WT_batch_min_t,   nlim_WT_batch_min_y(:,num_y.gn),   'Color', color_nlim.default, 'LineWidth',2);
%plot(clim_WT_batch_YPD_t, clim_WT_batch_YPD_y(:,num_y.gn), 'Color', color_clim.murphy);
ylabel(fig_batch.y_label.gn,'Color', xy_label_color)
xlabel(x_label,'Color', xy_label_color);
xlim(x_lim);
%ylim([-1 3])
ylim(y_lim_prot);
yticks([-1 1 3])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
axis square; 
box on;
hold off;
%text(labelx*max(xlim),  min(ylim)+labely*(max(ylim)-min(ylim)), figure_id(strcmp(figure_output,'Egn')), 'fontweight', 'bold', 'FontSize', fontsize);

% mt 
subplot(a, b, find(strcmp(figure_output,'Emt')))
hold on;
if exist('prot_sec_clim','var')    
papers = fieldnames(prot_sec_clim.mt); 
for i = 1: length(papers)
plot(prot_sec_clim.time_prot.(papers{i}), ...
log2(prot_sec_clim.mt.(papers{i})./prot_sec_clim.mt.(papers{i})(1)), ...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end
if exist('prot_sec_nlim','var')    
papers = fieldnames(prot_sec_nlim.mt); 
for i = 1: length(papers)
plot(prot_sec_nlim.time_prot.(papers{i}), ...
log2(prot_sec_nlim.mt.(papers{i})./prot_sec_nlim.mt.(papers{i})(1)), ...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end
end

plot(clim_WT_batch_min_t,   clim_WT_batch_min_y(:,num_y.mt),   'Color', color_clim.default, 'LineWidth',2);
plot(nlim_WT_batch_min_t,   nlim_WT_batch_min_y(:,num_y.mt),   'Color', color_nlim.default, 'LineWidth',2);
%plot(clim_WT_batch_YPD_t, clim_WT_batch_YPD_y(:,num_y.mt), 'Color', color_clim.murphy);
ylabel(fig_batch.y_label.mt,'Color', xy_label_color)
xlabel(x_label,'Color', xy_label_color);
xlim(x_lim);

ylim(y_lim_prot);
yticks([-1 1 3])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
axis square; 
box on;
hold off;
%text(labelx*max(xlim),  min(ylim)+labely*(max(ylim)-min(ylim)), figure_id(strcmp(figure_output,'Emt')), 'fontweight', 'bold', 'FontSize', fontsize);

% as 
subplot(a, b, find(strcmp(figure_output,'Eas')))
hold on;
if exist('prot_sec_clim','var')    
papers = fieldnames(prot_sec_clim.as); 
for i = 1: length(papers)
plot(prot_sec_clim.time_prot.(papers{i}), ...
log2(prot_sec_clim.as.(papers{i})./prot_sec_clim.as.(papers{i})(1)), ...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end 
if exist('prot_sec_nlim','var')    
papers = fieldnames(prot_sec_nlim.as); 
for i = 1: length(papers)
plot(prot_sec_nlim.time_prot.(papers{i}), ...
log2(prot_sec_nlim.as.(papers{i})./prot_sec_nlim.as.(papers{i})(1)), ...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end
end 

plot(clim_WT_batch_min_t,   clim_WT_batch_min_y(:,num_y.as),   'Color', color_clim.default, 'LineWidth',2);
plot(nlim_WT_batch_min_t,   nlim_WT_batch_min_y(:,num_y.as),   'Color', color_nlim.default, 'LineWidth',2);
%plot(clim_WT_batch_YPD_t, clim_WT_batch_YPD_y(:,num_y.as), 'Color', color_clim.murphy);
ylabel(fig_batch.y_label.as,'Color', xy_label_color)
xlabel(x_label,'Color', xy_label_color);
xlim(x_lim);
%ylim([-1 1])
ylim(y_lim_prot)
yticks([-1 1 3])

%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)

axis square; 
box on;
hold off;
%text(labelx*max(xlim),  min(ylim)+labely*(max(ylim)-min(ylim)), figure_id(strcmp(figure_output,'Eas')), 'fontweight', 'bold', 'FontSize', fontsize);

% at 
%{

subplot(a, b, find(strcmp(figure_output,'Eat')))

hold on;
if exist('prot_sec_clim','var')    
papers = fieldnames(prot_sec_clim.at); 
for i = 1: length(papers)
plot(prot_sec_clim.time_prot.(papers{i}), ...
log2(prot_sec_clim.at.(papers{i})./prot_sec_clim.at.(papers{i})(1)), ...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end 
if exist('prot_sec_nlim','var')    
papers = fieldnames(prot_sec_nlim.at); 
for i = 1: length(papers)
plot(prot_sec_nlim.time_prot.(papers{i}), ...
log2(prot_sec_nlim.at.(papers{i})./prot_sec_nlim.at.(papers{i})(1)), ...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end
end
plot(clim_WT_batch_min_t,   clim_WT_batch_min_y(:,num_y.at),   'Color', color_clim.default);
plot(nlim_WT_batch_min_t,   nlim_WT_batch_min_y(:,num_y.at),   'Color', color_nlim.default);
%plot(clim_WT_batch_YPD_t, clim_WT_batch_YPD_y(:,num_y.at), 'Color', color_clim.murphy);

ylabel(fig_batch.y_label.at,'Color', xy_label_color)
xlabel(x_label,'Color', xy_label_color);
xlim(x_lim);
ylim([-2,2])
patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)

axis square; 
box on;
hold off;
%text(labelx*max(xlim),  min(ylim)+labely*(max(ylim)-min(ylim)), figure_id(strcmp(figure_output,'Eat')), 'fontweight', 'bold', 'FontSize', fontsize);
%}


% lp
%
subplot(a, b, find(strcmp(figure_output,'Elp')))
hold on;
if exist('prot_sec_clim','var')    
papers = fieldnames(prot_sec_clim.lp); 
for i = 1: length(papers)
plot(prot_sec_clim.time_prot.(papers{i}), ...
log2(prot_sec_clim.lp.(papers{i})./prot_sec_clim.lp.(papers{i})(1)), ...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end
if exist('prot_sec_nlim','var')    
papers = fieldnames(prot_sec_nlim.lp); 
for i = 1: length(papers)
plot(prot_sec_nlim.time_prot.(papers{i}), ...
log2(prot_sec_nlim.lp.(papers{i})./prot_sec_nlim.lp.(papers{i})(1)), ...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end
end 
plot(clim_WT_batch_min_t,   clim_WT_batch_min_y(:,num_y.lp),   'Color', color_clim.default, 'LineWidth',2);
plot(nlim_WT_batch_min_t,   nlim_WT_batch_min_y(:,num_y.lp),   'Color', color_nlim.default, 'LineWidth',2);
%plot(clim_WT_batch_YPD_t, clim_WT_batch_YPD_y(:,num_y.lp), 'Color', color_clim.murphy);

ylabel(fig_batch.y_label.lp,'Color', xy_label_color)
xlabel(x_label,'Color', xy_label_color);
xlim(x_lim);
ylim(y_lim_prot)
yticks([-1 1 3])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)

axis square; 
box on;
hold off;

% lo
%
subplot(a, b, find(strcmp(figure_output,'Elo')))
hold on;
if exist('prot_sec_clim','var')    
papers = fieldnames(prot_sec_clim.lo); 
for i = 1: length(papers)
plot(prot_sec_clim.time_prot.(papers{i}), ...
log2(prot_sec_clim.lo.(papers{i})./prot_sec_clim.lo.(papers{i})(1)), ...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end
if exist('prot_sec_nlim','var')    
papers = fieldnames(prot_sec_nlim.lo); 
for i = 1: length(papers)
plot(prot_sec_nlim.time_prot.(papers{i}), ...
log2(prot_sec_nlim.lo.(papers{i})./prot_sec_nlim.lo.(papers{i})(1)), ...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end
end 
plot(clim_WT_batch_min_t,   clim_WT_batch_min_y(:,num_y.lo),   'Color', color_clim.default, 'LineWidth',2);
plot(nlim_WT_batch_min_t,   nlim_WT_batch_min_y(:,num_y.lo),   'Color', color_nlim.default, 'LineWidth',2);
%plot(clim_WT_batch_YPD_t, clim_WT_batch_YPD_y(:,num_y.lp), 'Color', color_clim.murphy);

ylabel(fig_batch.y_label.lo,'Color', xy_label_color)
xlabel(x_label,'Color', xy_label_color);
xlim(x_lim);
ylim(y_lim_prot)
yticks([-1 1 3])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)

axis square; 
box on;
hold off;


%{
% sp
%
subplot(a, b, find(strcmp(figure_output,'Elp')))
hold on;
if exist('prot_sec_clim','var')    
papers = fieldnames(prot_sec_clim.sp); 
for i = 1: length(papers)
plot(prot_sec_clim.time_prot.(papers{i}), ...
log2(prot_sec_clim.sp.(papers{i})./prot_sec_clim.sp.(papers{i})(1)), ...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end
if exist('prot_sec_nlim','var')    
papers = fieldnames(prot_sec_nlim.sp); 
for i = 1: length(papers)
plot(prot_sec_nlim.time_prot.(papers{i}), ...
log2(prot_sec_nlim.sp.(papers{i})./prot_sec_nlim.sp.(papers{i})(1)), ...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end
end 
plot(clim_WT_batch_min_t,   clim_WT_batch_min_y(:,num_y.sp),   'Color', color_clim.default);
plot(nlim_WT_batch_min_t,   nlim_WT_batch_min_y(:,num_y.sp),   'Color', color_nlim.default);
%plot(clim_WT_batch_YPD_t, clim_WT_batch_YPD_y(:,num_y.sp), 'Color', color_clim.murphy);

ylabel(fig_batch.y_label.sp,'Color', xy_label_color)
xlabel(x_label,'Color', xy_label_color);
xlim(x_lim);
ylim(y_lim_prot)
yticks([-1 1 3])
patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)

axis square; 
box on;
hold off;



% sd
%
subplot(a, b, find(strcmp(figure_output,'Esd')))
hold on;
if exist('prot_sec_clim','var')    
papers = fieldnames(prot_sec_clim.sd); 
for i = 1: length(papers)
plot(prot_sec_clim.time_prot.(papers{i}), ...
log2(prot_sec_clim.sd.(papers{i})./prot_sec_clim.sd.(papers{i})(1)), ...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end
if exist('prot_sec_nlim','var')    
papers = fieldnames(prot_sec_nlim.sd); 
for i = 1: length(papers)
plot(prot_sec_nlim.time_prot.(papers{i}), ...
log2(prot_sec_nlim.sd.(papers{i})./prot_sec_nlim.sd.(papers{i})(1)), ...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end
end 
plot(clim_WT_batch_min_t,   clim_WT_batch_min_y(:,num_y.sd),   'Color', color_clim.default);
plot(nlim_WT_batch_min_t,   nlim_WT_batch_min_y(:,num_y.sd),   'Color', color_nlim.default);
%plot(clim_WT_batch_YPD_t, clim_WT_batch_YPD_y(:,num_y.sd), 'Color', color_clim.murphy);

ylabel(fig_batch.y_label.sd,'Color', xy_label_color)
xlabel(x_label,'Color', xy_label_color);
xlim(x_lim);
ylim(y_lim_prot)
yticks([-1 1 3])
patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)

axis square; 
box on;
hold off;
%}

% r
%
subplot(a, b, find(strcmp(figure_output,'R')))
hold on;
if exist('prot_sec_clim','var')    
papers = fieldnames(prot_sec_clim.r); 
for i = 1: length(papers)
plot(prot_sec_clim.time_prot.(papers{i}), ...
log2(prot_sec_clim.r.(papers{i})./prot_sec_clim.r.(papers{i})(1)), ...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end
if exist('prot_sec_nlim','var')    
papers = fieldnames(prot_sec_nlim.r); 
for i = 1: length(papers)
plot(prot_sec_nlim.time_prot.(papers{i}), ...
log2(prot_sec_nlim.r.(papers{i})./prot_sec_nlim.r.(papers{i})(1)), ...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end
end 
plot(clim_WT_batch_min_t,   clim_WT_batch_min_y(:,num_y.r),   'Color', color_clim.default, 'LineWidth',2);
plot(nlim_WT_batch_min_t,   nlim_WT_batch_min_y(:,num_y.r),   'Color', color_nlim.default, 'LineWidth',2);
%plot(clim_WT_batch_YPD_t, clim_WT_batch_YPD_y(:,num_y.r), 'Color', color_clim.murphy);

ylabel(fig_batch.y_label.r,'Color', xy_label_color)
xlabel(x_label,'Color', xy_label_color);
xlim(x_lim);
ylim(y_lim_prot)
yticks([-1 1 3])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)

axis square; 
box on;
hold off;

% z
%
subplot(a, b, find(strcmp(figure_output,'Z')))
hold on;
if exist('prot_sec_clim','var')    
papers = fieldnames(prot_sec_clim.z); 
for i = 1: length(papers)
plot(prot_sec_clim.time_prot.(papers{i}), ...
log2(prot_sec_clim.z.(papers{i})./prot_sec_clim.z.(papers{i})(1)), ...
shape_clim.(papers{i}),...
'Color', color_clim.(papers{i}), ...
'MarkerFaceColor', color_clim.(papers{i}))
end 
end
if exist('prot_sec_nlim','var')    
papers = fieldnames(prot_sec_nlim.z); 
for i = 1: length(papers)
plot(prot_sec_nlim.time_prot.(papers{i}), ...
log2(prot_sec_nlim.z.(papers{i})./prot_sec_nlim.z.(papers{i})(1)), ...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end
end 
plot(clim_WT_batch_min_t,   clim_WT_batch_min_y(:,num_y.z),   'Color', color_clim.default, 'LineWidth',2);
plot(nlim_WT_batch_min_t,   nlim_WT_batch_min_y(:,num_y.z),   'Color', color_nlim.default, 'LineWidth',2);
%plot(clim_WT_batch_YPD_t, clim_WT_batch_YPD_y(:,num_y.r), 'Color', color_clim.murphy);

ylabel(fig_batch.y_label.z,'Color', xy_label_color)
xlabel(x_label,'Color', xy_label_color);
xlim(x_lim);
ylim(y_lim_prot)
yticks([-1 1 3])
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)

axis square; 
box on;
hold off;


%% Biomass composition



% clim_total_protein_con = sum(clim_WT_batch_YPD_y(1,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);
% nlim_total_protein_con   = sum(nlim_WT_batch_min_y(1,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);
% clim_total_protein_con   = sum(clim_WT_batch_min_y(1,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);
% protein
subplot(a, b, find(strcmp(figure_output,'protein')))
hold on;
if exist('protein_nlim','var')    
papers = fieldnames(protein_nlim);
for i = 1: length(papers)
plot(protein_nlim.(papers{i})(:,1),  protein_nlim.(papers{i})(:,2) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end   
if exist('protein_clim','var')    
papers = fieldnames(protein_clim);
for i = 1: length(papers)
plot(protein_clim.(papers{i})(:,1),  protein_clim.(papers{i})(:,2) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end 

plot(clim_WT_batch_min_t,   100*(clim_total_protein_con.*pt_uM_to_gPerL)./par.rho_cell,   'Color', color_clim.default, 'LineWidth',2);
plot(nlim_WT_batch_min_t,   100*(nlim_total_protein_con.*pt_uM_to_gPerL)./par.rho_cell,   'Color', color_nlim.default, 'LineWidth',2);
%plot(clim_WT_batch_YPD_t, 100*(clim_total_protein_con.*pt_uM_to_gPerL)./par.rho_cell, 'Color', color_clim.murphy);

hold off 
ylabel('protein (%)')
xlabel(x_label);
xlim(x_lim);
ylim(ylim_prot)
yticks([0 25 50])%_%
%_% patch_background(x_left, x_middle, max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), min(ylim)+labely*(max(ylim)-min(ylim)), figure_id(strcmp(figure_output,'Eat')), 'fontweight', 'bold', 'FontSize', fontsize);
axis square; 
box on;
hold off;



% lipid 
%
subplot(a, b, find(strcmp(figure_output,'lipid')))
hold on;
if exist('lipid_nlim','var')    
papers = fieldnames(lipid_nlim);
for i = 1: length(papers)
plot(lipid_nlim.(papers{i})(:,1),  lipid_nlim.(papers{i})(:,2) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end      
if exist('lipid_clim','var')    
papers = fieldnames(lipid_clim);
for i = 1: length(papers)
plot(lipid_clim.(papers{i})(:,1),  lipid_clim.(papers{i})(:,2) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end  
%_% plot(clim_WT_batch_min_t,   100.*clim_WT_batch_min_y(:,num_y.lp_e).*lp_uM_to_gPerL./(clim_WT_batch_min_y(:,num_y.cells)+10^-10)+lipid_c,   'Color', color_clim.default); %_% not yet
%_% plot(nlim_WT_batch_min_t,   100.*nlim_WT_batch_min_y(:,num_y.lp_e).*lp_uM_to_gPerL./(nlim_WT_batch_min_y(:,num_y.cells)+10^-10)+lipid_c,   'Color', color_nlim.default); %_% not yet
plot(clim_WT_batch_min_t,   100.*clim_WT_batch_min_y(:,num_y.lp_e).*lp_uM_to_gPerL./par.rho_cell +lipid_c , 'Color', color_clim.default, 'LineWidth',2);
plot(nlim_WT_batch_min_t,   100.*nlim_WT_batch_min_y(:,num_y.lp_e).*lp_uM_to_gPerL./par.rho_cell +lipid_c , 'Color', color_nlim.default, 'LineWidth',2); 

%plot(clim_WT_batch_YPD_t, 100.*clim_WT_batch_YPD_y(:,num_y.lp_e).*lp_uM_to_gPerL./(clim_WT_batch_YPD_y(:,num_y.cells)+10^-10)+lipid_clim_c, 'Color', color_clim.murphy);

ylabel('lipid (%)','Color', xy_label_color)
xlabel(x_label,'Color', xy_label_color);
xlim(x_lim);
ylim(ylim_lipid)
yticks([0 25 50])%_%
%_% patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
axis square; 
box on;
hold off;



% carbo 
%{
subplot(a, b, find(strcmp(figure_output,'carbo')))
hold on;
if exist('carbo_nlim','var')    
papers = fieldnames(carbo_nlim);
for i = 1: length(papers)
plot(carbo_nlim.(papers{i})(:,1),  carbo_nlim.(papers{i})(:,2) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end      
if exist('carbo_clim','var')    
papers = fieldnames(carbo_clim);
for i = 1: length(papers)
plot(carbo_clim.(papers{i})(:,1),  carbo_clim.(papers{i})(:,2) ,...
shape_nlim.(papers{i}),...
'Color', color_nlim.(papers{i}), ...
'MarkerFaceColor', color_nlim.(papers{i}))
end 
end  
plot(clim_WT_batch_min_t,   100.*clim_WT_batch_min_y(:,num_y.sc).*gl_uM_to_gPerL./par.rho_cell,   'Color', color_clim.default);
plot(nlim_WT_batch_min_t,   100.*nlim_WT_batch_min_y(:,num_y.sc).*gl_uM_to_gPerL./par.rho_cell,   'Color', color_nlim.default);
%plot(clim_WT_batch_YPD_t, 100.*clim_WT_batch_YPD_y(:,num_y.sc).*gl_uM_to_gPerL./par.rho_cell, 'Color', color_clim.murphy);

ylabel('storga carb (%)','Color', xy_label_color)
xlabel(x_label,'Color', xy_label_color);
xlim(x_lim);
ylim(ylim_clim)
patch_background(x_left, x_middle,  max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
axis square; 
box on;
hold off;


% RNA 
%RNA_clim_sim = 100*R_clim_sim./2; % assume to be propotional to RNA
RNA_clim_sim   = 100*R_clim_sim./2;   % assume to be propotional to RNA
RNA_nlim_sim   = 100*R_nlim_sim./2;   % assume to be propotional to RNA

subplot(a, b, find(strcmp(figure_output,'RNA')))
hold on;
plot(clim_WT_batch_min_t,   RNA_clim_sim,   'Color', color_clim.default);
plot(nlim_WT_batch_min_t,   RNA_nlim_sim,   'Color', color_nlim.default);
%plot(clim_WT_batch_YPD_t, RNA_clim_sim, 'Color', color_clim.murphy);

hold off 
ylabel('RNA (%)')
xlabel(x_label);
xlim(x_lim);
ylim(ylim_RNA)
patch_background(x_left, x_middle, max(ylim), max(x_lim), min(ylim), color1, color2, alpha1, alpha2)
%text(labelx*max(xlim), min(ylim)+labely*(max(ylim)-min(ylim)), figure_id(strcmp(figure_output,'Eat')), 'fontweight', 'bold', 'FontSize', fontsize);
axis square; 
box on;
hold off;
%}
set(gcf, 'Position', position)





%% plot legends 
% legends

figure
subplot(2,1,1)
if exist('paperinfo_clim','var')
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
num2str(paperinfo_clim.gl.(papers{i}),'%.1f'),...
'g/L', ...
', (NH4)2SO4=',...
num2str(paperinfo_clim.n.(papers{i}),'%.1f'),...
'g/L');
end

legend(legs, Interpreter="none")
legend box off
end

if exist('paperinfo_nlim','var')
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
num2str(paperinfo_nlim.gl.(papers{i}),'%.1f'),...
'g/L', ...
', (NH4)2SO4=',...
num2str(paperinfo_nlim.n.(papers{i}),'%.1f'),...
'g/L');
end

legend(legs, Interpreter="none")
legend box off
end 




end
