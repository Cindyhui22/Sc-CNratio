function plot_all_batch_w_exp(bart_WT_batch_min_t, zamp_WT_batch_min_t, murphy_WT_batch_YPD_t, ...
bart_WT_batch_min_y, zamp_WT_batch_min_y, murphy_WT_batch_YPD_y, ...
bart_WT_batch_min_met_reac, zamp_WT_batch_min_met_reac, murphy_WT_batch_YPD_met_reac, ...
bart_WT_batch_min_sig, zamp_WT_batch_min_sig, murphy_WT_batch_YPD_sig, ...
bart_WT_batch_min_prot_syn, zamp_WT_batch_min_prot_syn, murphy_WT_batch_YPD_prot_syn, ...
bart_WT_batch_min_rib, zamp_WT_batch_min_rib, murphy_WT_batch_YPD_rib, ...
bart_WT_batch_min_g_rate, zamp_WT_batch_min_g_rate, murphy_WT_batch_YPD_g_rate, ...
label, plt_clrs, par, num_y, num_prot, num_flux, x_label)

plot_num;                           
data_unit_conv;
paper_figure_label; 

% data_bartolomeo; 
% data_murphy; 
% data_zampar; 
% data_batch;
% data_4_prec_4_prec;

% [glucose_clim, ethanol_clim,  cell_clim, ...
%     precursor_clim, atp_clim, aa_clim, prot_sec_clim,  ...
%     tag_clim, pl_clim, paper_info_clim, color, shape] =  SC_C_lim_batch(); 

garcia2014.glu_snf1_batch   = [184280 0.129824561
153308 0.138166846;
136723 0.145864662;
82471  0.188327963;
53649  0.202721088;
35932  0.280236305;
11550  0.491084855;
0      0.997457931;];

color.garcia2014 = [141, 160, 203]./255; 
shape.garcia2014 = 'o'; 


% overwrite color and shape
% color.bartolomeo   = plt_clrs.green; % purple
% color.murphy       = plt_clrs.green;    
% color.zampar       = plt_clrs.green; % green
% color.solis        = plt_clrs.green; 
% color.garcia2014   = plt_clrs.green; 
% 
% shape.bartolomeo    ='o';
% shape.murphy        ='f';
% shape.zampar        ='^';
% shape.solis         ='v';
% shape.garcia2014    ='o'; 


total_protein_con_bart = sum(bart_WT_batch_min_y(:,num_y.r:num_y.sd).* par.l(num_prot.r:num_prot.sd)',2);
total_protein_con_zamp  = sum(zamp_WT_batch_min_y(:,num_y.r:num_y.sd).* par.l(num_prot.r:num_prot.sd)',2);
total_protein_con_murphy = sum(murphy_WT_batch_YPD_y(:,num_y.r:num_y.sd).* par.l(num_prot.r:num_prot.sd)',2);
%% plot 
a = 5;
b = 6;
figure;

%% real metabolites

% extracellular amino acids 
subplot(a, b, 1)
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_y(:,num_y.aaex),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_y(:,num_y.aaex),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y(:,num_y.aaex), 'Color', plt_clrs.red);
ylabel(label.met{num_met.aaex})
xlabel(x_label);
ymax = max([max(bart_WT_batch_min_y(:,num_y.aaex)), max(zamp_WT_batch_min_y(:,num_y.aaex)), max(murphy_WT_batch_YPD_y(:,num_y.aaex))]); 
if ymax > 0 
ylim([0 1.05.*ymax])
end 
xlim([0 bart_WT_batch_min_t(end)])
%xline(t_gl_depl_bart, 'Color', plt_clrs.green)
%xline(t_eh_depl_bart, 'Color', plt_clrs.green)
%xline(t_gl_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_eh_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_gl_depl_murphy, 'Color', plt_clrs.red)
%xline(t_eh_depl_murphy, 'Color', plt_clrs.red)
axis square; 
box on;
hold off;

% glucose
subplot(a, b, 1 + b)
hold on;
if exist('glucose_clim','var')  
bart_glu    = plot(glucose_clim.bartolomeo(:,1), glucose_clim.bartolomeo(:,2), shape.bartolomeo, 'Color', color.bartolomeo,   'MarkerFaceColor',color.bartolomeo);
zampar_glu  = plot(glucose_clim.zampar(:,1),     glucose_clim.zampar(:,2),     shape.zampar, 'Color', color.zampar,   'MarkerFaceColor',color.zampar);
murphy_glu  = plot(glucose_clim.murphy(:,1),     glucose_clim.murphy(:,2),     shape.murphy, 'Color', color.murphy, 'MarkerFaceColor',color.murphy);
end 
%plot(bisschops.glucose(:,1) + 10,      bisschops.glucose(:,2),      'o', 'Color', plt_clrs.orange)
%plot(solisescalante.glucose(:,1) + 10, solisescalante.glucose(:,2), 'o', 'Color', plt_clrs.pink)
plot(bart_WT_batch_min_t,   bart_WT_batch_min_y(:,num_y.gl),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_y(:,num_y.gl),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y(:,num_y.gl), 'Color', plt_clrs.red);
ylabel(label.met{num_met.gl})
xlabel(x_label);
ymax = max([max(bart_WT_batch_min_y(:,num_y.gl)), max(zamp_WT_batch_min_y(:,num_y.gl)), max(murphy_WT_batch_YPD_y(:,num_y.gl))]); 
if ymax > 0 
ylim([0 1.05.*ymax])
end 
xlim([0 bart_WT_batch_min_t(end)])
xlim([0 bart_WT_batch_min_t(end)])
%xline(t_gl_depl_bart, 'Color', plt_clrs.green)
%xline(t_eh_depl_bart, 'Color', plt_clrs.green)
%xline(t_gl_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_eh_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_gl_depl_murphy, 'Color', plt_clrs.red)
%xline(t_eh_depl_murphy, 'Color', plt_clrs.red)
axis square; 
box on;
hold off;

% ethanol 
%
subplot(a, b, 1 + 2*b)
hold on;
if exist('glucose_clim','var')  
plot(ethanol_clim.bartolomeo(:,1), ethanol_clim.bartolomeo(:,2), shape.bartolomeo, 'Color', color.bartolomeo, 'MarkerFaceColor', color.bartolomeo)
plot(ethanol_clim.zampar(:,1),     ethanol_clim.zampar(:,2),     shape.zampar, 'Color', color.zampar, 'MarkerFaceColor', color.zampar)
end 

plot(bart_WT_batch_min_t,   bart_WT_batch_min_y(:,num_y.eh),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_y(:,num_y.eh),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y(:,num_y.eh), 'Color', plt_clrs.red);
ylabel(label.met{num_met.eh})
xlabel(x_label);
ymax = max([max(bart_WT_batch_min_y(:,num_y.eh)), max(zamp_WT_batch_min_y(:,num_y.eh)), max(murphy_WT_batch_YPD_y(:,num_y.eh))]); 
if ymax > 0 
ylim([0 1.05.*ymax])
end 
xlim([0 bart_WT_batch_min_t(end)])
xlim([0 bart_WT_batch_min_t(end)])
%xline(t_gl_depl_bart, 'Color', plt_clrs.green)
%xline(t_eh_depl_bart, 'Color', plt_clrs.green)
%xline(t_gl_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_eh_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_gl_depl_murphy, 'Color', plt_clrs.red)
%xline(t_eh_depl_murphy, 'Color', plt_clrs.red)
axis square; 
box on;
hold off;
%}

% biomass
subplot(a, b, 1 + 3*b)
hold on;
if exist('glucose_clim','var')  
plot(cell_clim.bartolomeo(:,1), cell_clim.bartolomeo(:,2), shape.bartolomeo, 'Color', color.bartolomeo,   'MarkerFaceColor', color.bartolomeo)
plot(cell_clim.zampar(:,1),    cell_clim.zampar(:,2),      shape.zampar, 'Color', color.zampar,   'MarkerFaceColor', color.zampar)
plot(cell_clim.murphy(:,1),     cell_clim.murphy(:,2),     shape.murphy, 'Color', color.murphy, 'MarkerFaceColor', color.murphy)
end 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_y(:,end),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_y(:,end),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y(:,end), 'Color', plt_clrs.red);
ylabel(label.rib_cells{2})
xlabel(x_label);
ymax = max([max(bart_WT_batch_min_y(:,num_y.cells)), max(zamp_WT_batch_min_y(:,num_y.cells)), max(murphy_WT_batch_YPD_y(:,num_y.cells))]); 
if ymax > 0 
ylim([0 1.05.*ymax])
end 
xlim([0 bart_WT_batch_min_t(end)])
xlim([0 bart_WT_batch_min_t(end)])
%xline(t_gl_depl_bart, 'Color', plt_clrs.green)
%xline(t_eh_depl_bart, 'Color', plt_clrs.green)
%xline(t_gl_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_eh_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_gl_depl_murphy, 'Color', plt_clrs.red)
%xline(t_eh_depl_murphy, 'Color', plt_clrs.red)
axis square; 
box on;
hold off;

% CO2 
%{
co2 = bart_WT_batch_min_y(:,19).*(par.V_c/par.V_e) .* (met_reac.flux(:,2) + 3*met_reac.flux(:,4));
subplot(a, b, 2 + 3*b)
hold on;
plot(bartolomeo.co2(:,1), bartolomeo.co2(:,2), '.', 'Color', plt_clrs.blue)
plot(t, co2, 'Color', plt_clrs.blue);
ylabel('CO_{2}')
xlabel(x_label);
ylim([0 1.05.*max(co2)])
xlim([0 bart_WT_batch_min_t(end)])
%xline(t_gl_depl_bart, 'Color', plt_clrs.green)
%xline(t_eh_depl_bart, 'Color', plt_clrs.green)
%xline(t_gl_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_eh_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_gl_depl_murphy, 'Color', plt_clrs.red)
%xline(t_eh_depl_murphy, 'Color', plt_clrs.red)
axis square; 
box on;
hold off;
%}


%% virtual metabolites 

% intracellular glucose 
%{
subplot(a, b, 1 + 2*b)
hold on;
plot(solisescalante.g6p(:,1) + 10, solisescalante.g6p(:,2), 'o', 'Color', plt_clrs.green)
plot(zampar.g6p(:,1), zampar.g6p(:,2), 'o', 'Color', plt_clrs.blue)
plot(bart_WT_batch_min_t,   bart_WT_batch_min_y(:,num_y.gl_in),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_y(:,num_y.gl_in),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y(:,num_y.gl_in), 'Color', plt_clrs.red);
ylabel(label.met{num_met.gl_in})
xlabel(x_label);
ymax = max([max(bart_WT_batch_min_y(:,num_y.gl_in)), max(zamp_WT_batch_min_y(:,num_y.gl_in)), max(murphy_WT_batch_YPD_y(:,num_y.gl_in))]); 
if ymax > 0 
ylim([0 1.05.*ymax])
end 
xlim([0 bart_WT_batch_min_t(end)])
%xline(t_gl_depl_bart, 'Color', plt_clrs.green)
%xline(t_eh_depl_bart, 'Color', plt_clrs.green)
%xline(t_gl_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_eh_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_gl_depl_murphy, 'Color', plt_clrs.red)
%xline(t_eh_depl_murphy, 'Color', plt_clrs.red)
axis square; 
box on;
hold off;
%}

% precursor 
%
subplot(a, b, 2 + b)
hold on;
if exist('glucose_clim','var')  
plot(precursor_clim.f16bp.zampar(:,1), precursor_clim.f16bp.zampar(:,2),   'o-', 'Color',color.zampar, 'MarkerFaceColor', color.zampar)
plot(precursor_clim.f16bp.solis(:,1), precursor_clim.f16bp.solis(:,2), shape.solis, 'Color',color.solis, 'MarkerFaceColor', color.solis)
end 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_y(:,num_y.pc),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_y(:,num_y.pc),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y(:,num_y.pc), 'Color', plt_clrs.red);
ylabel(label.met{num_met.pc})
xlabel(x_label);
ymax = max([max(bart_WT_batch_min_y(:,num_y.pc)), max(zamp_WT_batch_min_y(:,num_y.pc)), max(murphy_WT_batch_YPD_y(:,num_y.pc))]); 
if ymax > 0 
ylim([0 1.05.*ymax])
end 
xlim([0 bart_WT_batch_min_t(end)])
%xline(t_gl_depl_bart, 'Color', plt_clrs.green)
%xline(t_eh_depl_bart, 'Color', plt_clrs.green)
%xline(t_gl_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_eh_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_gl_depl_murphy, 'Color', plt_clrs.red)
%xline(t_eh_depl_murphy, 'Color', plt_clrs.red)
axis square; 
box on;
hold off;
%}

% intracellular amino acids 
%
subplot(a, b, 2)
hold on;
if exist('glucose_clim','var')  
plot(aa_clim.glutamine.zampar(:,1), aa_clim.glutamine.zampar(:,2), 'o-', 'MarkerFaceColor', color.zampar);
end 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_y(:,num_y.aa_in),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_y(:,num_y.aa_in),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y(:,num_y.aa_in), 'Color', plt_clrs.red);
ylabel(label.met{num_met.aa_in})
xlabel(x_label);
ymax = max([max(bart_WT_batch_min_y(:,num_y.aa_in)), max(zamp_WT_batch_min_y(:,num_y.aa_in)), max(murphy_WT_batch_YPD_y(:,num_y.aa_in))]); 
if ymax > 0 
ylim([0 1.05.*ymax])
end
xlim([0 bart_WT_batch_min_t(end)])
%xline(t_gl_depl_bart, 'Color', plt_clrs.green)
%xline(t_eh_depl_bart, 'Color', plt_clrs.green)
%xline(t_gl_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_eh_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_gl_depl_murphy, 'Color', plt_clrs.red)
%xline(t_eh_depl_murphy, 'Color', plt_clrs.red)
axis square; 
box on;
hold off;
%}

% ATP 
%
subplot(a, b, 2 + 2*b)
hold on;
if exist('glucose_clim','var')  
plot(atp_clim.zampar(:,1),         atp_clim.zampar(:,2),    shape.zampar, 'Color', color.zampar, 'MarkerFaceColor', color.zampar)
end 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_y(:,num_y.ae),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_y(:,num_y.ae),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y(:,num_y.ae), 'Color', plt_clrs.red);
ylabel(label.met{num_met.ae})
xlabel(x_label);
ymax = max([max(bart_WT_batch_min_y(:,num_y.ae)), max(zamp_WT_batch_min_y(:,num_y.ae)), max(murphy_WT_batch_YPD_y(:,num_y.ae))]); 
if ymax > 0 
ylim([0 1.05.*ymax])
end 
xlim([0 bart_WT_batch_min_t(end)])
%xline(t_gl_depl_bart, 'Color', plt_clrs.green)
%xline(t_eh_depl_bart, 'Color', plt_clrs.green)
%xline(t_gl_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_eh_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_gl_depl_murphy, 'Color', plt_clrs.red)
%xline(t_eh_depl_murphy, 'Color', plt_clrs.red)
axis square; 
box on;
hold off;
%}

%% signaling 
subplot(a, b, 3)
hold on;
plot(bart_WT_batch_min_t,   bart_WT_batch_min_sig.snf1(:,1)   ./ par.s_tot, 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_sig.snf1(:,1)   ./ par.s_tot, 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_sig.snf1(:,1) ./ par.s_tot, 'Color', plt_clrs.red);
ylabel(label.sig.snf1)
xlabel(x_label);
ymax = max([max(bart_WT_batch_min_sig.snf1./par.s_tot), max(zamp_WT_batch_min_sig.snf1./par.s_tot), max(murphy_WT_batch_YPD_sig.snf1./par.s_tot)]); 
if ymax > 0 
ylim([0 1.05.*ymax])
end 
%xline(t_gl_depl_bart,   'Color', plt_clrs.green)
%xline(t_eh_depl_bart,   'Color', plt_clrs.green)
%xline(t_gl_depl_zamp,   'Color', plt_clrs.blue)
%xline(t_eh_depl_zamp,   'Color', plt_clrs.blue)
%xline(t_gl_depl_murphy, 'Color', plt_clrs.red)
%xline(t_eh_depl_murphy, 'Color', plt_clrs.red)
xlim([0 bart_WT_batch_min_t(end)])
yline(par.s_basal)
axis square; 
box on;
hold off; 

subplot(a, b, 3 + b)
hold on;
if exist('glucose_clim','var')
plot(garcia2014.glu_snf1_batch(:,1),    garcia2014.glu_snf1_batch(:,2),  shape.garcia2014, 'Color', color.garcia2014, 'MarkerFaceColor', color.garcia2014)
end 
plot(bart_WT_batch_min_y(:,num_y.gl),   bart_WT_batch_min_sig.snf1(:,1)   ./ par.s_tot, 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_y(:,num_y.gl),   zamp_WT_batch_min_sig.snf1(:,1)   ./ par.s_tot, 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_y(:,num_y.gl), murphy_WT_batch_YPD_sig.snf1(:,1) ./ par.s_tot, 'Color', plt_clrs.red);
ylabel(label.sig.snf1)
xlabel('glucose');
%{
ymax = max([max(bart_WT_batch_min_sig.snf1./par.s_tot), max(zamp_WT_batch_min_sig.snf1./par.s_tot), max(murphy_WT_batch_YPD_sig.snf1./par.s_tot)]); 
if ymax > 0 
ylim([0 1.05.*ymax])
end 
%}
ylim([0 1])
%xlim([0 bart_WT_batch_min_t(end)])
yline(1)
yline(par.s_basal)
axis square; 
box on;
hold off; 

subplot(a, b, 3 + 2*b)
hold on;
plot(bart_WT_batch_min_t,   bart_WT_batch_min_sig.tor(:,1) ./ par.tau_tot   , 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_sig.tor(:,1) ./ par.tau_tot   , 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_sig.tor(:,1)  ./ par.tau_tot, 'Color', plt_clrs.red);
ylabel(label.sig.tor)
xlabel(x_label);
ymax = max([max(bart_WT_batch_min_sig.tor  ./ par.tau_tot), max(zamp_WT_batch_min_sig.tor ./ par.tau_tot), max(murphy_WT_batch_YPD_sig.tor ./ par.tau_tot)]); 
if ymax > 0 
ylim([0 1.05.*ymax])
end 
%xline(t_gl_depl_bart,   'Color', plt_clrs.green)
%xline(t_eh_depl_bart,   'Color', plt_clrs.green)
%xline(t_gl_depl_zamp,   'Color', plt_clrs.blue)
%xline(t_eh_depl_zamp,   'Color', plt_clrs.blue)
%xline(t_gl_depl_murphy, 'Color', plt_clrs.red)
%xline(t_eh_depl_murphy, 'Color', plt_clrs.red)
xlim([0 bart_WT_batch_min_t(end)])
yline(par.tau_basal)
axis square; 
box on;
hold off; 

%{
subplot(a, b, 3 + 3*b)
hold on;
plot(bart_WT_batch_min_rib.ras,   bart_WT_batch_min_sig.tor(:,1)   ./ par.tau_tot, 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_rib.ras,   zamp_WT_batch_min_sig.tor(:,1)   ./ par.tau_tot, 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_rib.ras, murphy_WT_batch_YPD_sig.tor(:,1) ./ par.tau_tot, 'Color', plt_clrs.red);
ylabel(label.sig.tor)
xlabel('[R_{as}]');
%}
%{
ymax = max([max(bart_WT_batch_min_sig.tor ./ par.tau_tot), max(zamp_WT_batch_min_sig.tor./ par.tau_tot), max(murphy_WT_batch_YPD_sig.tor./ par.tau_tot)]); 
if ymax > 0 
ylim([0 1.05.*ymax])
end 
%}
%ylim([0 1])
%xlim([0 bart_WT_batch_min_t(end)])
%yline(1)
%yline(par.tau_basal)
axis square; 
box on;
hold off; 

%% total protein concentration
%
bart_WT_batch_min_total_protein_conc   = sum(bart_WT_batch_min_y(:,num_y.r:num_y.sd)  .* par.l',2);
zamp_WT_batch_min_total_protein_conc   = sum(zamp_WT_batch_min_y(:,num_y.r:num_y.sd)  .* par.l',2);
murphy_WT_batch_YPD_total_protein_conc = sum(murphy_WT_batch_YPD_y(:,num_y.r:num_y.sd).* par.l',2);
%{
%{
r_frac  = bart_WT_batch_min_y(:,num_y.r)*par.l(1)   ./ bart_WT_batch_min_total_protein_conc; 
z_frac  = bart_WT_batch_min_y(:,num_y.z)*par.l(2)   ./ bart_WT_batch_min_total_protein_conc; 
gy_frac = bart_WT_batch_min_y(:,num_y.gy)*par.l(3)  ./ bart_WT_batch_min_total_protein_conc; 
fe_frac = bart_WT_batch_min_y(:,num_y.fe)*par.l(4)  ./ bart_WT_batch_min_total_protein_conc; 
gn_frac = bart_WT_batch_min_y(:,num_y.gn)*par.l(5)  ./ bart_WT_batch_min_total_protein_conc; 
mt_frac = bart_WT_batch_min_y(:,num_y.mt)*par.l(6)  ./ bart_WT_batch_min_total_protein_conc; 
as_frac = bart_WT_batch_min_y(:,num_y.as)*par.l(7)  ./ bart_WT_batch_min_total_protein_conc; 
sp_frac = bart_WT_batch_min_y(:,num_y.sp)*par.l(8)  ./ bart_WT_batch_min_total_protein_conc; 
sd_frac = bart_WT_batch_min_y(:,num_y.sd)*par.l(9)  ./ bart_WT_batch_min_total_protein_conc; 
at_frac = bart_WT_batch_min_y(:,num_y.at)*par.l(10) ./ bart_WT_batch_min_total_protein_conc; 
ht_frac = bart_WT_batch_min_y(:,num_y.ht)*par.l(11) ./ bart_WT_batch_min_total_protein_conc; 
go_frac = bart_WT_batch_min_y(:,num_y.go)*par.l(12) ./ bart_WT_batch_min_total_protein_conc; 

max_z  = max(z_frac);
%}

subplot(a,b,3 + 4*b)
hold on; 

%yyaxis left
plot(bart_WT_batch_min_t,   bart_WT_batch_min_total_protein_conc,   '-', 'Color', plt_clrs.green)
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_total_protein_conc,   '-', 'Color', plt_clrs.blue)
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_total_protein_conc, '-', 'Color', plt_clrs.red)
ylabel('total protein concentration')
ylim([0 1.05*max(bart_WT_batch_min_total_protein_conc)])

%{
yyaxis right
plot(t, r_frac,  '-', 'Color', 'r')
plot(t, z_frac,  '-', 'Color', 'b')
plot(t, gy_frac, '-', 'Color', 'g')
plot(t, fe_frac, '-', 'Color', 'c')
plot(t, gn_frac, '-', 'Color', 'y')
plot(t, mt_frac, '-', 'Color', 'm')
plot(t, as_frac, '-', 'Color', [255 128 0]/255)
plot(t, sp_frac, '-', 'Color', [0 102 0]/255)
plot(t, sd_frac, '-', 'Color', [153 51 255]/255)
plot(t, at_frac, '-', 'Color', [160 160 160]/255)
plot(t, ht_frac, '-', 'Color', [153 0 76]/255)
plot(t, go_frac, '-', 'Color', [255 153 51]/255)
ylabel('protein fractions') 
ylim([0 1.05*max(max_z)])

plt = gca;
plt.YAxis(1).Color = 'k';
plt.YAxis(2).Color = 'k';
%}

%%xline(t_gl_depl)
%%xline(t_eh_depl)
%legend({'total', 'R', 'Z', 'gy', 'fe', 'gn', 'mt', 'as', 'sp', 'sd', 'at', 'ht', 'go'})
box on; 
axis square; 

hold off; 
%}

%% proteins 

% R 
subplot(a, b, 4)
hold on;
if exist('glucose_clim','var')  
plot(prot_sec_clim.time_prot.bartolomeo, prot_sec_clim.r.bartolomeo , shape.bartolomeo, 'Color', color.bartolomeo, 'MarkerFaceColor', color.bartolomeo)
plot(prot_sec_clim.time_prot.murphy, prot_sec_clim.r.murphy, shape.murphy, 'Color', color.murphy,      'MarkerFaceColor', color.murphy)
end 
plot(bart_WT_batch_min_t,   par.l(num_prot.r).*bart_WT_batch_min_y(:,num_y.r)   ./ bart_WT_batch_min_total_protein_conc,   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   par.l(num_prot.r).*zamp_WT_batch_min_y(:,num_y.r)   ./ zamp_WT_batch_min_total_protein_conc,   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, par.l(num_prot.r).*murphy_WT_batch_YPD_y(:,num_y.r) ./ murphy_WT_batch_YPD_total_protein_conc, 'Color', plt_clrs.red);
ylabel(label.prot{num_prot.r})
xlabel(x_label);
ymax1 = max(par.l(num_prot.r).*bart_WT_batch_min_y(:,num_y.r)   ./ bart_WT_batch_min_total_protein_conc);
ymax2 = max(par.l(num_prot.r).*zamp_WT_batch_min_y(:,num_y.r)   ./ zamp_WT_batch_min_total_protein_conc);
ymax3 = max(par.l(num_prot.r).*murphy_WT_batch_YPD_y(:,num_y.r) ./ murphy_WT_batch_YPD_total_protein_conc);
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
xlim([0 bart_WT_batch_min_t(end)])
%xline(t_gl_depl_bart, 'Color', plt_clrs.green)
%xline(t_eh_depl_bart, 'Color', plt_clrs.green)
%xline(t_gl_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_eh_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_gl_depl_murphy, 'Color', plt_clrs.red)
%xline(t_eh_depl_murphy, 'Color', plt_clrs.red)
axis square; 
box on;
hold off;

% Z 
subplot(a, b, 4 + b)
hold on;
if exist('glucose_clim','var') 
plot(prot_sec_clim.time_prot.bartolomeo, prot_sec_clim.z.bartolomeo , shape.bartolomeo, 'Color', color.bartolomeo, 'MarkerFaceColor', color.bartolomeo)
plot(prot_sec_clim.time_prot.murphy, prot_sec_clim.z.murphy, shape.murphy, 'Color', color.murphy,      'MarkerFaceColor', color.murphy)
end 
plot(bart_WT_batch_min_t,   par.l(num_prot.z).*bart_WT_batch_min_y(:,num_y.z)   ./ bart_WT_batch_min_total_protein_conc,   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   par.l(num_prot.z).*zamp_WT_batch_min_y(:,num_y.z)   ./ zamp_WT_batch_min_total_protein_conc,   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, par.l(num_prot.z).*murphy_WT_batch_YPD_y(:,num_y.z) ./ murphy_WT_batch_YPD_total_protein_conc, 'Color', plt_clrs.red);
ylabel(label.prot{num_prot.z})
xlabel(x_label);
ymax1 = max(par.l(num_prot.z).*bart_WT_batch_min_y(:,num_y.z)   ./ bart_WT_batch_min_total_protein_conc);
ymax2 = max(par.l(num_prot.z).*zamp_WT_batch_min_y(:,num_y.z)   ./ zamp_WT_batch_min_total_protein_conc);
ymax3 = max(par.l(num_prot.z).*murphy_WT_batch_YPD_y(:,num_y.z) ./ murphy_WT_batch_YPD_total_protein_conc);
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
xlim([0 bart_WT_batch_min_t(end)])
%xline(t_gl_depl_bart, 'Color', plt_clrs.green)
%xline(t_eh_depl_bart, 'Color', plt_clrs.green)
%xline(t_gl_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_eh_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_gl_depl_murphy, 'Color', plt_clrs.red)
%xline(t_eh_depl_murphy, 'Color', plt_clrs.red)
axis square; 
box on;
hold off;

% gy
subplot(a, b, 4 + 2*b)
hold on;
if exist('glucose_clim','var') 
plot(prot_sec_clim.time_prot.bartolomeo, prot_sec_clim.gy.bartolomeo , shape.bartolomeo, 'Color', color.bartolomeo, 'MarkerFaceColor', color.bartolomeo)
plot(prot_sec_clim.time_prot.murphy, prot_sec_clim.gy.murphy, shape.murphy, 'Color', color.murphy,      'MarkerFaceColor', color.murphy)
end 
plot(bart_WT_batch_min_t,   par.l(num_prot.gy).*bart_WT_batch_min_y(:,num_y.gy)   ./ bart_WT_batch_min_total_protein_conc,   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   par.l(num_prot.gy).*zamp_WT_batch_min_y(:,num_y.gy)   ./ zamp_WT_batch_min_total_protein_conc,   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, par.l(num_prot.gy).*murphy_WT_batch_YPD_y(:,num_y.gy) ./ murphy_WT_batch_YPD_total_protein_conc, 'Color', plt_clrs.red);
ylabel(label.prot{num_prot.gy})
xlabel(x_label);
ymax1 = max(par.l(num_prot.gy).*bart_WT_batch_min_y(:,num_y.gy)   ./ bart_WT_batch_min_total_protein_conc);
ymax2 = max(par.l(num_prot.gy).*zamp_WT_batch_min_y(:,num_y.gy)   ./ zamp_WT_batch_min_total_protein_conc);
ymax3 = max(par.l(num_prot.gy).*murphy_WT_batch_YPD_y(:,num_y.gy) ./ murphy_WT_batch_YPD_total_protein_conc);
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
xlim([0 bart_WT_batch_min_t(end)])
%xline(t_gl_depl_bart, 'Color', plt_clrs.green)
%xline(t_eh_depl_bart, 'Color', plt_clrs.green)
%xline(t_gl_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_eh_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_gl_depl_murphy, 'Color', plt_clrs.red)
%xline(t_eh_depl_murphy, 'Color', plt_clrs.red)
axis square; 
box on;
hold off;

% fe 
subplot(a, b, 4 + 3*b)
hold on;
if exist('glucose_clim','var') 
plot(prot_sec_clim.time_prot.bartolomeo, prot_sec_clim.fe.bartolomeo , shape.bartolomeo, 'Color', color.bartolomeo, 'MarkerFaceColor', color.bartolomeo)
plot(prot_sec_clim.time_prot.murphy, prot_sec_clim.fe.murphy, shape.murphy, 'Color', color.murphy,      'MarkerFaceColor', color.murphy)
end 
plot(bart_WT_batch_min_t,   par.l(num_prot.fe).*bart_WT_batch_min_y(:,num_y.fe)   ./ bart_WT_batch_min_total_protein_conc,   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   par.l(num_prot.fe).*zamp_WT_batch_min_y(:,num_y.fe)   ./ zamp_WT_batch_min_total_protein_conc,   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, par.l(num_prot.fe).*murphy_WT_batch_YPD_y(:,num_y.fe) ./ murphy_WT_batch_YPD_total_protein_conc, 'Color', plt_clrs.red);
ylabel(label.prot{num_prot.fe})
xlabel(x_label);
ymax1 = max(par.l(num_prot.fe).*bart_WT_batch_min_y(:,num_y.fe)   ./ bart_WT_batch_min_total_protein_conc);
ymax2 = max(par.l(num_prot.fe).*zamp_WT_batch_min_y(:,num_y.fe)   ./ zamp_WT_batch_min_total_protein_conc);
ymax3 = max(par.l(num_prot.fe).*murphy_WT_batch_YPD_y(:,num_y.fe) ./ murphy_WT_batch_YPD_total_protein_conc);
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
xlim([0 bart_WT_batch_min_t(end)])
%xline(t_gl_depl_bart, 'Color', plt_clrs.green)
%xline(t_eh_depl_bart, 'Color', plt_clrs.green)
%xline(t_gl_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_eh_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_gl_depl_murphy, 'Color', plt_clrs.red)
%xline(t_eh_depl_murphy, 'Color', plt_clrs.red)
axis square; 
box on;
hold off;

% gn 
subplot(a, b, 5 + 2*b)
hold on;
if exist('glucose_clim','var') 
plot(prot_sec_clim.time_prot.bartolomeo, prot_sec_clim.gn.bartolomeo , shape.bartolomeo, 'Color', color.bartolomeo, 'MarkerFaceColor', color.bartolomeo)
plot(prot_sec_clim.time_prot.murphy, prot_sec_clim.gn.murphy, shape.murphy, 'Color', color.murphy,      'MarkerFaceColor', color.murphy)
end 
plot(bart_WT_batch_min_t,   par.l(num_prot.gn).*bart_WT_batch_min_y(:,num_y.gn)   ./ bart_WT_batch_min_total_protein_conc,   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   par.l(num_prot.gn).*zamp_WT_batch_min_y(:,num_y.gn)   ./ zamp_WT_batch_min_total_protein_conc,   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, par.l(num_prot.gn).*murphy_WT_batch_YPD_y(:,num_y.gn) ./ murphy_WT_batch_YPD_total_protein_conc, 'Color', plt_clrs.red);
ylabel(label.prot{num_prot.gn})
xlabel(x_label);
ymax1 = max(par.l(num_prot.gn).*bart_WT_batch_min_y(:,num_y.gn)   ./ bart_WT_batch_min_total_protein_conc);
ymax2 = max(par.l(num_prot.gn).*zamp_WT_batch_min_y(:,num_y.gn)   ./ zamp_WT_batch_min_total_protein_conc);
ymax3 = max(par.l(num_prot.gn).*murphy_WT_batch_YPD_y(:,num_y.gn) ./ murphy_WT_batch_YPD_total_protein_conc);
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
xlim([0 bart_WT_batch_min_t(end)])
%xline(t_gl_depl_bart, 'Color', plt_clrs.green)
%xline(t_eh_depl_bart, 'Color', plt_clrs.green)
%xline(t_gl_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_eh_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_gl_depl_murphy, 'Color', plt_clrs.red)
%xline(t_eh_depl_murphy, 'Color', plt_clrs.red)
axis square; 
box on;
hold off;

% mt 
subplot(a, b, 5)
hold on;
if exist('glucose_clim','var') 
plot(prot_sec_clim.time_prot.bartolomeo, prot_sec_clim.mt.bartolomeo , shape.bartolomeo, 'Color', color.bartolomeo, 'MarkerFaceColor', color.bartolomeo)
plot(prot_sec_clim.time_prot.murphy, prot_sec_clim.mt.murphy, shape.murphy, 'Color', color.murphy,      'MarkerFaceColor', color.murphy)
end 
plot(bart_WT_batch_min_t,   par.l(num_prot.mt).*bart_WT_batch_min_y(:,num_y.mt)   ./ bart_WT_batch_min_total_protein_conc,   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   par.l(num_prot.mt).*zamp_WT_batch_min_y(:,num_y.mt)   ./ zamp_WT_batch_min_total_protein_conc,   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, par.l(num_prot.mt).*murphy_WT_batch_YPD_y(:,num_y.mt) ./ murphy_WT_batch_YPD_total_protein_conc, 'Color', plt_clrs.red);
ylabel(label.prot{num_prot.mt})
xlabel(x_label);
ylim([0 1.05*max(ylim)])
xlim([0 bart_WT_batch_min_t(end)])
%xline(t_gl_depl_bart, 'Color', plt_clrs.green)
%xline(t_eh_depl_bart, 'Color', plt_clrs.green)
%xline(t_gl_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_eh_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_gl_depl_murphy, 'Color', plt_clrs.red)
%xline(t_eh_depl_murphy, 'Color', plt_clrs.red)
axis square; 
box on;
hold off;

% as 
subplot(a, b, 5 + b)
hold on;
if exist('glucose_clim','var') 
plot(prot_sec_clim.time_prot.bartolomeo, prot_sec_clim.as.bartolomeo , shape.bartolomeo, 'Color', color.bartolomeo, 'MarkerFaceColor', color.bartolomeo)
plot(prot_sec_clim.time_prot.murphy, prot_sec_clim.as.murphy, shape.murphy, 'Color', color.murphy,      'MarkerFaceColor', color.murphy)
end 
plot(bart_WT_batch_min_t,   par.l(num_prot.as).*bart_WT_batch_min_y(:,num_y.as)   ./ bart_WT_batch_min_total_protein_conc,   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   par.l(num_prot.as).*zamp_WT_batch_min_y(:,num_y.as)   ./ zamp_WT_batch_min_total_protein_conc,   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, par.l(num_prot.as).*murphy_WT_batch_YPD_y(:,num_y.as) ./ murphy_WT_batch_YPD_total_protein_conc, 'Color', plt_clrs.red);
ylabel(label.prot{num_prot.as})
xlabel(x_label);
ylim([0 1.05*max(ylim)])
xlim([0 bart_WT_batch_min_t(end)])
%xline(t_gl_depl_bart, 'Color', plt_clrs.green)
%xline(t_eh_depl_bart, 'Color', plt_clrs.green)
%xline(t_gl_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_eh_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_gl_depl_murphy, 'Color', plt_clrs.red)
%xline(t_eh_depl_murphy, 'Color', plt_clrs.red)
axis square; 
box on;
hold off;

% at 
subplot(a, b, 5 + 3*b)
hold on;
if exist('glucose_clim','var') 
plot(prot_sec_clim.time_prot.bartolomeo, prot_sec_clim.at.bartolomeo , shape.bartolomeo, 'Color', color.bartolomeo, 'MarkerFaceColor', color.bartolomeo)
plot(prot_sec_clim.time_prot.murphy, prot_sec_clim.at.murphy, shape.murphy, 'Color', color.murphy,      'MarkerFaceColor', color.murphy)
end 
plot(bart_WT_batch_min_t,   par.l(num_prot.at).*bart_WT_batch_min_y(:,num_y.at)   ./ bart_WT_batch_min_total_protein_conc,   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   par.l(num_prot.at).*zamp_WT_batch_min_y(:,num_y.at)   ./ zamp_WT_batch_min_total_protein_conc,   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, par.l(num_prot.at).*murphy_WT_batch_YPD_y(:,num_y.at) ./ murphy_WT_batch_YPD_total_protein_conc, 'Color', plt_clrs.red);
ylabel(label.prot{num_prot.at})
xlabel(x_label);
ylim([0 1.05*max(ylim)])
xlim([0 bart_WT_batch_min_t(end)])
%xline(t_gl_depl_bart, 'Color', plt_clrs.green)
%xline(t_eh_depl_bart, 'Color', plt_clrs.green)
%xline(t_gl_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_eh_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_gl_depl_murphy, 'Color', plt_clrs.red)
%xline(t_eh_depl_murphy, 'Color', plt_clrs.red)
axis square; 
box on;
hold off;

%% ribosomes 
%figure; 

% total ribosomes
subplot(a, b, 6)

%yyaxis left
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_y(:,num_y.r0),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_y(:,num_y.r0),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y(:,num_y.r0), 'Color', plt_clrs.red);
%plot(bartolomeo.time_prot, bartolomeo.rib, 'o', 'Color', plt_clrs.blue,'HandleVisibility','off')
ylabel(label.rib_cells{1})
ymax1 = max(bart_WT_batch_min_y(:,num_y.r0));
ymax2 = max(zamp_WT_batch_min_y(:,num_y.r0));
ymax3 = max(murphy_WT_batch_YPD_y(:,num_y.r0));
ylim([0 1.05*max([ymax1, ymax2, ymax3])])

%{
yyaxis right 
plot(bart_WT_batch_min_t, rib.rf  ./ bart_WT_batch_min_y(:,num_y.r0), '-', 'Color', plt_clrs.red);   % active free ribosome fraction
plot(bart_WT_batch_min_t, rib.rat ./ bart_WT_batch_min_y(:,num_y.r0), '-', 'Color', plt_clrs.green); % active translating ribosome fraction
plot(bart_WT_batch_min_t, rib.ras ./ bart_WT_batch_min_y(:,num_y.r0), '-', 'Color', plt_clrs.gray);  % stalled translating ribosome fraction
ylabel('ribosome fractions')
ylim([0 1])

legend({'total', 'free', 'translating', 'stalled'})
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%}

%xline(t_gl_depl_bart, 'Color', plt_clrs.green)
%xline(t_eh_depl_bart, 'Color', plt_clrs.green)
%xline(t_gl_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_eh_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_gl_depl_murphy, 'Color', plt_clrs.red)
%xline(t_eh_depl_murphy, 'Color', plt_clrs.red)
xlim([0 bart_WT_batch_min_t(end)])
xlabel(x_label);

axis square; 
box on;
hold off; 

% bound ribosomes          
subplot(a, b, 6 + b)
hold on; 

%yyaxis left 
bart_WT_batch_min_transl_rib   = bart_WT_batch_min_y(:,num_y.r0)   - bart_WT_batch_min_rib.rf; 
zamp_WT_batch_min_transl_rib   = zamp_WT_batch_min_y(:,num_y.r0)   - zamp_WT_batch_min_rib.rf; 
murphy_WT_batch_YPD_transl_rib = murphy_WT_batch_YPD_y(:,num_y.r0) - murphy_WT_batch_YPD_rib.rf; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_transl_rib./bart_WT_batch_min_y(:,num_y.r0),   'Color', plt_clrs.green)
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_transl_rib./zamp_WT_batch_min_y(:,num_y.r0),   'Color', plt_clrs.blue)
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_transl_rib./murphy_WT_batch_YPD_y(:,num_y.r0), 'Color', plt_clrs.red)
%ymax1 = max(bart_WT_batch_min_transl_rib);
%ymax2 = max(zamp_WT_batch_min_transl_rib);
%ymax3 = max(murphy_WT_batch_YPD_transl_rib);
ylim([0 1])
ylabel('R_{b}/R_{0}')
yline(0.75)

%{
yyaxis right 
plot(bart_WT_batch_min_t, rib.rat./transl_rib, '-', 'Color', plt_clrs.red)
plot(bart_WT_batch_min_t, rib.ras./transl_rib, '-', 'Color', plt_clrs.green)
ylim([0 1.05])
ylabel('R_{at}, R_{as} fractions')

legend({'bound', 'active translating', 'stalled ribosomes'})
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%}

%xline(t_gl_depl_bart, 'Color', plt_clrs.green)
%xline(t_eh_depl_bart, 'Color', plt_clrs.green)
%xline(t_gl_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_eh_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_gl_depl_murphy, 'Color', plt_clrs.red)
%xline(t_eh_depl_murphy, 'Color', plt_clrs.red)
xlim([0 bart_WT_batch_min_t(end)])
xlabel(x_label)
axis square; 
box on;
hold off;

%% tRNA 

% total tRNA
subplot(a, b, 6 + 2*b)
hold on;

%{
%yyaxis left 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_prot_syn.t0,   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_prot_syn.t0,   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_prot_syn.t0, 'Color', plt_clrs.red);
ylabel('t_{0}')
ylim([0 1.05.*max(bart_WT_batch_min_prot_syn.t0)])
%}

%yyaxis right 
plot(bart_WT_batch_min_t, bart_WT_batch_min_prot_syn.tc./bart_WT_batch_min_prot_syn.t0, '-', 'Color', plt_clrs.green)   % charged tRNA fraction
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_prot_syn.tc./zamp_WT_batch_min_prot_syn.t0, '-', 'Color', plt_clrs.blue)
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_prot_syn.tc./murphy_WT_batch_YPD_prot_syn.t0, '-', 'Color', plt_clrs.red)
%plot(zamp_WT_batch_min_t, prot_syn.tu./prot_syn.t0, '-', 'Color', plt_clrs.green) % uncharged tRNA fraction
ylabel('t_{c}/t_{0}')
ylim([0 1.05])

%{
legend({'total', 'charged', 'uncharged'})
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%}

%xline(t_gl_depl_bart, 'Color', plt_clrs.green)
%xline(t_eh_depl_bart, 'Color', plt_clrs.green)
%xline(t_gl_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_eh_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_gl_depl_murphy, 'Color', plt_clrs.red)
%xline(t_eh_depl_murphy, 'Color', plt_clrs.red)
yline(0.8)
xlim([0 bart_WT_batch_min_t(end)])
xlabel(x_label);
axis square; 
box on;
hold off; 

%yyaxis right 

%{
% charged tRNA
subplot(a, b, 6 + b)
hold on;
plot(t, prot_syn.tc, 'Color', plt_clrs.blue);
ylabel('t_{c}')
xlabel(x_label);
ylim([0 1.05.*max(prot_syn.tc)])
%%xline(t_gl_depl)
%%xline(t_eh_depl)
axis square; 
box on;
hold off; 

% uncharged tRNA
subplot(a, b, 6 + 2*b)
hold on;
plot(t, prot_syn.tu, 'Color', plt_clrs.blue);
ylabel('t_{u}')
xlabel(x_label);
ylim([0 1.05.*max(prot_syn.tu)])
%%xline(t_gl_depl)
%%xline(t_eh_depl)
axis square; 
box on;
hold off; 
%}

%% growth rate

subplot(a, b, 2 + 3*b)
hold on;
plot(bart_WT_batch_min_t,   bart_WT_batch_min_g_rate,   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_g_rate,   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_g_rate, 'Color', plt_clrs.red);
ylabel('\lambda')
xlabel(x_label);
ymax1 = max(bart_WT_batch_min_g_rate);
ymax2 = max(zamp_WT_batch_min_g_rate);
ymax3 = max(murphy_WT_batch_YPD_g_rate);
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
%xline(t_gl_depl_bart, 'Color', plt_clrs.green)
%xline(t_eh_depl_bart, 'Color', plt_clrs.green)
%xline(t_gl_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_eh_depl_zamp, 'Color', plt_clrs.blue)
%xline(t_gl_depl_murphy, 'Color', plt_clrs.red)
%xline(t_eh_depl_murphy, 'Color', plt_clrs.red)
xlim([0 bart_WT_batch_min_t(end)])
axis square; 
box on;
hold off; 


%% newly introduced variables

% lipid
subplot(a, b, 1 + 4*b)
hold on;

plot(bart_WT_batch_min_t,   bart_WT_batch_min_y(:,num_y.lp_e),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_y(:,num_y.lp_e),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y(:,num_y.lp_e), 'Color', plt_clrs.red);
ylabel('lipid (uM)')
xlabel(x_label);
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim([0 bart_WT_batch_min_t(end)])
axis square; 
box on;
hold off; 


% NH4
subplot(a, b, 2 + 4*b)
hold on;
plot(bart_WT_batch_min_t,   bart_WT_batch_min_y(:,num_y.nh4),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_y(:,num_y.nh4),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y(:,num_y.nh4), 'Color', plt_clrs.red);
ylabel('nh4 (uM)')
xlabel(x_label);
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim([0 bart_WT_batch_min_t(end)])
axis square; 
box on;
hold off; 


% E_lp
subplot(a, b, 4 + 4*b)
hold on;
if exist('glucose_clim','var')
plot(prot_sec_clim.time_prot.bartolomeo, prot_sec_clim.lp.bartolomeo , shape.bartolomeo, 'Color', color.bartolomeo, 'MarkerFaceColor', color.bartolomeo)
plot(prot_sec_clim.time_prot.murphy, prot_sec_clim.lp.murphy, shape.murphy, 'Color', color.murphy,      'MarkerFaceColor', color.murphy)
end 
plot(bart_WT_batch_min_t,   par.l(num_prot.lp).*bart_WT_batch_min_y(:,num_y.lp)   ./ bart_WT_batch_min_total_protein_conc,   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   par.l(num_prot.lp).*zamp_WT_batch_min_y(:,num_y.lp)   ./ zamp_WT_batch_min_total_protein_conc,   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, par.l(num_prot.lp).*murphy_WT_batch_YPD_y(:,num_y.lp) ./ murphy_WT_batch_YPD_total_protein_conc, 'Color', plt_clrs.red);
ylabel(label.prot{num_prot.lp})
xlabel(x_label);
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim([0 bart_WT_batch_min_t(end)])
axis square; 
box on;
hold off; 


% E_lo
subplot(a, b, 5 + 4*b)
hold on;
if exist('glucose_clim','var')
plot(prot_sec_clim.time_prot.bartolomeo, prot_sec_clim.lo.bartolomeo , shape.bartolomeo, 'Color', color.bartolomeo, 'MarkerFaceColor', color.bartolomeo)
plot(prot_sec_clim.time_prot.murphy, prot_sec_clim.lo.murphy, shape.murphy, 'Color', color.murphy,      'MarkerFaceColor', color.murphy)
end 
plot(bart_WT_batch_min_t,   par.l(num_prot.lo).*bart_WT_batch_min_y(:,num_y.lo)   ./ bart_WT_batch_min_total_protein_conc,   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   par.l(num_prot.lo).*zamp_WT_batch_min_y(:,num_y.lo)   ./ zamp_WT_batch_min_total_protein_conc,   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, par.l(num_prot.lo).*murphy_WT_batch_YPD_y(:,num_y.lo) ./ murphy_WT_batch_YPD_total_protein_conc, 'Color', plt_clrs.red);
ylabel(label.prot{num_prot.lo})
xlabel(x_label);
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim([0 bart_WT_batch_min_t(end)])
axis square; 
box on;
hold off;

% E_sp
subplot(a, b, 6 + 3*b)
hold on;
if exist('glucose_clim','var')
plot(prot_sec_clim.time_prot.bartolomeo, prot_sec_clim.sp.bartolomeo , shape.bartolomeo, 'Color', color.bartolomeo, 'MarkerFaceColor', color.bartolomeo)
plot(prot_sec_clim.time_prot.murphy, prot_sec_clim.sp.murphy, shape.murphy, 'Color', color.murphy,      'MarkerFaceColor', color.murphy)
end 
plot(bart_WT_batch_min_t,   par.l(num_prot.sp).*bart_WT_batch_min_y(:,num_y.sp)   ./ bart_WT_batch_min_total_protein_conc,   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   par.l(num_prot.sp).*zamp_WT_batch_min_y(:,num_y.sp)   ./ zamp_WT_batch_min_total_protein_conc,   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, par.l(num_prot.sp).*murphy_WT_batch_YPD_y(:,num_y.sp) ./ murphy_WT_batch_YPD_total_protein_conc, 'Color', plt_clrs.red);
ylabel('sp')
xlabel(x_label);
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim([0 bart_WT_batch_min_t(end)])
axis square; 
box on;
hold off; 

% E_sd
subplot(a, b, 6 + 4*b)
hold on;
if exist('glucose_clim','var')
plot(prot_sec_clim.time_prot.bartolomeo, prot_sec_clim.sd.bartolomeo , shape.bartolomeo, 'Color', color.bartolomeo, 'MarkerFaceColor', color.bartolomeo)
plot(prot_sec_clim.time_prot.murphy, prot_sec_clim.sd.murphy, shape.murphy, 'Color', color.murphy,      'MarkerFaceColor', color.murphy)
end 
plot(bart_WT_batch_min_t,   par.l(num_prot.sd).*bart_WT_batch_min_y(:,num_y.sd)   ./ bart_WT_batch_min_total_protein_conc,   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   par.l(num_prot.sd).*zamp_WT_batch_min_y(:,num_y.sd)   ./ zamp_WT_batch_min_total_protein_conc,   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, par.l(num_prot.sd).*murphy_WT_batch_YPD_y(:,num_y.sd) ./ murphy_WT_batch_YPD_total_protein_conc, 'Color', plt_clrs.red);
ylabel('sd')
xlabel(x_label);
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim([0 bart_WT_batch_min_t(end)])
axis square; 
box on;
hold off; 


% Biomass (SC)
subplot(a, b, 3 + 3*b)
hold on;
plot(bart_WT_batch_min_t,   bart_WT_batch_min_y(:,num_y.sc) , 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_y(:,num_y.sc),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y(:,num_y.sc), 'Color', plt_clrs.red);
ylabel('storage carb (uM)')
xlabel(x_label);
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim([0 bart_WT_batch_min_t(end)])
axis square; 
box on;
hold off; 

% Biomass (Lipid)
subplot(a, b, 3 + 4*b)
hold on;
plot(bart_WT_batch_min_t,   bart_WT_batch_min_y(:,num_y.lp_e)  *lp_uM_to_gPerL , 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_y(:,num_y.lp_e)  *lp_uM_to_gPerL,   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y(:,num_y.lp_e)*lp_uM_to_gPerL, 'Color', plt_clrs.red);
ylabel('lipid (g/L)')
xlabel(x_label);
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim([0 bart_WT_batch_min_t(end)])
axis square; 
box on;
hold off; 

set(gcf,'Position', [440 86 1001 712])
sgtitle('batch fermentation with experiment')

end 
