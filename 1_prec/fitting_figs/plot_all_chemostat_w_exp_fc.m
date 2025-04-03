function plot_all_chemostat_w_exp_fc(D, ...
WT_y_steady_1gL, WT_y_steady_10gL, WT_y_steady_08gL, ...
WT_met_reac_steady_1gL, WT_met_reac_steady_10gL, WT_met_reac_steady_08gL, ...
WT_other_met_reac_steady_1gL,WT_other_met_reac_steady_10gL,WT_other_met_reac_steady_08gL,...
WT_sig_steady_1gL, WT_sig_steady_10gL, WT_sig_steady_08gL, ...
WT_prot_syn_steady_1gL, WT_prot_syn_steady_10gL, WT_prot_syn_steady_08gL, ...
WT_rib_steady_1gL, WT_rib_steady_10gL, WT_rib_steady_08gL, ...
WT_g_rate_steady_1gL, WT_g_rate_steady_10gL, WT_g_rate_steady_08gL, ...
par, label, plt_clrs, num_y, num_prot, num_flux, yeast_type, media_type)


%% load data
data_unit_conv; 
paper_figure_label; 
% [D_sc, cell_sc, glucose_sc, Jgy_sc, Jeh_sc,... 
% precursor_sc, atp_sc, aa_sc,  prot_sec_sc, ...
% lipid_sc, protein_sc,carbo_sc, RNA_sc, paperinfo_sc,...
% color_sc, shape_sc, gas_sc] = SC_C_lim_chem(); 
% 
% [D_rt, cell_rt, glucose_rt, Jgy_rt, nitrogen_rt, lipid_rt,...
%     protein_rt, carbo_rt, RNA_rt, mRNA_sec_rt, paperinfo_rt, ...
%     color_rt, shape_rt] ...
%     = Rt_C_lim_chem();


% media_type = 'nlim';  
% yeast_type = 'sc_rt';
switch media_type
case {'clim'} 
if strcmp(yeast_type,'sc') ...
|| strcmp(yeast_type,'sc_rt') ...         
[D_sc, cell_sc, glucose_sc, Jgy_sc, Jeh_sc,... 
precursor_sc, atp_sc, aa_sc,  prot_sec_sc, ...
lipid_sc, protein_sc,carbo_sc, RNA_sc, paperinfo_sc,...
color_sc, shape_sc, gas_sc] = SC_C_lim_chem(); 
end 
if strcmp(yeast_type,'rt') ...
|| strcmp(yeast_type,'sc_rt') ...                
[D_rt, cell_rt, glucose_rt, Jgy_rt, nh4_rt, lipid_rt,...
protein_rt, carbo_rt, RNA_rt, mRNA_sec_rt, paperinfo_rt, ...
color_rt, shape_rt] ...
= Rt_C_lim_chem();
[~, ~, ~, ~, ~, ~, ~, ~,  ~, ~, ~,~, ~, ~,...
color_sc, shape_sc] = SC_C_lim_chem(); 
end 
case {'nlim'} 
if strcmp(yeast_type,'sc') ...
|| strcmp(yeast_type,'sc_rt') ...    
[D_sc, cell_sc, glucose_sc, Jgy_sc, Jeh_sc, nh4_sc, gas_sc, ... 
precursor_sc, atp_sc, aa_sc,  prot_sec_sc, ...
lipid_sc, protein_sc,carbo_sc, RNA_sc, paper_info_sc,...                    
color_sc, shape_sc, Jnh4_sc]  = SC_N_lim_chem(); 
end 
if strcmp(yeast_type,'rt') ...
|| strcmp(yeast_type,'sc_rt') ...              
[D_rt, cell_rt, glucose_rt, Jgy_rt, nh4_rt, lipid_rt, protein_rt,...
carbo_rt,RNA_rt, mRNA_sec_rt, paper_info_rt, ...
color_rt, shape_rt] = Rt_N_lim_chem();  
[~, ~, ~, ~, ~, ~, ~, ... 
~, ~, ~,  ~, ...
~, ~,~, ~, ~,...                    
color_sc, shape_sc]  = SC_N_lim_chem();             
end 
otherwise 
disp('media_type should be either clim or nlim')
end 

color_sc.default =  color_sc.kumar_hg;
color_rt.default  = color_sc.hackett; 
if strcmp(yeast_type,'sc_rt')
% sc green rt blue
% overwrite colors and shape
color_sc.default    = plt_clrs.green;
color_sc.kumar_hg   = plt_clrs.green;
color_sc.xia        = plt_clrs.green;
color_sc.hackett    = plt_clrs.green; 
color_sc.kumar_lg   = plt_clrs.green;
color_sc.boer       = plt_clrs.green;

shape_sc.kumar_hg   = 'o'; 
shape_sc.xia        = '^'; 

shape_sc.hackett    = 'o'; 
shape_sc.kumar_lg   = 'd'; 
shape_sc.boer       = 's';


shape_sc.lange     = 'v'; 
color_sc.lange     = plt_clrs.green;
shape_sc.ertugay   = 'v'; 
color_sc.ertugay   = plt_clrs.green;
shape_sc.Yu4       = 'v'; 
color_sc.Yu4       = plt_clrs.green;

% rt
color_rt.default    = plt_clrs.blue;
color_rt.shen       = plt_clrs.blue;
shape_rt.shen       = 'o'; 
color_rt.shen2017   = plt_clrs.blue;
shape_rt.shen2017   = '^';  
color_rt.Dinh       = plt_clrs.blue; 
shape_rt.Dinh       = 'o'; 
end

% overwrite colors and shape
% sc green rt blue
% overwrite colors and shape
color_sc.default    = plt_clrs.green;
color_sc.kumar_hg   = plt_clrs.green;
color_sc.xia        = plt_clrs.green;
color_sc.hackett    = plt_clrs.green; 
color_sc.kumar_lg   = plt_clrs.green;
color_sc.boer       = plt_clrs.green;

shape_sc.kumar_hg   = 'o'; 
shape_sc.xia        = '^'; 

shape_sc.hackett    = 'o'; 
shape_sc.kumar_lg   = 'd'; 
shape_sc.boer       = 's';


shape_sc.lange     = 'v'; 
color_sc.lange     = plt_clrs.green;
shape_sc.ertugay   = 'v'; 
color_sc.ertugay   = plt_clrs.green;
shape_sc.Yu4       = 'v'; 
color_sc.Yu4       = plt_clrs.green;

% rt
color_rt.default    = plt_clrs.blue;
color_rt.shen       = plt_clrs.blue;
shape_rt.shen       = 'o'; 
color_rt.shen2017   = plt_clrs.blue;
shape_rt.shen2017   = '^';  
color_rt.Dinh       = plt_clrs.blue; 
shape_rt.Dinh       = 'o'; 


norm_indx = 1; 
total_protein_con_10gL = sum(WT_y_steady_10gL(:,num_y.r:num_y.sd).* par.l(num_prot.r:num_prot.sd)',2);
total_protein_con_1gL  = sum(WT_y_steady_1gL(:,num_y.r:num_y.sd).* par.l(num_prot.r:num_prot.sd)',2);
total_protein_con_08gL = sum(WT_y_steady_08gL(:,num_y.r:num_y.sd).* par.l(num_prot.r:num_prot.sd)',2);
% proteins frac (sim)
WT_y_steady_10gL(:,[num_y.r: num_y.sd]) = (WT_y_steady_10gL(:,[num_y.r: num_y.sd]).*par.l')./ total_protein_con_10gL;
WT_y_steady_1gL(:,[num_y.r: num_y.sd])  = (WT_y_steady_1gL(:,[num_y.r: num_y.sd]).*par.l')./ total_protein_con_1gL;
WT_y_steady_08gL(:,[num_y.r: num_y.sd]) = (WT_y_steady_08gL(:,[num_y.r: num_y.sd]).*par.l')./ total_protein_con_08gL;

% prot_sec_sc.r.hackett   = prot_sec_sc.r.hackett./ prot_sec_sc.r.hackett(1);
% prot_sec_sc.r.xia       = prot_sec_sc.r.xia./ prot_sec_sc.r.xia(1);
% prot_sec_sc.z.hackett   = prot_sec_sc.z.hackett./ prot_sec_sc.z.hackett(1);
% prot_sec_sc.z.xia       = prot_sec_sc.z.xia./prot_sec_sc.z.xia(1);
% prot_sec_sc.e.hackett   = prot_sec_sc.e.hackett./ prot_sec_sc.e.hackett(1);
% prot_sec_sc.e.xia       = prot_sec_sc.e.xia./prot_sec_sc.e.xia(1);
% prot_sec_sc.gy.hackett  = prot_sec_sc.gy.hackett./ prot_sec_sc.gy.hackett(1);
% prot_sec_sc.gy.xia      = prot_sec_sc.gy.xia./prot_sec_sc.gy.xia(1);
% prot_sec_sc.fe.hackett  = prot_sec_sc.fe.hackett./ prot_sec_sc.fe.hackett(1);
% prot_sec_sc.fe.xia      = prot_sec_sc.fe.xia./prot_sec_sc.fe.xia(1);
% prot_sec_sc.gn.hackett  = prot_sec_sc.gn.hackett./ prot_sec_sc.gn.hackett(1);
% prot_sec_sc.gn.xia      = prot_sec_sc.gn.xia./prot_sec_sc.gn.xia(1);
% prot_sec_sc.mt.hackett  = prot_sec_sc.mt.hackett./ prot_sec_sc.mt.hackett(1);
% prot_sec_sc.mt.xia      = prot_sec_sc.mt.xia./prot_sec_sc.mt.xia(1);
% prot_sec_sc.sp.hackett  = prot_sec_sc.sp.hackett./ prot_sec_sc.sp.hackett(1);
% prot_sec_sc.sp.xia      = prot_sec_sc.sp.xia./prot_sec_sc.sp.xia(1);
% prot_sec_sc.sd.hackett  = prot_sec_sc.sd.hackett./ prot_sec_sc.sd.hackett(1);
% prot_sec_sc.sd.xia      = prot_sec_sc.sd.xia./prot_sec_sc.sd.xia(1);
% prot_sec_sc.as.hackett  = prot_sec_sc.as.hackett./ prot_sec_sc.as.hackett(1);
% prot_sec_sc.as.xia      = prot_sec_sc.as.xia./prot_sec_sc.as.xia(1);
% prot_sec_sc.at.hackett  = prot_sec_sc.at.hackett./ prot_sec_sc.at.hackett(1);
% prot_sec_sc.at.xia      = prot_sec_sc.at.xia./prot_sec_sc.at.xia(1);
% prot_sec_sc.lp.hackett  = prot_sec_sc.lp.hackett./ prot_sec_sc.lp.hackett(1);
%prot_sec_sc.lp.xia      = prot_sec_sc.lp.xia./prot_sec_sc.lp.xia(1);

x_label ='D (h^{-1})'; 
plot_num; 

%% plotting stuff

D_crit_simu = 0.17; 
D_crit_simu_indx = find(D-D_crit_simu>0,1); 

%% real metabolites
a = 5;
b = 7;
figure; 

%{
for k = 1:num.met
subplot(a, b, 1 + (b*(k-1)))
hold on;
plot(D, WT_y_steady_1gL(:,k), '-', 'Color', color_sc.default); 
ylabel(label.met{k})
xlabel(x_label_chem);
if 1.05*max(WT_y_steady_1gL(:,k)) > 0
ylim([0 1.1.*max(WT_y_steady_1gL(:,k))]) 
end
%xline(d_crit)
axis square; 
box on;
hold off; 
end 
%}

% extracellular amino acids 
subplot(a, b, 1)
hold on; 
plot(D, WT_y_steady_1gL(:,num_y.aaex), '-', 'Color', color_sc.default);
plot(D, WT_y_steady_10gL(:,num_y.aaex), '-', 'Color', color_sc.kumar_hg);
plot(D, WT_y_steady_08gL(:,num_y.aaex), '-', 'Color', color_rt.default);
ylabel(label.met{num_met.aaex})
xlabel(x_label);
ymax1 = max(WT_y_steady_1gL(:,num_y.aaex)); 
ymax2 = max(WT_y_steady_10gL(:,num_y.aaex)); 
ymax3 = max(WT_y_steady_08gL(:,num_y.aaex)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end 
%xline(d_crit)
xlim([D(1)-0.05 D(end)+0.05])
axis square; 
box on;
hold off;

% glucose
subplot(a, b, 1 + b)
hold on;
plot(D, WT_y_steady_1gL(:,num_y.gl), '-', 'Color', color_sc.default);
plot(D, WT_y_steady_10gL(:,num_y.gl), '-', 'Color', color_sc.kumar_hg);
plot(D, WT_y_steady_08gL(:,num_y.gl), '-', 'Color', color_rt.default);
ylabel(label.met{num_met.gl})
xlabel(x_label);
ymax1 = max(WT_y_steady_1gL(:,num_y.gl)); 
ymax2 = max(WT_y_steady_10gL(:,num_y.gl)); 
ymax3 = max(WT_y_steady_08gL(:,num_y.gl)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end %xline(d_crit)
xlim([D(1)-0.05 D(end)+0.05])
axis square; 
box on;
hold off;

% J_gy
subplot(a, b, 2 + b); 
hold on; 
papers = fieldnames(Jgy_sc);
for i = 1: length(papers)
plot(D_sc.(papers{i}),  Jgy_sc.(papers{i}) ,...
shape_sc.(papers{i}),...
'Color', color_sc.(papers{i}), ...
'MarkerFaceColor', color_sc.(papers{i}))
end    
papers = fieldnames(Jgy_rt); 
for i = 1: length(papers)
plot(D_rt.(papers{i}),  Jgy_rt.(papers{i}) ,...
shape_rt.(papers{i}),...
'Color', color_rt.(papers{i}), ...
'MarkerFaceColor', color_rt.(papers{i}))
end 
plot(D, WT_met_reac_steady_1gL.flux(:,num_flux.gy), '-', 'Color', color_sc.default);
plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.gy), '-', 'Color', color_sc.kumar_hg);
plot(D, WT_met_reac_steady_08gL.flux(:,num_flux.gy), '-', 'Color', color_rt.default);
ylabel('J_{gy}')    
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.flux(:,num_flux.gy)); 
ymax2 = max(WT_met_reac_steady_10gL.flux(:,num_flux.gy)); 
ymax3 = max(WT_met_reac_steady_08gL.flux(:,num_flux.gy)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end 
%xline(d_crit)
axis square; 
box on; 
hold off; 

% ethanol 
%
subplot(a, b, 1 + 3*b)
hold on;
plot(D, WT_y_steady_1gL(:,num_y.eh), '-', 'Color', color_sc.default);
plot(D, WT_y_steady_10gL(:,num_y.eh), '-', 'Color', color_sc.kumar_hg);
plot(D, WT_y_steady_08gL(:,num_y.eh), '-', 'Color', color_rt.default);
ylabel(label.met{num_met.eh})
xlabel(x_label);
ymax1 = max(WT_y_steady_1gL(:,num_y.eh)); 
ymax2 = max(WT_y_steady_10gL(:,num_y.eh)); 
ymax3 = max(WT_y_steady_08gL(:,num_y.eh)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end 
%xline(d_crit)
xlim([D(1)-0.05 D(end)+0.05])
axis square; 
box on;
hold off;
%}

% J_fe - J_gn
subplot(a, b, 2 + 3*b); 
hold on;  
papers = fieldnames(Jeh_sc);
for i = 1: length(papers)
plot(D_sc.(papers{i}),  Jeh_sc.(papers{i}) ,...
shape_sc.(papers{i}),...
'Color', color_sc.(papers{i}), ...
'MarkerFaceColor', color_sc.(papers{i}))
end  

plot(D, WT_met_reac_steady_1gL.flux(:,num_flux.fe) - WT_met_reac_steady_1gL.flux(:,num_flux.gn), '-', 'Color',  color_sc.kumar_lg);
plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.fe) - WT_met_reac_steady_10gL.flux(:,num_flux.gn), '-', 'Color',  color_sc.kumar_hg);
plot(D, WT_met_reac_steady_08gL.flux(:,num_flux.fe) - WT_met_reac_steady_08gL.flux(:,num_flux.gn), '-', 'Color',  color_sc.hackett);
ylabel('J_{fe} - J_{gn}')    
ymax1 = max(WT_met_reac_steady_1gL.flux(:,num_flux.fe) - WT_met_reac_steady_1gL.flux(:,num_flux.gn)); 
ymax2 = max(WT_met_reac_steady_10gL.flux(:,num_flux.fe) - WT_met_reac_steady_10gL.flux(:,num_flux.gn)); 
ymax3 = max(WT_met_reac_steady_08gL.flux(:,num_flux.fe) - WT_met_reac_steady_08gL.flux(:,num_flux.gn)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end 
xlabel(x_label)
%xline(d_crit)
xlim([D(1)-0.05 D(end)+0.05])
axis square; 
box on; 
hold off; 

% biomass
subplot(a, b, 1 + 4*b)
hold on;
papers = fieldnames(cell_sc); 
for i = 1: length(papers)
plot(D_sc.(papers{i}),  cell_sc.(papers{i}) ,...
shape_sc.(papers{i}),...
'Color', color_sc.(papers{i}), ...
'MarkerFaceColor', color_sc.(papers{i}))
end 
papers = fieldnames(cell_rt); 
for i = 1: length(papers)
plot(D_rt.(papers{i}),  cell_rt.(papers{i}) ,...
shape_rt.(papers{i}),...
'Color', color_rt.(papers{i}), ...
'MarkerFaceColor', color_rt.(papers{i}))
end 
plot(D, WT_y_steady_1gL(:,end), '-', 'Color', color_sc.default);
plot(D, WT_y_steady_10gL(:,end), '-', 'Color', color_sc.kumar_hg);
plot(D, WT_y_steady_08gL(:,end), '-', 'Color', color_rt.default);
ylabel(label.rib_cells{2})
xlabel(x_label);
if max(ylim) > 0 
ylim([0, 1.05*max(ylim)])
end 
%xline(d_crit)
xlim([D(1)-0.05 D(end)+0.05])
axis square; 
box on;
hold off;

% biomass yield
% gl_in_1gL  = 1.0  * gl_gPerL_to_uM; 
% gl_in_10gL = 10.0 * gl_gPerL_to_uM; 
% gl_in_08gL = 0.8  * gl_gPerL_to_uM; 
% biomass_yield_1gL  = (WT_y_steady_1gL(:,end)  .* (1/gdw_to_cell))./((gl_in_1gL*(ones(size(WT_y_steady_1gL(:,num_y.gl),1),1))) - ((WT_y_steady_1gL(:,num_y.gl) .* 1/gl_gPerL_to_uM)));
% biomass_yield_10gL = (WT_y_steady_10gL(:,end) .* (1/gdw_to_cell))./((gl_in_10gL*(ones(size(WT_y_steady_10gL(:,num_y.gl),1),1))) - ((WT_y_steady_10gL(:,num_y.gl) .* 1/gl_gPerL_to_uM)));
% biomass_yield_08gL = (WT_y_steady_08gL(:,end) .* (1/gdw_to_cell))./((gl_in_08gL*(ones(size(WT_y_steady_08gL(:,num_y.gl),1),1))) - ((WT_y_steady_08gL(:,num_y.gl) .* 1/gl_gPerL_to_uM)));
% subplot(a, b, 2 + 4*b)
% hold on; 
% plot(D, biomass_yield_1gL, '-', 'Color', color_sc.default);
% plot(D, biomass_yield_10gL, '-', 'Color', color_sc.kumar_hg);
% plot(D, biomass_yield_08gL, '-', 'Color', color_rt.default);
% ylabel({'biomass yield'; '(gdw cells/ g glucose'; 'consumed)'})
% xlabel(x_label);
% ymax1 = max(biomass_yield_1gL); 
% ymax2 = max(biomass_yield_10gL); 
% ymax3 = max(biomass_yield_08gL); 
% if max([ymax1, ymax2, ymax3]) > 0 
%     ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
% end 
% %xline(d_crit)
% xlim([D(1)-0.05 D(end)+0.05])
% axis square; 
% box on; 
% hold off; 


% CO2 production rate
co2_1gL = (WT_met_reac_steady_1gL.flux(:,num_flux.fe) + 3*WT_met_reac_steady_1gL.flux(:,num_flux.mt));
co2_10gL = (WT_met_reac_steady_10gL.flux(:,num_flux.fe) + 3*WT_met_reac_steady_10gL.flux(:,num_flux.mt));
co2_08gL = (WT_met_reac_steady_08gL.flux(:,num_flux.fe) + 3*WT_met_reac_steady_08gL.flux(:,num_flux.mt));
subplot(a, b, 3 + 3*b)
hold on;
papers = fieldnames(gas_sc.co2);
for i = 1: length(papers)
plot(D_sc.(papers{i}),  gas_sc.co2.(papers{i}) ,...
shape_sc.(papers{i}),...
'Color', color_sc.(papers{i}), ...
'MarkerFaceColor', color_sc.(papers{i}))
end 

plot(D, co2_1gL, '-', 'Color', color_sc.default);
plot(D, co2_10gL, '-', 'Color', color_sc.kumar_hg);
plot(D, co2_08gL, '-', 'Color', color_rt.default);
ylabel('CO_{2}')
xlabel(x_label);
ymax1 = max(co2_1gL); 
ymax2 = max(co2_10gL); 
ymax3 = max(co2_08gL);
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end 
%xline(d_crit)
xlim([D(1)-0.05 D(end)+0.05])
axis square; 
box on;
hold off;
%}

%% virtual metabolites 



% precursor 
%
subplot(a, b, 2 + 2*b)
hold on;
precs = fieldnames(precursor_sc);
for k = 1: length(precs)
prec = precs{k};
papers = fieldnames(precursor_sc.(prec)); 
for i = 1: length(papers)
plot(D_sc.(papers{i}),  precursor_sc.(prec).(papers{i})./precursor_sc.(prec).(papers{i})(1) ,...
shape_sc.(papers{i}),...
'Color', color_sc.(papers{i}), ...
'MarkerFaceColor', color_sc.(papers{i}))
end 
end 
plot(D, WT_y_steady_1gL(:,num_y.pc), '-', 'Color', color_sc.default);
plot(D, WT_y_steady_10gL(:,num_y.pc), '-', 'Color', color_sc.kumar_hg);
plot(D, WT_y_steady_08gL(:,num_y.pc), '-', 'Color', color_rt.default);
ylabel(label.met{num_met.pc})
xlabel(x_label);
ymax1 = max(WT_y_steady_1gL(:,num_y.pc)); 
ymax2 = max(WT_y_steady_10gL(:,num_y.pc)); 
ymax3 = max(WT_y_steady_08gL(:,num_y.pc)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end 
%xline(d_crit)
xlim([D(1)-0.05 D(end)+0.05])
axis square; 
box on;
hold off;
%}

% intracellular amino acids 
%
subplot(a, b, 3)
hold on;
papers = fieldnames(aa_sc.glutamine); 
for i = 1: length(papers)
plot(D_sc.(papers{i}),  aa_sc.glutamine.(papers{i})./aa_sc.glutamine.(papers{i})(1) ,...
shape_sc.(papers{i}),...
'Color', color_sc.(papers{i}), ...
'MarkerFaceColor', color_sc.(papers{i}))
end 
plot(D, WT_y_steady_1gL(:,num_y.aa_in), '-', 'Color', color_sc.default);
plot(D, WT_y_steady_10gL(:,num_y.aa_in), '-', 'Color', color_sc.kumar_hg);
plot(D, WT_y_steady_08gL(:,num_y.aa_in), '-', 'Color', color_rt.default);
ylabel(label.met{num_met.aa_in})
xlabel(x_label);
ymax1 = max(WT_y_steady_1gL(:,num_y.aa_in)); 
ymax2 = max(WT_y_steady_10gL(:,num_y.aa_in)); 
ymax3 = max(WT_y_steady_08gL(:,num_y.aa_in)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end 
%xline(d_crit)
xlim([D(1)-0.05 D(end)+0.05])
axis square; 
box on;
hold off;
%}

%{
% Jas = lambda*[Aa] + rho*lambda
kumar_lg_Jas = (kumar_lg.dr .* kumar_lg.glutamine) + (par.pro_den.*kumar_lg.dr); 
kumar_hg_Jas = (kumar_hg.dr .* kumar_hg.glutamine) + (par.pro_den.*kumar_hg.dr); 
%
subplot(a, b, 2)
hold on;
plot(kumar_lg.dr, kumar_lg_Jas, 'o', 'MarkerFaceColor', color_sc.kumar_lg);
plot(kumar_hg.dr, kumar_hg_Jas, 'o', 'MarkerFaceColor', color_sc.kumar_hg);
plot(D, WT_met_reac_steady_1gL.flux(:,num_flux.as), '-', 'Color', color_sc.default); 
plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.as), '-', 'Color', color_sc.kumar_hg); 
ylabel('J_{as}')
xlabel(x_label);
ymax1 = max(WT_met_reac_steady_1gL.flux(:,num_flux.as)); 
ymax2 = max(WT_met_reac_steady_10gL.flux(:,num_flux.as)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end 
%xline(d_crit)
xlim([D(1)-0.05 D(end)+0.05])
axis square; 
box on;
hold off;
%}

% ATP 
%
subplot(a, b, 3 + 2*b)
hold on;
papers = fieldnames(atp_sc); 
for i = 1: length(papers)
plot(D_sc.(papers{i}),  atp_sc.(papers{i}) ,...
shape_sc.(papers{i}),...
'Color', color_sc.(papers{i}), ...
'MarkerFaceColor', color_sc.(papers{i}))
end  

plot(D, WT_y_steady_1gL(:,num_y.ae), '-', 'Color', color_sc.default);
plot(D, WT_y_steady_10gL(:,num_y.ae), '-', 'Color', color_sc.kumar_hg);
plot(D, WT_y_steady_08gL(:,num_y.ae), '-', 'Color', color_rt.default);
ylabel(label.met(num_met.ae))
xlabel(x_label);
ymax1 = max(WT_y_steady_1gL(:,num_y.ae)); 
ymax2 = max(WT_y_steady_10gL(:,num_y.ae)); 
ymax3 = max(WT_y_steady_08gL(:,num_y.ae));
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end 
%xline(d_crit)
xlim([D(1)-0.05 D(end)+0.05])
axis square; 
box on;
hold off;

% NH4 
%
subplot(a, b, 1 + 2*b)
hold on;
% plot(kumar_lg.dr, kumar_lg.nh4, 'o', 'MarkerFaceColor', color_sc.kumar_lg)
% plot(kumar_hg.dr, kumar_hg.nh4, 'o', 'MarkerFaceColor', color_sc.kumar_hg)
% plot(hackett.dr,  hackett.nh4,  'o', 'MarkerFaceColor', color_sc.hackett)
plot(D, WT_y_steady_1gL(:,num_y.nh4), '-', 'Color', color_sc.default);
plot(D, WT_y_steady_10gL(:,num_y.nh4), '-', 'Color', color_sc.kumar_hg);
plot(D, WT_y_steady_08gL(:,num_y.nh4), '-', 'Color', color_rt.default);
ylabel(label.met(num_met.nh4))
xlabel(x_label);
if max(ylim) > 0 
ylim([0 1.1.*max(ylim)])
end 
%xline(d_crit)
xlim([D(1)-0.05 D(end)+0.05])
axis square; 
box on;
hold off;


% lipid 
%
subplot(a, b, 3 + 1*b)
hold on;
% plot(kumar_lg.dr, kumar_lg.nh4, 'o', 'MarkerFaceColor', color_sc.kumar_lg)
% plot(kumar_hg.dr, kumar_hg.nh4, 'o', 'MarkerFaceColor', color_sc.kumar_hg)
% plot(hackett.dr,  hackett.nh4,  'o', 'MarkerFaceColor', color_sc.hackett)
plot(D, WT_y_steady_1gL(:,num_y.lp_e), '-', 'Color', color_sc.default);
plot(D, WT_y_steady_10gL(:,num_y.lp_e), '-', 'Color', color_sc.kumar_hg);
plot(D, WT_y_steady_08gL(:,num_y.lp_e), '-', 'Color', color_rt.default);
ylabel('lipid (uM)')
xlabel(x_label);
if max(ylim) > 0 
ylim([0 1.1.*max(ylim)])
end 
%xline(d_crit)
xlim([D(1)-0.05 D(end)+0.05])
axis square; 
box on;
hold off;
%}

% Jmt
%Jmt_1gL = (1/3) .* kumar_lg.co2_prod;

subplot(a, b, 4 + 2*b)
hold on; 
%plot(kumar_lg.dr, Jmt_1gL, 'o', 'MarkerFaceColor', color_sc.kumar_lg)
% if exist('gas_sc','var')
%     Jmt_10gL = (1/3) .* gas_sc.o2.kumar_hg;
%     plot(D_sc.kumar_hg, Jmt_10gL, 'o', 'MarkerFaceColor', color_sc.kumar_hg)
% end
plot(D, WT_met_reac_steady_1gL.flux(:,num_flux.mt), '-', 'Color', color_sc.default);
plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.mt), '-', 'Color', color_sc.kumar_hg);
plot(D, WT_met_reac_steady_08gL.flux(:,num_flux.mt), '-', 'Color', color_rt.default);
ylabel('J_{mt}')
xlabel(x_label);
ymax1 = max(WT_met_reac_steady_1gL.flux(:,num_flux.mt)); 
ymax2 = max(WT_met_reac_steady_10gL.flux(:,num_flux.mt)); 
ymax3 = max(WT_met_reac_steady_08gL.flux(:,num_flux.mt)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end 
%xline(d_crit)
xlim([D(1)-0.05 D(end)+0.05])
axis square; 
box on;
hold off;

%% signaling 
subplot(a, b, 4)
hold on;
plot(D, WT_sig_steady_1gL.snf1(:,1), '-', 'Color', color_sc.default);
plot(D, WT_sig_steady_10gL.snf1(:,1), '-', 'Color', color_sc.kumar_hg);
plot(D, WT_sig_steady_08gL.snf1(:,1), '-', 'Color', color_rt.default);
ylabel(label.sig.snf1)
xlabel(x_label);
ymax1 = max(WT_sig_steady_1gL.snf1(:,1)); 
ymax2 = max(WT_sig_steady_10gL.snf1(:,1)); 
ymax3 = max(WT_sig_steady_08gL.snf1(:,1)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end 
%xline(d_crit)
yline(par.s_basal)
xlim([D(1)-0.05 D(end)+0.05])
axis square; 
box on;
hold off; 

subplot(a, b, 4 + b)
hold on;
plot(D, WT_sig_steady_1gL.tor(:,1) ./ par.tau_tot, '-', 'Color', color_sc.default);
plot(D, WT_sig_steady_10gL.tor(:,1) ./ par.tau_tot, '-', 'Color', color_sc.kumar_hg);
plot(D, WT_sig_steady_08gL.tor(:,1) ./ par.tau_tot, '-', 'Color', color_rt.default);
ylabel(label.sig.tor)
xlabel(x_label);
ymax1 = max(WT_sig_steady_1gL.tor(:,1) ./ par.tau_tot); 
ymax2 = max(WT_sig_steady_10gL.tor(:,1) ./ par.tau_tot); 
ymax3 = max(WT_sig_steady_08gL.tor(:,1) ./ par.tau_tot); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end 
%xline(d_crit)
yline(par.tau_basal)
xlim([D(1)-0.05 D(end)+0.05])
axis square; 
box on;
hold off; 

%% total protein concentration
%
% total_protein_con_1gL = sum(WT_y_steady_1gL(:,num_y.r:num_y.lo).*par.l',2);
% total_protein_con_10gL = sum(WT_y_steady_1gL(:,num_y.r:num_y.lo).*par.l',2);
% total_protein_con_08gL = sum(WT_y_steady_08gL(:,num_y.r:num_y.lo).*par.l',2);

%{
yyaxis left
plot(D, r_frac,  '-', 'Color', 'r')
plot(D, z_frac,  '-', 'Color', 'b')
plot(D, gy_frac, '-', 'Color', 'g')
plot(D, fe_frac, '-', 'Color', 'c')
plot(D, gn_frac, '-', 'Color', 'y')
plot(D, mt_frac, '-', 'Color', 'm')
plot(D, as_frac, '-', 'Color', [255 128 0]/255)
plot(D, sp_frac, '-', 'Color', [0 102 0]/255)
plot(D, sd_frac, '-', 'Color', [153 51 255]/255)
plot(D, at_frac, '-', 'Color', [160 160 160]/255)
plot(D, ht_frac, '-', 'Color', [153 0 76]/255)
plot(D, go_frac, '-', 'Color', [255 153 51]/255)
ylabel('protein fractions') 
%ylim([0 1.05*max(max_z)])
ylim([0 1])


yyaxis right
%}
hold on 
plot(D, total_protein_con_1gL, '-', 'Color', color_sc.default)
plot(D, total_protein_con_10gL, '-', 'Color', color_sc.kumar_hg)
plot(D, total_protein_con_08gL, '-', 'Color', color_rt.default)
ylabel('total protein concentration')
ymax1 = max(total_protein_con_1gL); 
ymax2 = max(total_protein_con_10gL); 
ymax3 = max(total_protein_con_08gL); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end 

%{
plt = gca;
plt.YAxis(1).Color = 'k';
plt.YAxis(2).Color = 'k';
%}

%xline(d_crit)
xlim([D(1)-0.05 D(end)+0.05])
%legend({'total', 'R', 'Z', 'gy', 'fe', 'gn', 'mt', 'as', 'sp', 'sd', 'at', 'ht', 'go'})
box on; 
axis square; 

hold off; 
%}
%% proteins 
%figure; 
%{
for k = 1:num.prot
subplot(a, b, 3 + (b*(k-1)))
hold on;
plot(D, WT_y_steady_1gL(:,k + num.met), '-', 'Color', color_sc.default);
ylabel(label.prot{k})
xlabel(x_label_chem);
ylim([0 1.1.*max(WT_y_steady_1gL(:,k + num.met))])
%xline(d_crit)
axis square; 
box on;
hold off; 
end 
%}


% R 
subplot(a, b, 5)
hold on;
papers = fieldnames(prot_sec_sc.r); 
for i = 1: length(papers)
plot(D_sc.(papers{i}),  prot_sec_sc.r.(papers{i})./prot_sec_sc.gy.(papers{i})(1) ,...
shape_sc.(papers{i}),...
'Color', color_sc.(papers{i}), ...
'MarkerFaceColor', color_sc.(papers{i}))
end 
%     papers = fieldnames(prot_sec_rt.r); 
%     for i = 1: length(papers)
%         plot(D_rt.(papers{i}),  prot_sec_rt.r.(papers{i})./prot_sec_sc.r.(papers{i})(1) ,...
%             shape_rt.(papers{i}),...
%             'Color', color_rt.(papers{i}), ...
%             'MarkerFaceColor', color_rt.(papers{i}))
%     end  

plot(D, WT_y_steady_1gL(:,num_y.r)  ./ WT_y_steady_1gL(1,num_y.r), '-', 'Color', color_sc.default);
plot(D, WT_y_steady_10gL(:,num_y.r) ./ WT_y_steady_10gL(1,num_y.r), '-', 'Color', color_sc.kumar_hg);
plot(D, WT_y_steady_08gL(:,num_y.r) ./ WT_y_steady_08gL(1,num_y.r), '-', 'Color', color_rt.default);
ylabel(label.prot{num_prot.r})
xlabel(x_label);
ymax1 = max(WT_y_steady_1gL(:,num_y.r)  ./ WT_y_steady_1gL(1,num_y.r)); 
ymax2 = max(WT_y_steady_10gL(:,num_y.r) ./ WT_y_steady_10gL(1,num_y.r)); 
ymax3 = max(WT_y_steady_08gL(:,num_y.r) ./ WT_y_steady_08gL(1,num_y.r)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end 
xlim([D(1)-0.05 D(end)+0.05])
%xline(d_crit)
axis square; 
box on;
hold off;

% Z 
subplot(a, b, 5 + b)
hold on;
papers = fieldnames(prot_sec_sc.z); 
for i = 1: length(papers)
plot(D_sc.(papers{i}),  prot_sec_sc.z.(papers{i})./prot_sec_sc.z.(papers{i})(1) ,...
shape_sc.(papers{i}),...
'Color', color_sc.(papers{i}), ...
'MarkerFaceColor', color_sc.(papers{i}))
end 
%     papers = fieldnames(prot_sec_rt.z); 
%     for i = 1: length(papers)
%         plot(D_rt.(papers{i}),  prot_sec_rt.z.(papers{i})./prot_sec_sc.z.(papers{i})(1) ,...
%             shape_rt.(papers{i}),...
%             'Color', color_rt.(papers{i}), ...
%             'MarkerFaceColor', color_rt.(papers{i}))
%     end  

plot(D, WT_y_steady_1gL(:,num_y.z)  ./ WT_y_steady_1gL(1,num_y.z), '-', 'Color', color_sc.default);
plot(D, WT_y_steady_10gL(:,num_y.z) ./ WT_y_steady_10gL(1,num_y.z), '-', 'Color', color_sc.kumar_hg);
plot(D, WT_y_steady_08gL(:,num_y.z) ./ WT_y_steady_08gL(1,num_y.z), '-', 'Color', color_rt.default);
ylabel(label.prot{num_prot.z})
xlabel(x_label);
ymax1 = max(WT_y_steady_1gL(:,num_y.z)  ./ WT_y_steady_1gL(1,num_y.z)); 
ymax2 = max(WT_y_steady_10gL(:,num_y.z) ./ WT_y_steady_10gL(1,num_y.z)); 
ymax3 = max(WT_y_steady_08gL(:,num_y.z) ./ WT_y_steady_08gL(1,num_y.z)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end
xlim([D(1)-0.05 D(end)+0.05])
%xline(d_crit)
axis square; 
box on;
hold off;

% gy
subplot(a, b, 5 + 2*b)
hold on;
papers = fieldnames(prot_sec_sc.gy); 
for i = 1: length(papers)
plot(D_sc.(papers{i}),  prot_sec_sc.gy.(papers{i})./prot_sec_sc.gy.(papers{i})(1) ,...
shape_sc.(papers{i}),...
'Color', color_sc.(papers{i}), ...
'MarkerFaceColor', color_sc.(papers{i}))
end 
%     papers = fieldnames(prot_sec_rt.gy); 
%     for i = 1: length(papers)
%         plot(D_rt.(papers{i}),  prot_sec_rt.gy.(papers{i})./prot_sec_sc.gy.(papers{i})(1) ,...
%             shape_rt.(papers{i}),...
%             'Color', color_rt.(papers{i}), ...
%             'MarkerFaceColor', color_rt.(papers{i}))
%     end  

plot(D, WT_y_steady_1gL(:,num_y.gy)  ./ WT_y_steady_1gL(1,num_y.gy), '-', 'Color', color_sc.default);
plot(D, WT_y_steady_10gL(:,num_y.gy) ./ WT_y_steady_10gL(1,num_y.gy), '-', 'Color', color_sc.kumar_hg);
plot(D, WT_y_steady_08gL(:,num_y.gy) ./ WT_y_steady_08gL(1,num_y.gy), '-', 'Color', color_rt.default);
ylabel(label.prot{num_prot.gy})
xlabel(x_label);
ymax1 = max(WT_y_steady_1gL(:,num_y.gy)  ./ WT_y_steady_1gL(1,num_y.gy)); 
ymax2 = max(WT_y_steady_10gL(:,num_y.gy) ./ WT_y_steady_10gL(1,num_y.gy)); 
ymax3 = max(WT_y_steady_08gL(:,num_y.gy) ./ WT_y_steady_08gL(1,num_y.gy)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end
xlim([D(1)-0.05 D(end)+0.05])
%xline(d_crit)
axis square; 
box on;
hold off;

% fe 
subplot(a, b, 5 + 3*b)
hold on;
papers = fieldnames(prot_sec_sc.fe); 
for i = 1: length(papers)
plot(D_sc.(papers{i}),  prot_sec_sc.fe.(papers{i})./prot_sec_sc.fe.(papers{i})(1) ,...
shape_sc.(papers{i}),...
'Color', color_sc.(papers{i}), ...
'MarkerFaceColor', color_sc.(papers{i}))
end 
%     papers = fieldnames(prot_sec_rt.fe); 
%     for i = 1: length(papers)
%         plot(D_rt.(papers{i}),  prot_sec_rt.fe.(papers{i})./prot_sec_sc.fe.(papers{i})(1) ,...
%             shape_rt.(papers{i}),...
%             'Color', color_rt.(papers{i}), ...
%             'MarkerFaceColor', color_rt.(papers{i}))
%     end  
plot(D, WT_y_steady_1gL(:,num_y.fe)  ./ WT_y_steady_1gL(1,num_y.fe), '-', 'Color', color_sc.default);
plot(D, WT_y_steady_10gL(:,num_y.fe) ./ WT_y_steady_10gL(1,num_y.fe), '-', 'Color', color_sc.kumar_hg);
plot(D, WT_y_steady_08gL(:,num_y.fe) ./ WT_y_steady_08gL(1,num_y.fe), '-', 'Color', color_rt.default);
ylabel(label.prot{num_prot.fe})
xlabel(x_label);
ymax1 = max(WT_y_steady_1gL(:,num_y.fe)  ./ WT_y_steady_1gL(1,num_y.fe)); 
ymax2 = max(WT_y_steady_10gL(:,num_y.fe) ./ WT_y_steady_10gL(1,num_y.fe)); 
ymax3 = max(WT_y_steady_08gL(:,num_y.fe) ./ WT_y_steady_08gL(1,num_y.fe)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end
xlim([D(1)-0.05 D(end)+0.05])
%xline(d_crit)
axis square; 
box on;
hold off;

% gn 
subplot(a, b, 5 + 4*b)
hold on;
papers = fieldnames(prot_sec_sc.gn); 
for i = 1: length(papers)
plot(D_sc.(papers{i}),  prot_sec_sc.gn.(papers{i})./prot_sec_sc.gn.(papers{i})(1) ,...
shape_sc.(papers{i}),...
'Color', color_sc.(papers{i}), ...
'MarkerFaceColor', color_sc.(papers{i}))
end 
%     papers = fieldnames(prot_sec_rt.gn); 
%     for i = 1: length(papers)
%         plot(D_rt.(papers{i}),  prot_sec_rt.gn.(papers{i})./prot_sec_sc.gn.(papers{i})(1) ,...
%             shape_rt.(papers{i}),...
%             'Color', color_rt.(papers{i}), ...
%             'MarkerFaceColor', color_rt.(papers{i}))
%     end  
plot(D, WT_y_steady_1gL(:,num_y.gn)  ./ WT_y_steady_1gL(1,num_y.gn), '-', 'Color', color_sc.default);
plot(D, WT_y_steady_10gL(:,num_y.gn) ./ WT_y_steady_10gL(1,num_y.gn), '-', 'Color', color_sc.kumar_hg);
plot(D, WT_y_steady_08gL(:,num_y.gn) ./ WT_y_steady_08gL(1,num_y.gn), '-', 'Color', color_rt.default);
ylabel(label.prot{num_prot.gn})
xlabel(x_label);
ymax1 = max(WT_y_steady_1gL(:,num_y.gn)  ./ WT_y_steady_1gL(1,num_y.gn)); 
ymax2 = max(WT_y_steady_10gL(:,num_y.gn) ./ WT_y_steady_10gL(1,num_y.gn)); 
ymax3 = max(WT_y_steady_08gL(:,num_y.gn) ./ WT_y_steady_08gL(1,num_y.gn)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end
xlim([D(1)-0.05 D(end)+0.05])
%xline(d_crit)
axis square; 
box on;
hold off;

% mt 
subplot(a, b, 6)
hold on;
papers = fieldnames(prot_sec_sc.mt); 
for i = 1: length(papers)
plot(D_sc.(papers{i}),  prot_sec_sc.mt.(papers{i})./prot_sec_sc.mt.(papers{i})(1) ,...
shape_sc.(papers{i}),...
'Color', color_sc.(papers{i}), ...
'MarkerFaceColor', color_sc.(papers{i}))
end 
%     papers = fieldnames(prot_sec_rt.mt); 
%     for i = 1: length(papers)
%         plot(D_rt.(papers{i}),  prot_sec_rt.mt.(papers{i})./prot_sec_sc.mt.(papers{i})(1) ,...
%             shape_rt.(papers{i}),...
%             'Color', color_rt.(papers{i}), ...
%             'MarkerFaceColor', color_rt.(papers{i}))
%     end  
plot(D, WT_y_steady_1gL(:,num_y.mt)  ./ WT_y_steady_1gL(1,num_y.mt), '-', 'Color', color_sc.default);
plot(D, WT_y_steady_10gL(:,num_y.mt) ./ WT_y_steady_10gL(1,num_y.mt), '-', 'Color', color_sc.kumar_hg);
plot(D, WT_y_steady_08gL(:,num_y.mt) ./ WT_y_steady_08gL(1,num_y.mt), '-', 'Color', color_rt.default);
ylabel(label.prot{num_prot.mt})
xlabel(x_label);
ymax1 = max(WT_y_steady_1gL(:,num_y.mt)  ./ WT_y_steady_1gL(1,num_y.mt)); 
ymax2 = max(WT_y_steady_10gL(:,num_y.mt) ./ WT_y_steady_10gL(1,num_y.mt)); 
ymax3 = max(WT_y_steady_08gL(:,num_y.mt) ./ WT_y_steady_08gL(1,num_y.mt)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end
xlim([D(1)-0.05 D(end)+0.05])
%xline(d_crit)
axis square; 
box on;
hold off;

% as 
subplot(a, b, 6 + b)
hold on;
papers = fieldnames(prot_sec_sc.as); 
for i = 1: length(papers)
plot(D_sc.(papers{i}),  prot_sec_sc.as.(papers{i})./prot_sec_sc.as.(papers{i})(1) ,...
shape_sc.(papers{i}),...
'Color', color_sc.(papers{i}), ...
'MarkerFaceColor', color_sc.(papers{i}))
end 
%     papers = fieldnames(prot_sec_rt.as); 
%     for i = 1: length(papers)
%         plot(D_rt.(papers{i}),  prot_sec_rt.as.(papers{i})./prot_sec_sc.as.(papers{i})(1) ,...
%             shape_rt.(papers{i}),...
%             'Color', color_rt.(papers{i}), ...
%             'MarkerFaceColor', color_rt.(papers{i}))
%     end  
plot(D, WT_y_steady_1gL(:,num_y.as)  ./ WT_y_steady_1gL(1,num_y.as), '-', 'Color', color_sc.default);
plot(D, WT_y_steady_10gL(:,num_y.as) ./ WT_y_steady_10gL(1,num_y.as), '-', 'Color', color_sc.kumar_hg);
plot(D, WT_y_steady_08gL(:,num_y.as) ./ WT_y_steady_08gL(1,num_y.as), '-', 'Color', color_rt.default);
ylabel(label.prot{num_prot.as})
xlabel(x_label);
ymax1 = max(WT_y_steady_1gL(:,num_y.as)  ./ WT_y_steady_1gL(1,num_y.as)); 
ymax2 = max(WT_y_steady_10gL(:,num_y.as) ./ WT_y_steady_10gL(1,num_y.as)); 
ymax3 = max(WT_y_steady_08gL(:,num_y.as) ./ WT_y_steady_08gL(1,num_y.as)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end
xlim([D(1)-0.05 D(end)+0.05])
%xline(d_crit)
axis square; 
box on;
hold off;

% at 
subplot(a, b, 6 + 4*b)
hold on;
papers = fieldnames(prot_sec_sc.at); 
for i = 1: length(papers)
plot(D_sc.(papers{i}),  prot_sec_sc.at.(papers{i})./prot_sec_sc.at.(papers{i})(1) ,...
shape_sc.(papers{i}),...
'Color', color_sc.(papers{i}), ...
'MarkerFaceColor', color_sc.(papers{i}))
end 
%     papers = fieldnames(prot_sec_rt.at); 
%     for i = 1: length(papers)
%         plot(D_rt.(papers{i}),  prot_sec_rt.at.(papers{i})./prot_sec_sc.at.(papers{i})(1) ,...
%             shape_rt.(papers{i}),...
%             'Color', color_rt.(papers{i}), ...
%             'MarkerFaceColor', color_rt.(papers{i}))
%     end  
plot(D, WT_y_steady_1gL(:,num_y.at)  ./ WT_y_steady_1gL(1,num_y.at), '-', 'Color', color_sc.default);
plot(D, WT_y_steady_10gL(:,num_y.at) ./ WT_y_steady_10gL(1,num_y.at), '-', 'Color', color_sc.kumar_hg);
plot(D, WT_y_steady_08gL(:,num_y.at) ./ WT_y_steady_08gL(1,num_y.at), '-', 'Color', color_rt.default);
ylabel(label.prot{num_prot.at})
xlabel(x_label);
ymax1 = max(WT_y_steady_1gL(:,num_y.at)  ./ WT_y_steady_1gL(1,num_y.at)); 
ymax2 = max(WT_y_steady_10gL(:,num_y.at) ./ WT_y_steady_10gL(1,num_y.at)); 
ymax3 = max(WT_y_steady_08gL(:,num_y.at) ./ WT_y_steady_08gL(1,num_y.at)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end 
xlim([D(1)-0.05 D(end)+0.05])
%xline(d_crit)
axis square; 
box on;
hold off;


% lp
subplot(a, b, 6 +2*b)
hold on;
papers = fieldnames(prot_sec_sc.lp); 
for i = 1: length(papers)
plot(D_sc.(papers{i}),  prot_sec_sc.lp.(papers{i})./prot_sec_sc.lp.(papers{i})(1) ,...
shape_sc.(papers{i}),...
'Color', color_sc.(papers{i}), ...
'MarkerFaceColor', color_sc.(papers{i}))
end 
%     papers = fieldnames(prot_sec_rt.lp); 
%     for i = 1: length(papers)
%         plot(D_rt.(papers{i}),  prot_sec_rt.lp.(papers{i})./prot_sec_sc.lp.(papers{i})(1) ,...
%             shape_rt.(papers{i}),...
%             'Color', color_rt.(papers{i}), ...
%             'MarkerFaceColor', color_rt.(papers{i}))
%     end  
plot(D, WT_y_steady_1gL(:,num_y.lp)  ./ WT_y_steady_1gL(1,num_y.lp), '-', 'Color', color_sc.default);
plot(D, WT_y_steady_10gL(:,num_y.lp) ./ WT_y_steady_10gL(1,num_y.lp), '-', 'Color', color_sc.kumar_hg);
plot(D, WT_y_steady_08gL(:,num_y.lp) ./ WT_y_steady_08gL(1,num_y.lp), '-', 'Color', color_rt.default);
ylabel(label.prot{num_prot.lp})
xlabel(x_label);
if max(ylim) > 0 
ylim([0 1.1.*max(ylim)])
end 
xlim([D(1)-0.05 D(end)+0.05])
%xline(d_crit)
axis square; 
box on;
hold off;



% lo 
subplot(a, b, 6 + 3*b)
hold on;
papers = fieldnames(prot_sec_sc.lo); 
for i = 1: length(papers)
plot(D_sc.(papers{i}),  prot_sec_sc.lo.(papers{i})./prot_sec_sc.lo.(papers{i})(1) ,...
shape_sc.(papers{i}),...
'Color', color_sc.(papers{i}), ...
'MarkerFaceColor', color_sc.(papers{i}))
end 
%     papers = fieldnames(prot_sec_rt.lo); 
%     for i = 1: length(papers)
%         plot(D_rt.(papers{i}),  prot_sec_rt.lo.(papers{i})./prot_sec_sc.lo.(papers{i})(1) ,...
%             shape_rt.(papers{i}),...
%             'Color', color_rt.(papers{i}), ...
%             'MarkerFaceColor', color_rt.(papers{i}))
%     end 

plot(D, WT_y_steady_1gL(:,num_y.lo)  ./ WT_y_steady_1gL(1,num_y.lo), '-', 'Color', color_sc.default);
plot(D, WT_y_steady_10gL(:,num_y.lo) ./ WT_y_steady_10gL(1,num_y.lo), '-', 'Color', color_sc.kumar_hg);
plot(D, WT_y_steady_08gL(:,num_y.lo) ./ WT_y_steady_08gL(1,num_y.lo), '-', 'Color', color_rt.default);
ylabel(label.prot{num_prot.lo})
xlabel(x_label);
if max(ylim) > 0 
ylim([0 1.1.*max(ylim)])
end
xlim([D(1)-0.05 D(end)+0.05])
%xline(d_crit)
axis square; 
box on;
hold off;


% sp
subplot(a, b, 7 +3*b)
hold on;
papers = fieldnames(prot_sec_sc.sp); 
for i = 1: length(papers)
plot(D_sc.(papers{i}),  prot_sec_sc.sp.(papers{i})./prot_sec_sc.sp.(papers{i})(1) ,...
shape_sc.(papers{i}),...
'Color', color_sc.(papers{i}), ...
'MarkerFaceColor', color_sc.(papers{i}))
end 
%     papers = fieldnames(prot_sec_rt.lp); 
%     for i = 1: length(papers)
%         plot(D_rt.(papers{i}),  prot_sec_rt.lp.(papers{i})./prot_sec_sc.lp.(papers{i})(1) ,...
%             shape_rt.(papers{i}),...
%             'Color', color_rt.(papers{i}), ...
%             'MarkerFaceColor', color_rt.(papers{i}))
%     end  
plot(D, WT_y_steady_1gL(:,num_y.sp)  ./ WT_y_steady_1gL(1,num_y.sp), '-', 'Color', color_sc.default);
plot(D, WT_y_steady_10gL(:,num_y.sp) ./ WT_y_steady_10gL(1,num_y.sp), '-', 'Color', color_sc.kumar_hg);
plot(D, WT_y_steady_08gL(:,num_y.sp) ./ WT_y_steady_08gL(1,num_y.sp), '-', 'Color', color_rt.default);
ylabel('sp')
xlabel(x_label);
if max(ylim) > 0 
ylim([0 1.1.*max(ylim)])
end 
xlim([D(1)-0.05 D(end)+0.05])
%xline(d_crit)
axis square; 
box on;
hold off;



% sd 
subplot(a, b, 7 + 4*b)
hold on;
papers = fieldnames(prot_sec_sc.sd); 
for i = 1: length(papers)
plot(D_sc.(papers{i}),  prot_sec_sc.sd.(papers{i})./prot_sec_sc.sd.(papers{i})(1) ,...
shape_sc.(papers{i}),...
'Color', color_sc.(papers{i}), ...
'MarkerFaceColor', color_sc.(papers{i}))
end 
%     papers = fieldnames(prot_sec_rt.sd); 
%     for i = 1: length(papers)
%         plot(D_rt.(papers{i}),  prot_sec_rt.sd.(papers{i})./prot_sec_sc.sd.(papers{i})(1) ,...
%             shape_rt.(papers{i}),...
%             'Color', color_rt.(papers{i}), ...
%             'MarkerFaceColor', color_rt.(papers{i}))
%     end 

plot(D, WT_y_steady_1gL(:,num_y.sd)  ./ WT_y_steady_1gL(1,num_y.sd), '-', 'Color', color_sc.default);
plot(D, WT_y_steady_10gL(:,num_y.sd) ./ WT_y_steady_10gL(1,num_y.sd), '-', 'Color', color_sc.kumar_hg);
plot(D, WT_y_steady_08gL(:,num_y.sd) ./ WT_y_steady_08gL(1,num_y.sd), '-', 'Color', color_rt.default);
ylabel('sd')
xlabel(x_label);
if max(ylim) > 0 
ylim([0 1.1.*max(ylim)])
end
xlim([D(1)-0.05 D(end)+0.05])
%xline(d_crit)
axis square; 
box on;
hold off;


%% proteins produced 
%{
beta_tot_1gL = ones(length(D),1);
beta_tot_10gL = ones(length(D),1);
for k = 1:length(D)
beta_tot_1gL(k) = sum(WT_prot_syn_steady_1gL.beta(k,:).*par.l'); 
beta_tot_10gL(k) = sum(WT_prot_syn_steady_10gL.beta(k,:).*par.l'); 
end 

subplot(a, b, 4 + 3*b)
hold on; 
plot(D, beta_tot_1gL, '-', 'Color', color_sc.default)
plot(D, beta_tot_10gL, '-', 'Color', color_sc.kumar_hg)
ylabel({'proteins produced', '(uM aa h^{-1})'})
xlabel(x_label);
ymax1 = max(beta_tot_1gL); 
ymax2 = max(beta_tot_10gL); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end 
xlim([D(1)-0.05 D(end)+0.05])
%xline(d_crit)
axis square; 
box on;
hold off;
%}

%% ribosomes 
%figure; 

% total ribosomes
subplot(a, b, 7)
hold on;
%yyaxis left
plot(D, WT_y_steady_1gL(:,num_y.r0), '-', 'Color', color_sc.default);
plot(D, WT_y_steady_10gL(:,num_y.r0), '-', 'Color', color_sc.kumar_hg);
plot(D, WT_y_steady_08gL(:,num_y.r0), '-', 'Color', color_rt.default);
ylabel(label.rib_cells{1})
ymax1 = max(WT_y_steady_1gL(:,num_y.r0)); 
ymax2 = max(WT_y_steady_10gL(:,num_y.r0)); 
ymax3 = max(WT_y_steady_08gL(:,num_y.r0)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end 
%{
yyaxis right 
plot(D, WT_rib_steady_1gL.rf  ./ WT_y_steady_1gL(:,num_y.r0), '-', 'Color', color_rt.default);   % free ribosome fraction
plot(D, WT_rib_steady_1gL.rat ./ WT_y_steady_1gL(:,num_y.r0), '-', 'Color', color_sc.kumar_hg); % active translating ribosome fraction
plot(D, WT_rib_steady_1gL.ras ./ WT_y_steady_1gL(:,num_y.r0), '-', 'Color', plt_clrs.gray);  % stalled translating ribosome fraction
ylabel('t_{c}, t_{u} fractions')
ylim([0 1])

legend({'total','free','translating','stalled'})
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%}

%xline(d_crit,'HandleVisibility','off')
xlim([D(1)-0.05 D(end)+0.05])
xlabel(x_label);

axis square; 
box on;
hold off; 

% translating ribosomes
%{
kumar_trans_rib = par.pro_den .* kumar_lg.dr ./ par.k_e; 
subplot(a, b, 7 + b)
hold on; 
plot(kumar_lg.dr, kumar_trans_rib, 'o', 'MarkerFaceColor', color_sc.kumar_hg)
plot(D, WT_rib_steady_1gL.rat, '-', 'Color', color_sc.default); % active translating ribosome fraction
ylabel('R_{at}')
%xline(d_crit,'HandleVisibility','off')
xlim([D(1)-0.05 D(end)+0.05])

axis square; 
box on;
hold off; 
%}

% free ribosomes    
%{
subplot(a, b, 7 + b)
hold on; 

yyaxis left 
plot(D, WT_rib_steady_1gL.rf, '-', 'Color', color_sc.default)
ylim([0 1.05*max(WT_rib_steady_1gL.rf)])
ylabel('R_{f}')

yyaxis right 
plot(D, WT_rib_steady_1gL.raf./WT_rib_steady_1gL.rf, '-', 'Color', color_rt.default)
plot(D, WT_rib_steady_1gL.ri./WT_rib_steady_1gL.rf, '-', 'Color', color_sc.kumar_hg)
ylim([0 1.05])
ylabel('R_{af}, R_{i} fractions')

legend({'total free', 'active free', 'inactive free'})
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%xline(d_crit,'HandleVisibility','off')
xlim([D(1)-0.05 D(end)+0.05])
xlabel(x_label)
axis square; 
box on;
hold off;
%}

% bound ribosomes          
subplot(a, b, 7 + b)
hold on; 

%yyaxis left 
transl_rib_1gL  = WT_y_steady_1gL(:,num_y.r0)  - WT_rib_steady_1gL.rf; 
transl_rib_10gL = WT_y_steady_10gL(:,num_y.r0) - WT_rib_steady_10gL.rf; 
transl_rib_08gL = WT_y_steady_08gL(:,num_y.r0) - WT_rib_steady_08gL.rf; 
plot(D, transl_rib_1gL/WT_y_steady_1gL(:,num_y.r0), '-', 'Color', color_sc.default)
plot(D, transl_rib_10gL/WT_y_steady_10gL(:,num_y.r0), '-', 'Color', color_sc.kumar_hg)
plot(D, transl_rib_08gL/WT_y_steady_08gL(:,num_y.r0), '-', 'Color', color_rt.default)
ylim([0 1])
yline(0.8)
%{
ymax1 = max(transl_rib_1gL/WT_y_steady_1gL(:,num_y.r0)); 
ymax2 = max(transl_rib_10gL/WT_y_steady_10gL(:,num_y.r0)); 
ymax3 = max(transl_rib_08gL/WT_y_steady_08gL(:,num_y.r0)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end 
%}
ylabel('R_{b}/R_{0}')

%{
yyaxis right 
plot(D, WT_rib_steady_1gL.rat./transl_rib_1gL, '-', 'Color', color_rt.default)
plot(D, WT_rib_steady_1gL.ras./transl_rib_1gL, '-', 'Color', color_sc.kumar_hg)
ylim([0 1.05])
ylabel('R_{at}, R_{as} fractions')

legend({'bound', 'active translating', 'stalled ribosomes'})
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%}
%xline(d_crit,'HandleVisibility','off')
xlim([D(1)-0.05 D(end)+0.05])
xlabel(x_label)
axis square; 
box on;
hold off;

%% tRNA 

% total tRNA
subplot(a, b, 7 + 2*b)
hold on;

%yyaxis left 
plot(D, WT_prot_syn_steady_1gL.t0, '-', 'Color', color_sc.default);
plot(D, WT_prot_syn_steady_10gL.t0, '-', 'Color', color_sc.kumar_hg);
plot(D, WT_prot_syn_steady_08gL.t0, '-', 'Color', color_rt.default);
ylabel('t_{0}')
ymax1 = max(WT_prot_syn_steady_1gL.t0); 
ymax2 = max(WT_prot_syn_steady_10gL.t0); 
ymax3 = max(WT_prot_syn_steady_08gL.t0);
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end 
%{
yyaxis right 
plot(D, WT_prot_syn_steady_1gL.tc./WT_prot_syn_steady_1gL.t0, '-', 'Color', color_rt.default)   % charged tRNA fraction
plot(D, WT_prot_syn_steady_1gL.tu./WT_prot_syn_steady_1gL.t0, '-', 'Color', color_sc.kumar_hg) % uncharged tRNA fraction
ylim([0 1.05])

legend({'total', 'charged', 'uncharged'})
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%}

%xline(d_crit,'HandleVisibility','off')
xlim([D(1)-0.05 D(end)+0.05])
xlabel(x_label);
axis square; 
box on;
hold off; 

%yyaxis right 

%{
% charged tRNA
subplot(a, b, 6 + b)
hold on;
plot(D, WT_prot_syn_steady_1gL.tc, '-', 'Color', color_sc.default);
ylabel('t_{c}')
xlabel(x_label_chem);
ylim([0 1.1.*max(WT_prot_syn_steady_1gL.tc)])
%xline(d_crit)
xline(t_eh_depl)
axis square; 
box on;
hold off; 

% uncharged tRNA
subplot(a, b, 6 + 2*b)
hold on;
plot(D, WT_prot_syn_steady_1gL.tu, '-', 'Color', color_sc.default);
ylabel('t_{u}')
xlabel(x_label_chem);
ylim([0 1.1.*max(WT_prot_syn_steady_1gL.tu)])
%xline(d_crit)
xline(t_eh_depl)
axis square; 
box on;
hold off; 
%}

%% growth rate

subplot(a, b, 3 + 4*b)
hold on;
plot(D, WT_g_rate_steady_1gL, '-', 'Color', color_sc.default);
plot(D, WT_g_rate_steady_10gL, '-', 'Color', color_sc.kumar_hg);
plot(D, WT_g_rate_steady_08gL, '-', 'Color', color_rt.default);
ylabel('\lambda')
xlabel(x_label);
ymax1 = max(WT_g_rate_steady_1gL); 
ymax2 = max(WT_g_rate_steady_10gL); 
ymax3 = max(WT_g_rate_steady_08gL); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.1.*max([ymax1, ymax2, ymax3])])
end 
%xline(d_crit)
%xlim([D(1)-0.05 D(end)+0.05])
axis square; 
box on;
hold off; 



% Biomass SC  % 
subplot(a, b, 4 + 3*b)
hold on;
if exist('carbo_sc','var')
papers = fieldnames(carbo_sc.sc); 
for i = 1: length(papers)
plot(D_sc.(papers{i}),  carbo_sc.sc.(papers{i}) ,...
shape_sc.(papers{i}),...
'Color', color_sc.(papers{i}), ...
'MarkerFaceColor', color_sc.(papers{i}))
end 
end 
plot(D,   100.*WT_y_steady_1gL(:,num_y.sc).*gl_uM_to_gPerL./par.rho_cell , 'Color', color_sc.kumar_lg);
plot(D,   100.*WT_y_steady_10gL(:,num_y.sc).*gl_uM_to_gPerL./par.rho_cell,   'Color', color_sc.kumar_hg);
plot(D,   100.*WT_y_steady_08gL(:,num_y.sc).*gl_uM_to_gPerL./par.rho_cell, 'Color', color_rt.default);
ylabel('SC (%)')
xlabel(x_label);
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
axis square; 
box on;
hold off; 


% Biomass SC uM
subplot(a, b, 2 + 4*b)
hold on;
if exist('carbo_sc','var')
papers = fieldnames(carbo_sc.sc); 
for i = 1: length(papers)
plot(D_sc.(papers{i}),  carbo_sc.sc.(papers{i}).*par.rho_cell.* gl_gPerL_to_uM./100 ,...
shape_sc.(papers{i}),...
'Color', color_sc.(papers{i}), ...
'MarkerFaceColor', color_sc.(papers{i}))
end 
end 
plot(D,   WT_y_steady_1gL(:,num_y.sc) , 'Color', color_sc.kumar_lg);
plot(D,   WT_y_steady_10gL(:,num_y.sc),   'Color', color_sc.kumar_hg);
plot(D,   WT_y_steady_08gL(:,num_y.sc), 'Color', color_rt.default);
ylabel('SC (uM)')
xlabel(x_label);
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
axis square; 
box on;
hold off; 


% Biomass (Lipid)
subplot(a, b, 4 + 4*b)
hold on;
plot(D,   WT_y_steady_1gL(:,num_y.lp_e)  *lp_uM_to_gPerL , 'Color', color_sc.kumar_lg);
plot(D,   WT_y_steady_10gL(:,num_y.lp_e)  *lp_uM_to_gPerL,   'Color', color_sc.kumar_hg);
plot(D,   WT_y_steady_08gL(:,num_y.lp_e) *lp_uM_to_gPerL, 'Color', color_rt.default);
ylabel('Lipid (g/L)')
xlabel(x_label);
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
axis square; 
box on;
hold off; 


sgtitle('chemostat culture')



%% ATP production 
%{
figure; 
a = 3;
b = 4; 

ferm_prot_frac_1gL  = (par.l(3).*WT_y_steady_1gL(:,num_y.gy)   + par.l(4).*WT_y_steady_1gL(:,num_y.fe))  ./ total_protein_con_1gL; 
ferm_prot_frac_10gL = (par.l(3).*WT_y_steady_10gL(:,num_y.gy)  + par.l(4).*WT_y_steady_10gL(:,num_y.fe)) ./ total_protein_con_10gL;
ferm_prot_frac_08gL = (par.l(3).*WT_y_steady_08gL(:,num_y.gy)  + par.l(4).*WT_y_steady_08gL(:,num_y.fe)) ./ total_protein_con_08gL;

resp_prot_frac_1gL  = (par.l(3).*WT_y_steady_1gL(:,num_y.gy)   + par.l(4).*WT_y_steady_1gL(:,num_y.mt))  ./ total_protein_con_1gL; 
resp_prot_frac_10gL = (par.l(3).*WT_y_steady_10gL(:,num_y.gy)  + par.l(4).*WT_y_steady_10gL(:,num_y.mt)) ./ total_protein_con_10gL;
resp_prot_frac_08gL = (par.l(3).*WT_y_steady_08gL(:,num_y.gy)  + par.l(4).*WT_y_steady_08gL(:,num_y.mt)) ./ total_protein_con_08gL;

ferm_flux_1gL  = WT_met_reac_steady_1gL.flux(:,num_flux.fe);
ferm_flux_10gL = WT_met_reac_steady_10gL.flux(:,num_flux.fe);
ferm_flux_08gL = WT_met_reac_steady_08gL.flux(:,num_flux.fe);

resp_flux_1gL  = WT_met_reac_steady_1gL.flux(:,num_flux.mt);
resp_flux_10gL = WT_met_reac_steady_10gL.flux(:,num_flux.mt);
resp_flux_08gL = WT_met_reac_steady_08gL.flux(:,num_flux.mt);

gluconeo_flux_1gL  = WT_met_reac_steady_1gL.flux(:,num_flux.gn);
gluconeo_flux_10gL = WT_met_reac_steady_10gL.flux(:,num_flux.gn);
gluconeo_flux_08gL = WT_met_reac_steady_08gL.flux(:,num_flux.gn);

as_flux_1gL  = WT_met_reac_steady_1gL.flux(:,num_flux.as);
as_flux_10gL = WT_met_reac_steady_10gL.flux(:,num_flux.as);
as_flux_08gL = WT_met_reac_steady_08gL.flux(:,num_flux.as);

sp_flux_1gL  = WT_met_reac_steady_1gL.flux(:,num_flux.sp);
sp_flux_10gL = WT_met_reac_steady_10gL.flux(:,num_flux.sp);
sp_flux_08gL = WT_met_reac_steady_08gL.flux(:,num_flux.sp);

sd_flux_1gL  = WT_met_reac_steady_1gL.flux(:,num_flux.sd);
sd_flux_10gL = WT_met_reac_steady_10gL.flux(:,num_flux.sd);
sd_flux_08gL = WT_met_reac_steady_08gL.flux(:,num_flux.sd);

prot_syn_flux_1gL  = WT_prot_syn_steady_1gL.prot_syn_rate(:); 
prot_syn_flux_10gL = WT_prot_syn_steady_10gL.prot_syn_rate(:);
prot_syn_flux_08gL = WT_prot_syn_steady_08gL.prot_syn_rate(:);

% atp_other_flux_1gL  = WT_other_met_reac_steady_1gL.flux.eo; 
% atp_other_flux_10gL = WT_other_met_reac_steady_10gL.flux.eo; 
% atp_other_flux_08gL = WT_other_met_reac_steady_08gL.flux.eo; 

% threshold ATP production flux (mmol/h/gdw converted to umol/h/cell)
atp_prod_thresh = 21.785095320623917 * 1000 * one_over_gdw_to_1_over_cell * (1/par.V_c);

% total protein concentration
%{
total_protein_con_1gL  = sum(WT_y_steady_1gL(:,num_y.r:num_y.ht).*par.l',2);
total_protein_con_10gL = sum(WT_y_steady_10gL(:,num_y.r:num_y.ht).*par.l',2);
total_protein_con_08gL = sum(WT_y_steady_08gL(:,num_y.r:num_y.ht).*par.l',2);
%}

% ATP production flux

total_ATP_prod_1gL  = par.q_gy.*WT_met_reac_steady_1gL.flux(:,num_flux.gy) ... 
+ par.q_mt.*WT_met_reac_steady_1gL.flux(:,num_flux.mt) ...
+ par.q_sd.*WT_met_reac_steady_1gL.flux(:,num_flux.sd) ...
+ par.q_gn.*WT_met_reac_steady_1gL.flux(:,num_flux.gn);

total_ATP_prod_10gL  = par.q_gy.*WT_met_reac_steady_10gL.flux(:,num_flux.gy) ... 
+ par.q_mt.*WT_met_reac_steady_10gL.flux(:,num_flux.mt) ...
+ par.q_sd.*WT_met_reac_steady_10gL.flux(:,num_flux.sd) ...
+ par.q_gn.*WT_met_reac_steady_10gL.flux(:,num_flux.gn);

total_ATP_prod_08gL  = par.q_gy.*WT_met_reac_steady_08gL.flux(:,num_flux.gy) ... 
+ par.q_mt.*WT_met_reac_steady_08gL.flux(:,num_flux.mt) ...
+ par.q_sd.*WT_met_reac_steady_08gL.flux(:,num_flux.sd) ...
+ par.q_gn.*WT_met_reac_steady_08gL.flux(:,num_flux.gn);

% ATP utilization     

total_ATP_use_1gL  = par.q_sp.*WT_met_reac_steady_1gL.flux(:,num_flux.sp) ... 
+ par.q_as.*WT_met_reac_steady_1gL.flux(:,num_flux.as) ...
+ par.q_at.*WT_met_reac_steady_1gL.flux(:,num_flux.at) ...
+ par.q_p .* prot_syn_flux_1gL ...
+ WT_other_met_reac_steady_1gL.flux.eo; 

total_ATP_use_10gL  = par.q_sp.*WT_met_reac_steady_10gL.flux(:,num_flux.sp) ... 
+ par.q_as.*WT_met_reac_steady_10gL.flux(:,num_flux.as) ...
+ par.q_at.*WT_met_reac_steady_10gL.flux(:,num_flux.at) ...
+ par.q_p .* prot_syn_flux_10gL ...
+ WT_other_met_reac_steady_10gL.flux.eo; 

total_ATP_use_08gL  = par.q_sp.*WT_met_reac_steady_08gL.flux(:,num_flux.sp) ... 
+ par.q_as.*WT_met_reac_steady_08gL.flux(:,num_flux.as) ...
+ par.q_at.*WT_met_reac_steady_08gL.flux(:,num_flux.at) ...
+ par.q_p .* prot_syn_flux_08gL ...
+ WT_other_met_reac_steady_08gL.flux.eo; 

% find ATP flux when respiration starts to decrease 
if find(max(resp_flux_1gL)) > 0
indx_max = find(resp_flux_1gL == max(resp_flux_1gL));
atp_val_1  = total_ATP_prod_1gL(indx_max);
else
atp_val_1 = 0; 
end 

if find(max(resp_flux_10gL)) > 0
indx_max = find(resp_flux_10gL == max(resp_flux_10gL));
atp_val_10  = total_ATP_prod_10gL(indx_max);
else
atp_val_10 = 0; 
end 

if find(max(resp_flux_08gL)) > 0
indx_max = find(resp_flux_08gL == max(resp_flux_08gL));
atp_val_08  = total_ATP_prod_08gL(indx_max);
else
atp_val_08 = 0; 
end 

% ATP flux vs. dilution rate
subplot(a, b, 1)
hold on;
plot(D, total_ATP_prod_1gL,  'Color', color_sc.kumar_lg,  'HandleVisibility', 'off');
plot(D, total_ATP_prod_10gL, 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off');
plot(D, total_ATP_prod_08gL, 'Color', color_sc.hackett,   'HandleVisibility', 'off');
ylabel('ATP flux');
xlabel('dilution rate (h^{-1})')
axis square; 
box on;
hold off;

xlabel('dilution rate (h^{-1})')
axis square; 
box on;
hold off;

% metabolites
subplot(a, b, 2+b)
hold on;
plot(0, 0, '-',  'Color', 'k')
plot(0, 0, '--', 'Color', 'k')
legend({'glucose'; 'ethanol'})

plot(total_ATP_prod_1gL,  WT_y_steady_1gL(:,num_y.gl),  'Color', color_sc.kumar_lg,  'HandleVisibility', 'off');
plot(total_ATP_prod_10gL, WT_y_steady_10gL(:,num_y.gl), 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off');
plot(total_ATP_prod_08gL, WT_y_steady_08gL(:,num_y.gl), 'Color', color_sc.hackett,   'HandleVisibility', 'off');

plot(total_ATP_prod_1gL,  WT_y_steady_1gL(:,num_y.eh),  '--', 'Color', color_sc.kumar_lg, 'HandleVisibility', 'off');
plot(total_ATP_prod_10gL, WT_y_steady_10gL(:,num_y.eh), '--', 'Color', color_sc.kumar_hg,  'HandleVisibility', 'off');
plot(total_ATP_prod_08gL, WT_y_steady_08gL(:,num_y.eh), '--', 'Color', color_sc.hackett,   'HandleVisibility', 'off');
xline(atp_prod_thresh, 'HandleVisibility', 'off');

%xline(atp_val_1,  'Color', color_sc.kumar_lg,  'HandleVisibility', 'off')
%xline(atp_val_10, 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off')
%xline(atp_val_08, 'Color', color_sc.hackett,   'HandleVisibility', 'off')

hold off;
xlabel('ATP production flux')
ylabel('metabolites');
axis square; 
box on;

subplot(a, b, 2)
hold on;
plot(0, 0, '-',  'Color', 'k')
plot(0, 0, '--', 'Color', 'k')
legend({'glucose'; 'ethanol'})

plot(D, WT_y_steady_1gL(:,num_y.gl),  'Color', color_sc.kumar_lg, 'HandleVisibility', 'off');
plot(D, WT_y_steady_10gL(:,num_y.gl), 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off');
plot(D, WT_y_steady_08gL(:,num_y.gl), 'Color', color_sc.hackett, 'HandleVisibility', 'off');

plot(D, WT_y_steady_1gL(:,num_y.eh), '--', 'Color', color_sc.kumar_lg, 'HandleVisibility', 'off');
plot(D, WT_y_steady_10gL(:,num_y.eh), '--', 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off');
plot(D, WT_y_steady_08gL(:,num_y.eh), '--', 'Color', color_sc.hackett, 'HandleVisibility', 'off');
ylabel('metabolites');
xlabel('dilution rate (h^{-1})')
axis square; 
box on;
hold off;

% fluxes 
subplot(a, b, 3+b) 
hold on; 
plot(0, 0, '-',  'Color', 'k')
plot(0, 0, '--', 'Color', 'k')
legend({'J_{gy}'; 'J_{fe} - J_{gn}'})

plot(total_ATP_prod_1gL,  WT_met_reac_steady_1gL.flux(:,num_flux.gy),  'Color', color_sc.kumar_lg, 'HandleVisibility', 'off');
plot(total_ATP_prod_10gL, WT_met_reac_steady_10gL.flux(:,num_flux.gy), 'Color', color_sc.kumar_hg,  'HandleVisibility', 'off');
plot(total_ATP_prod_08gL, WT_met_reac_steady_08gL.flux(:,num_flux.gy), 'Color', color_sc.hackett,   'HandleVisibility', 'off');

plot(total_ATP_prod_1gL,  WT_met_reac_steady_1gL.flux(:,num_flux.fe)  - WT_met_reac_steady_1gL.flux(:,num_flux.gn),  '--', 'Color', color_sc.kumar_lg, 'HandleVisibility', 'off');
plot(total_ATP_prod_10gL, WT_met_reac_steady_10gL.flux(:,num_flux.fe) - WT_met_reac_steady_10gL.flux(:,num_flux.gn), '--', 'Color', color_sc.kumar_hg,  'HandleVisibility', 'off');
plot(total_ATP_prod_08gL, WT_met_reac_steady_08gL.flux(:,num_flux.fe) - WT_met_reac_steady_08gL.flux(:,num_flux.gn), '--', 'Color', color_sc.hackett,   'HandleVisibility', 'off');
xline(atp_prod_thresh, 'HandleVisibility', 'off');
%xline(atp_val_1,  'Color', color_sc.kumar_lg,  'HandleVisibility', 'off')
%xline(atp_val_10, 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off')
%xline(atp_val_08, 'Color', color_sc.hackett,   'HandleVisibility', 'off')
hold off; 
xlabel('ATP production flux')
ylabel('fluxes');
axis square; 
box on; 


subplot(a, b, 3+2*b)
hold on; 
plot(0, 0, '-',  'Color', 'k')
plot(0, 0, '--', 'Color', 'k')
legend({'J_{gy} + J_{fe}'; 'J_{gy} + J_{mt}'})

plot(total_ATP_prod_1gL,  par.q_gy * ferm_flux_1gL,  'Color', color_sc.kumar_lg, 'HandleVisibility', 'off');
plot(total_ATP_prod_10gL, par.q_gy * ferm_flux_10gL, 'Color', color_sc.kumar_hg,  'HandleVisibility', 'off');
plot(total_ATP_prod_08gL, par.q_gy * ferm_flux_08gL, 'Color', color_sc.hackett,   'HandleVisibility', 'off');

plot(total_ATP_prod_1gL,  (par.q_gy + par.q_mt) * resp_flux_1gL,  '--', 'Color', color_sc.kumar_lg, 'HandleVisibility', 'off');
plot(total_ATP_prod_10gL, (par.q_gy + par.q_mt) * resp_flux_10gL, '--', 'Color', color_sc.kumar_hg,  'HandleVisibility', 'off');
plot(total_ATP_prod_08gL, (par.q_gy + par.q_mt) * resp_flux_08gL, '--', 'Color', color_sc.hackett,   'HandleVisibility', 'off');
xline(atp_prod_thresh, 'HandleVisibility', 'off');
%xline(atp_val_1,  'Color', color_sc.kumar_lg,  'HandleVisibility', 'off')
%xline(atp_val_10, 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off')
%xline(atp_val_08, 'Color', color_sc.hackett,   'HandleVisibility', 'off')
hold off; 
xlabel('ATP production flux')
ylabel('energy production fluxes');
axis square; 
box on; 

% fermentation and respiration proteins vs. total ATP production flux 
subplot(a, b, 1+2*b)
hold on;
plot(0, 0, '-',  'Color', 'k')
plot(0, 0, '--', 'Color', 'k')
legend({'E_{gy} + E_{fe} fraction'; 'E_{gy} + E_{mt} fraction'})

plot(total_ATP_prod_1gL,  ferm_prot_frac_1gL,  'Color', color_sc.kumar_lg,  'HandleVisibility', 'off');
plot(total_ATP_prod_10gL, ferm_prot_frac_10gL, 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off');
plot(total_ATP_prod_08gL, ferm_prot_frac_08gL, 'Color', color_sc.hackett,   'HandleVisibility', 'off');

plot(total_ATP_prod_1gL,  resp_prot_frac_1gL,  '--', 'Color', color_sc.kumar_lg,  'HandleVisibility', 'off');
plot(total_ATP_prod_10gL, resp_prot_frac_10gL, '--', 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off');
plot(total_ATP_prod_08gL, resp_prot_frac_08gL, '--', 'Color', color_sc.hackett,   'HandleVisibility', 'off');

xline(atp_prod_thresh, 'HandleVisibility', 'off');
%xline(atp_val_1,  'Color', color_sc.kumar_lg,  'HandleVisibility', 'off')
%xline(atp_val_10, 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off')
%xline(atp_val_08, 'Color', color_sc.hackett,   'HandleVisibility', 'off')

hold off;
ylabel('energy metabolism fractions')
xlabel('ATP flux');

ymax1 = max(ferm_prot_frac_1gL);
ymax2 = max(ferm_prot_frac_10gL);
ymax3 = max(ferm_prot_frac_08gL);
ymax4 = max(resp_prot_frac_1gL);
ymax5 = max(resp_prot_frac_10gL);
ymax6 = max(resp_prot_frac_08gL);
ylim([0 1.1*max([ymax1, ymax2, ymax3, ymax4, ymax5, ymax6])])

axis square; 
box on;

% fermentation and respiration proteins vs. total ATP production flux 
subplot(a, b, 2+2*b)
hold on;
plot(total_ATP_prod_1gL,  ferm_prot_frac_1gL  + resp_prot_frac_1gL  - (par.l(3).*WT_y_steady_1gL(:,num_y.gy))./ total_protein_con_1gL, 'Color', color_sc.kumar_lg,  'HandleVisibility', 'off')
plot(total_ATP_prod_10gL, ferm_prot_frac_10gL + resp_prot_frac_10gL - (par.l(3).*WT_y_steady_1gL(:,num_y.gy))./ total_protein_con_1gL, 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off')
plot(total_ATP_prod_08gL, ferm_prot_frac_08gL + resp_prot_frac_08gL - (par.l(3).*WT_y_steady_1gL(:,num_y.gy))./ total_protein_con_1gL, 'Color', color_sc.hackett,   'HandleVisibility', 'off')

xline(atp_prod_thresh, 'HandleVisibility', 'off');

hold off;
ylabel('total energy metabolism fraction')
xlabel('ATP flux');
ylim([0 1])

axis square; 
box on;

% efficiencies 
subplot(a, b, 1+b) 
hold on; 
plot(0, 0, '-',  'Color', 'k')
plot(0, 0, '--', 'Color', 'k')
legend({'fermentation'; 'respiration'})

plot(total_ATP_prod_1gL,  ferm_flux_1gL  / ferm_prot_frac_1gL,  'Color', color_sc.kumar_lg,  'HandleVisibility', 'off');
plot(total_ATP_prod_10gL, ferm_flux_10gL / ferm_prot_frac_10gL, 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off');
plot(total_ATP_prod_08gL, ferm_flux_08gL / ferm_prot_frac_08gL, 'Color', color_sc.hackett,   'HandleVisibility', 'off');

plot(total_ATP_prod_1gL,  resp_flux_1gL  / resp_prot_frac_1gL,  '--', 'Color', color_sc.kumar_lg,  'HandleVisibility', 'off');
plot(total_ATP_prod_10gL, resp_flux_10gL / resp_prot_frac_10gL, '--', 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off');
plot(total_ATP_prod_08gL, resp_flux_08gL / resp_prot_frac_08gL, '--', 'Color', color_sc.hackett,   'HandleVisibility', 'off');
xline(atp_prod_thresh, 'HandleVisibility', 'off');
%xline(atp_val_1,  'Color', color_sc.kumar_lg,  'HandleVisibility', 'off')
%xline(atp_val_10, 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off')
%xline(atp_val_08, 'Color', color_sc.hackett,   'HandleVisibility', 'off')
hold off; 
ylabel('energy metabolism efficiencies')
xlabel('ATP flux');

ymax1 = max(ferm_flux_1gL  / ferm_prot_frac_1gL);
ymax2 = max(ferm_flux_10gL / ferm_prot_frac_10gL);
ymax3 = max(ferm_flux_08gL / ferm_prot_frac_08gL);
ymax4 = max(resp_flux_1gL  / resp_prot_frac_1gL);
ymax5 = max(resp_flux_10gL / resp_prot_frac_10gL);
ymax6 = max(resp_flux_08gL / resp_prot_frac_08gL);
ylim([0 1.1*max([ymax1, ymax2, ymax3, ymax4, ymax5, ymax6])])

axis square; 
box on;



subplot(a, b, 4)
hold on; 
plot(0, 0, '-',  'Color', 'k')
plot(0, 0, '--', 'Color', 'k')
plot(0, 0, '-.', 'Color', 'k')
legend({'total ATP usage/total ATP production'; '(qgy + +qmy)*J_{mt}/total ATP prod'; 'qgy*J_{fe}'})

plot(D, total_ATP_use_1gL  ./ total_ATP_prod_1gL,  'Color', color_sc.kumar_lg,  'HandleVisibility', 'off');
plot(D, total_ATP_use_10gL ./ total_ATP_prod_10gL, 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off');
plot(D, total_ATP_use_08gL ./ total_ATP_prod_08gL, 'Color', color_sc.hackett,   'HandleVisibility', 'off');

plot(D, ((par.q_gy + par.q_mt) * resp_flux_1gL)  ./ total_ATP_prod_1gL,  '--', 'Color', color_sc.kumar_lg,  'HandleVisibility', 'off');
plot(D, ((par.q_gy + par.q_mt) * resp_flux_10gL) ./ total_ATP_prod_10gL, '--', 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off');
plot(D, ((par.q_gy + par.q_mt) * resp_flux_08gL) ./ total_ATP_prod_08gL, '--', 'Color', color_sc.hackett,   'HandleVisibility', 'off');

plot(D, ((par.q_gy) * ferm_flux_1gL)  ./ total_ATP_prod_1gL,  '-.', 'Color', color_sc.kumar_lg,  'HandleVisibility', 'off');
plot(D, ((par.q_gy) * ferm_flux_10gL) ./ total_ATP_prod_10gL, '-.', 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off');
plot(D, ((par.q_gy) * ferm_flux_08gL) ./ total_ATP_prod_08gL, '-.', 'Color', color_sc.hackett,   'HandleVisibility', 'off');

%xline(atp_prod_thresh, 'HandleVisibility', 'off');
%xline(atp_val_1,  'Color', color_sc.kumar_lg,  'HandleVisibility', 'off')
%xline(atp_val_10, 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off')
%xline(atp_val_08, 'Color', color_sc.hackett,   'HandleVisibility', 'off')
hold off; 
xlabel(x_label)
ylabel('energy fractions');
axis square; 
box on; 

subplot(a, b, 4+b)
hold on; 
plot(0, 0, '-',  'Color', 'k')
plot(0, 0, '-.', 'Color', 'k')
legend({'ATP from respiration'; 'total ATP usage'})

plot(D, ((par.q_gy + par.q_mt) * resp_flux_1gL),  'Color', color_sc.kumar_lg,  'HandleVisibility', 'off');
plot(D, ((par.q_gy + par.q_mt) * resp_flux_10gL), 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off');
plot(D, ((par.q_gy + par.q_mt) * resp_flux_08gL), 'Color', color_sc.hackett,   'HandleVisibility', 'off');
plot(D, total_ATP_use_1gL,  '-.', 'Color', color_sc.kumar_lg,  'HandleVisibility', 'off')
plot(D, total_ATP_use_10gL, '-.', 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off')
plot(D, total_ATP_use_08gL, '-.', 'Color', color_sc.kumar_lg,  'HandleVisibility', 'off')
%xline(atp_prod_thresh, 'HandleVisibility', 'off');
%xline(atp_val_1,  'Color', color_sc.kumar_lg,  'HandleVisibility', 'off')
%xline(atp_val_10, 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off')
%xline(atp_val_08, 'Color', color_sc.hackett,   'HandleVisibility', 'off')
hold off; 
xlabel(x_label)
ylabel('energy from respiration');
axis square; 
box on; 

subplot(a, b, 4+2*b)
hold on; 
plot(0, 0, '-',  'Color', 'k')
plot(0, 0, '-.', 'Color', 'k')
legend({'ATP from fermentation'; 'total ATP usage'})

plot(D, ((par.q_gy) * ferm_flux_1gL),  'Color', color_sc.kumar_lg,  'HandleVisibility', 'off');
plot(D, ((par.q_gy) * ferm_flux_10gL), 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off');
plot(D, ((par.q_gy) * ferm_flux_08gL), 'Color', color_sc.hackett,   'HandleVisibility', 'off');
plot(D, total_ATP_use_1gL,  '-.', 'Color', color_sc.kumar_lg,  'HandleVisibility', 'off')
plot(D, total_ATP_use_10gL, '-.', 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off')
plot(D, total_ATP_use_08gL, '-.', 'Color', color_sc.kumar_lg,  'HandleVisibility', 'off')
%xline(atp_prod_thresh, 'HandleVisibility', 'off');
%xline(atp_val_1,  'Color', color_sc.kumar_lg,  'HandleVisibility', 'off')
%xline(atp_val_10, 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off')
%xline(atp_val_08, 'Color', color_sc.hackett,   'HandleVisibility', 'off')
hold off; 
xlabel(x_label)
ylabel('energy from fermentation');
axis square; 
box on; 

subplot(a, b, 3)
hold on; 
plot(0, 0, '-',  'Color', 'k')
plot(0, 0, '--', 'Color', 'k')
legend({'total ATP usage - total ATP produced'})

plot(D, -total_ATP_use_1gL  + total_ATP_prod_1gL,  'Color', color_sc.kumar_lg,  'HandleVisibility', 'off');
plot(D, -total_ATP_use_10gL + total_ATP_prod_10gL, 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off');
plot(D, -total_ATP_use_08gL + total_ATP_prod_08gL, 'Color', color_sc.hackett,   'HandleVisibility', 'off');
%xline(atp_prod_thresh, 'HandleVisibility', 'off');
%xline(atp_val_1,  'Color', color_sc.kumar_lg,  'HandleVisibility', 'off')
%xline(atp_val_10, 'Color', color_sc.kumar_hg, 'HandleVisibility', 'off')
%xline(atp_val_08, 'Color', color_sc.hackett,   'HandleVisibility', 'off')
hold off; 
xlabel(x_label)
ylabel('energy balance');
axis square; 
box on; 

%}


% bound ribosomes          
figure;

subplot(1,2,1)
hold on; 
%yyaxis left 
transl_rib_1gL  = WT_y_steady_1gL(:,num_y.r0)  - WT_rib_steady_1gL.rf; 
transl_rib_10gL = WT_y_steady_10gL(:,num_y.r0) - WT_rib_steady_10gL.rf; 
transl_rib_08gL = WT_y_steady_08gL(:,num_y.r0) - WT_rib_steady_08gL.rf; 
plot(D, 1- transl_rib_1gL./WT_y_steady_1gL(:,num_y.r0), '-', 'Color', color_sc.default)
plot(D, 1- transl_rib_10gL./WT_y_steady_10gL(:,num_y.r0), '-', 'Color', color_sc.kumar_hg)
plot(D, 1- transl_rib_08gL./WT_y_steady_08gL(:,num_y.r0), '-', 'Color', color_rt.default)
hold off 
ylim([0 1])
ylabel('(1-(R_{f}/R_{0}))/R_{0}')
xlim([D(1)-0.05 D(end)+0.05])
xlabel(x_label);
axis square; 
box on;
hold off; 

subplot(1,2,2)
hold on; 

%yyaxis left 
transl_rib_1gL  = WT_y_steady_1gL(:,num_y.r0)  - WT_rib_steady_1gL.rf; 
transl_rib_10gL = WT_y_steady_10gL(:,num_y.r0) - WT_rib_steady_10gL.rf; 
transl_rib_08gL = WT_y_steady_08gL(:,num_y.r0) - WT_rib_steady_08gL.rf; 
plot(D, WT_rib_steady_1gL.rf  ./ WT_y_steady_1gL(:,num_y.r0),  '-', 'Color', color_sc.default)
plot(D, WT_rib_steady_10gL.rf ./ WT_y_steady_10gL(:,num_y.r0), '-', 'Color', color_sc.kumar_hg)
plot(D, WT_rib_steady_08gL.rf ./ WT_y_steady_08gL(:,num_y.r0), '-', 'Color', color_rt.default)
hold off 
ylim([0 1])
ylabel('R_{f}/R_{0}')
xlim([D(1)-0.05 D(end)+0.05])
xlabel(x_label);
axis square; 
box on;
hold off; 


%% protein faction --------------------------------------------------------------------------
figure; 
a = 3; b = 5; 
prot_ = fieldnames(num_prot); 
for i = 1: max(struct2array(num_prot))
subplot(a, b, i)
hold on;
%plot(hackett.dr, hackett.r, 'o', 'MarkerFaceColor', color_sc.hackett)
plot(D, WT_y_steady_1gL(:,num_y.(prot_{i}))  , '-', 'Color', color_sc.default);
plot(D, WT_y_steady_10gL(:,num_y.(prot_{i})) , '-', 'Color', color_sc.kumar_hg);
plot(D, WT_y_steady_08gL(:,num_y.(prot_{i})) , '-', 'Color', color_rt.default);
ylabel(label.prot{num_prot.(prot_{i})})
xlabel(x_label);
if max(ylim) > 0 
ylim([0 1.1.*max(ylim)])
end 
xlim([D(1)-0.05 D(end)+0.05])
%xline(d_crit)
axis square; 
box on;
hold off;
end 



% fraction (clim) vs fraction (nlim) at lowest dilution rate 
subplot(a, b, i+1)
hold on;
%plot(hackett.dr, hackett.lo, 'o', 'MarkerFaceColor', color_sc.hackett)
plot(100*WT_y_steady_1gL(1,num_y.r), 100*WT_y_steady_08gL(1,num_y.r) , 'o');
plot(100*WT_y_steady_1gL(1,num_y.z), 100*WT_y_steady_08gL(1,num_y.z) , 'o');
plot(100*WT_y_steady_1gL(1,num_y.gy), 100*WT_y_steady_08gL(1,num_y.gy) , 'o');
plot(100*WT_y_steady_1gL(1,num_y.fe), 100*WT_y_steady_08gL(1,num_y.fe) , 'o');

plot(100*WT_y_steady_1gL(1,num_y.gn), 100*WT_y_steady_08gL(1,num_y.gn) , 'o');
plot(100*WT_y_steady_1gL(1,num_y.mt), 100*WT_y_steady_08gL(1,num_y.mt) , 'o');
plot(100*WT_y_steady_1gL(1,num_y.as), 100*WT_y_steady_08gL(1,num_y.as) , 'o');
plot(100*WT_y_steady_1gL(1,num_y.at), 100*WT_y_steady_08gL(1,num_y.at) , 'o');
plot(100*WT_y_steady_1gL(1,num_y.lp), 100*WT_y_steady_08gL(1,num_y.lp) , 'o');
plot(100*WT_y_steady_1gL(1,num_y.lo), 100*WT_y_steady_08gL(1,num_y.lo) , 'o');

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

% fraction (clim) vs fraction (nlim) at 0.2 dilution rate like Yu2021
D_exp = 0.2 ; 
[norm_D, norm_indx] = min(abs(D - D_exp)); 

subplot(a, b, i+2)
hold on;
%plot(hackett.dr, hackett.lo, 'o', 'MarkerFaceColor', color_sc.hackett)
plot(100*WT_y_steady_1gL(norm_indx,num_y.r), 100*WT_y_steady_08gL(norm_indx,num_y.r) , 'o');
plot(100*WT_y_steady_1gL(norm_indx,num_y.z), 100*WT_y_steady_08gL(norm_indx,num_y.z) , 'o');
plot(100*WT_y_steady_1gL(norm_indx,num_y.gy), 100*WT_y_steady_08gL(norm_indx,num_y.gy) , 'o');
plot(100*WT_y_steady_1gL(norm_indx,num_y.fe), 100*WT_y_steady_08gL(norm_indx,num_y.fe) , 'o');

plot(100*WT_y_steady_1gL(norm_indx,num_y.gn), 100*WT_y_steady_08gL(norm_indx,num_y.gn) , 'o');
plot(100*WT_y_steady_1gL(norm_indx,num_y.mt), 100*WT_y_steady_08gL(norm_indx,num_y.mt) , 'o');
plot(100*WT_y_steady_1gL(norm_indx,num_y.as), 100*WT_y_steady_08gL(norm_indx,num_y.as) , 'o');
plot(100*WT_y_steady_1gL(norm_indx,num_y.at), 100*WT_y_steady_08gL(norm_indx,num_y.at) , 'o');
plot(100*WT_y_steady_1gL(norm_indx,num_y.lp), 100*WT_y_steady_08gL(norm_indx,num_y.lp) , 'o');
plot(100*WT_y_steady_1gL(norm_indx,num_y.lo), 100*WT_y_steady_08gL(norm_indx,num_y.lo) , 'o');
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
legend(label.prot(num_prot.r:num_prot.lo))
sgtitle('protein fraction')



set(gcf, 'Position', [1441 266 932 532])



%% protein concentration --------------------------------------------------------------------------
figure; 
a = 4; b = 4; 
prot_ = fieldnames(num_prot); 
for i = 1: max(struct2array(num_prot))
subplot(a, b, i)
hold on;
%plot(hackett.dr, hackett.r, 'o', 'MarkerFaceColor', color_sc.hackett)
plot(D, WT_y_steady_1gL(:,num_y.(prot_{i})).*total_protein_con_1gL  , '-', 'Color', color_sc.default);
plot(D, WT_y_steady_10gL(:,num_y.(prot_{i})).*total_protein_con_10gL , '-', 'Color', color_sc.kumar_hg);
plot(D, WT_y_steady_08gL(:,num_y.(prot_{i})).*total_protein_con_08gL  , '-', 'Color', color_rt.default);
ylabel(label.prot{num_prot.(prot_{i})})
xlabel(x_label);
if max(ylim) > 0 
ylim([0 1.1.*max(ylim)])
end 
xlim([D(1)-0.05 D(end)+0.05])
%xline(d_crit)
axis square; 
box on;
hold off;
end 


% total concentration 
subplot(a, b, i+1)
hold on;
%plot(hackett.dr, hackett.lo, 'o', 'MarkerFaceColor', color_sc.hackett)
plot(D, total_protein_con_1gL  , '-', 'Color', color_sc.default);
plot(D, total_protein_con_10gL , '-', 'Color', color_sc.kumar_hg);
plot(D, total_protein_con_08gL, '-', 'Color', color_rt.default);
ylabel('sum prot con')
xlabel(x_label);
if max(ylim) > 0 
ylim([0 1.1.*max(ylim)])
end
xlim([D(1)-0.05 D(end)+0.05])
axis square; 
box on;
hold off;

sgtitle('protein concentration')

set(gcf, 'Position', [1441 161 698 637])


%% protein alpha --------------------------------------------------------------------------
figure; 
a = 4; b = 4; 
prot_ = fieldnames(num_prot); 
for i = 1: max(struct2array(num_prot))
subplot(a, b, i)
hold on;
%plot(hackett.dr, hackett.r, 'o', 'MarkerFaceColor', color_sc.hackett)
plot(D, WT_prot_syn_steady_1gL.alpha(:,num_prot.(prot_{i}))  , '-', 'Color', color_sc.default);
plot(D, WT_prot_syn_steady_10gL.alpha(:,num_prot.(prot_{i})) , '-', 'Color', color_sc.kumar_hg);
plot(D, WT_prot_syn_steady_08gL.alpha(:,num_prot.(prot_{i}))  , '-', 'Color', color_rt.default);
ylabel(label.prot{num_prot.(prot_{i})})
xlabel(x_label);
if max(ylim) > 0 
ylim([0 1.1.*max(ylim)])
end 
xlim([D(1)-0.05 D(end)+0.05])
%xline(d_crit)
axis square; 
box on;
hold off;
end 

sgtitle('protein \alpha')

set(gcf, 'Position', [1441 161 698 637])



end 
