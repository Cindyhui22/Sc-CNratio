function plot_all_chemostat_intermediates(D, ...
WT_y_steady_1gL, WT_y_steady_10gL, ...
WT_met_reac_steady_1gL, WT_met_reac_steady_10gL, ...
WT_other_met_reac_steady_1gL, WT_other_met_reac_steady_10gL, ...
WT_prot_syn_steady_1gL, WT_prot_syn_steady_10gL, ...
WT_rib_steady_1gL, WT_rib_steady_10gL, ...
label, plt_clrs, par, num_y, num_flux, x_label)
data_chemostat;  
data_hackett_prot; 
plot_num; 
paper_figure_label; 
% tRNA contributions
%tc_contri = bart_WT_batch_min_prot_syn.tc./par.K_bart_WT_batch_min_rib_tc; 
%tu_contri = bart_WT_batch_min_prot_syn.tu./par.K_bart_WT_batch_min_rib_tu;
%rat_tc_contri = tc_contri./(1 + tc_contri + tu_contri);
%ras_tu_contri = tu_contri./(1 + tc_contri + tu_contri);

%% real metabolites 
figure; 
a = 3+3;
b = 6;
fig_index = 1:b;

%--------------------------------------------------------------------------
%                               glucose 
%--------------------------------------------------------------------------

% glucose
subplot(a, b, fig_index(1));
hold on; 
plot(D, WT_y_steady_1gL(:,num_y.gl), 'Color', plt_clrs.green);
plot(D, WT_y_steady_10gL(:,num_y.gl), 'Color', plt_clrs.blue);
ylabel(label.met{num_met.gl})    
xlabel(x_label)
ymax1 = max(WT_y_steady_1gL(:,num_y.gl)); 
ymax2 = max(WT_y_steady_10gL(:,num_y.gl)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on;
hold off; 

% J_gy
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.flux(:,num_flux.gy), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.gy), 'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.gy})    
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.flux(:,num_flux.gy)); 
ymax2 = max(WT_met_reac_steady_10gL.flux(:,num_flux.gy)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_gy
subplot(a, b, fig_index(3)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.prot(:,num_flux.gy), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.prot(:,num_flux.gy), 'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.gy})
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.prot(:,num_flux.gy)); 
ymax2 = max(WT_met_reac_steady_10gL.prot(:,num_flux.gy)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Gl/(km+Gl) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.substrate(:,num_flux.gy), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.substrate(:,num_flux.gy), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.gy}); 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.substrate(:,num_flux.gy)); 
ymax2 = max(WT_met_reac_steady_10gL.substrate(:,num_flux.gy)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.sig(:,num_flux.gy), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.sig(:,num_flux.gy), 'Color', plt_clrs.blue);
ylabel(label.met_reac.sig{num_flux.gy}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.sig(:,num_flux.gy)); 
ymax2 = max(WT_met_reac_steady_10gL.sig(:,num_flux.gy)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.atp(:,num_flux.gy), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.atp(:,num_flux.gy), 'Color', plt_clrs.blue);
ylabel(label.met_reac.atp{num_flux.gy}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.atp(:,num_flux.gy)); 
ymax2 = max(WT_met_reac_steady_10gL.atp(:,num_flux.gy)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

fig_index = fig_index + b; 
%--------------------------------------------------------------------------
%                              ethanol
%--------------------------------------------------------------------------
% ethanol
subplot(a, b, fig_index(1));
hold on; 
plot(D, WT_y_steady_1gL(:,num_y.eh), 'Color', plt_clrs.green);
plot(D, WT_y_steady_10gL(:,num_y.eh), 'Color', plt_clrs.blue);
ylabel(label.met{num_met.eh})    
xlabel(x_label)
ymax1 = max(WT_y_steady_1gL(:,num_y.eh)); 
ymax2 = max(WT_y_steady_10gL(:,num_y.eh)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on;
hold off; 

% J_fe 
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.flux(:,num_flux.fe), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.fe})    
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.flux(:,num_flux.fe)); 
ymax2 = max(WT_met_reac_steady_10gL.flux(:,num_flux.fe)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_fe
subplot(a, b, fig_index(3)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.prot(:,num_flux.fe), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.prot(:,num_flux.fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.fe})
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.prot(:,num_flux.fe)); 
ymax2 = max(WT_met_reac_steady_10gL.prot(:,num_flux.fe)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.substrate(:,num_flux.fe), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.substrate(:,num_flux.fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.fe}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.substrate(:,num_flux.fe)); 
ymax2 = max(WT_met_reac_steady_10gL.substrate(:,num_flux.fe)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

%{
% signaling contribution 
subplot(a, b, fig_index(5)); 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.sig(:,num_flux.fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.sig{num_flux.fe}) 
xlabel(x_label)
ylim([0 1.05.*max(bart_WT_batch_min_met_reac.sig(:,num_flux.fe))])
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 


% ATP contribution 
subplot(a, b, fig_index(6)); 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.atp(:,num_flux.fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.atp{num_flux.fe}) 
xlabel(x_label)
ylim([0 1.05.*max(bart_WT_batch_min_met_reac.atp(:,num_flux.fe))])
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
%}
fig_index = fig_index + b; 

%--------------------------------------------------------------------------
% ethanol
%{
subplot(a, b, fig_index(1));
plot(bart_WT_batch_min_t, y(:,4), 'Color', plt_clrs.blue);
ylabel(label.met(4))    
xlabel(x_label)
axis square; 
box on;
%}

% J_gn 
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.flux(:,num_flux.gn), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.gn), 'Color', plt_clrs.blue);
ylabel('J_{gn}')    
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.flux(:,num_flux.gn)); 
ymax2 = max(WT_met_reac_steady_10gL.flux(:,num_flux.gn)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_gn
subplot(a, b, fig_index(3)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.prot(:,num_flux.gn), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.prot(:,num_flux.gn), 'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.gn})
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.prot(:,num_flux.gn)); 
ymax2 = max(WT_met_reac_steady_10gL.prot(:,num_flux.gn)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Eh/(km+Eh) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.substrate(:,num_flux.gn), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.substrate(:,num_flux.gn), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.gn}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.substrate(:,num_flux.gn)); 
ymax2 = max(WT_met_reac_steady_10gL.substrate(:,num_flux.gn)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.sig(:,num_flux.gn), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.sig(:,num_flux.gn), 'Color', plt_clrs.blue);
ylabel(label.met_reac.sig{num_flux.gn}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.sig(:,num_flux.gn)); 
ymax2 = max(WT_met_reac_steady_10gL.sig(:,num_flux.gn)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

%{
% ATP contribution 
subplot(a, b, fig_index(6)); 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.atp(:,num_flux.gn), 'Color', plt_clrs.blue);
ylabel(label.met_reac.atp{num_flux.gn}) 
xlabel(x_label)
ylim([0 1.05.*max(bart_WT_batch_min_met_reac.atp(:,num_flux.gn))])
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
%}

%--------------------------------------------------------------------------
%                              Lipid
%--------------------------------------------------------------------------
fig_index = fig_index + b; 
% lipid
subplot(a, b, fig_index(1));
hold on; 
plot(D,   WT_y_steady_1gL(:,num_y.lp_e),   'Color', plt_clrs.green);
plot(D,   WT_y_steady_10gL(:,num_y.lp_e),   'Color', plt_clrs.blue);
ylabel(label.met{num_met.lp_e})    
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on;
hold off; 

% J_lp,fe
subplot(a, b, fig_index(2)); 
hold on; 
plot(D,   WT_met_reac_steady_1gL.flux(:,num_flux.lp_fe),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.flux(:,num_flux.lp_fe),   'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.lp_fe})    
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_lp
subplot(a, b, fig_index(3)); 
hold on; 
plot(D,   WT_met_reac_steady_1gL.prot(:,num_flux.lp_fe),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.prot(:,num_flux.lp_fe),   'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.lp_fe})
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D,   WT_met_reac_steady_1gL.substrate(:,num_flux.lp_fe),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.substrate(:,num_flux.lp_fe),   'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.lp_fe}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 


% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D,   WT_met_reac_steady_1gL.sig(:,num_flux.lp_fe),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.sig(:,num_flux.lp_fe),   'Color', plt_clrs.blue);
ylabel(label.met_reac.sig{num_flux.lp_fe}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

%
% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on 
plot(D,   WT_met_reac_steady_1gL.atp(:,num_flux.lp_fe),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.atp(:,num_flux.lp_fe),   'Color', plt_clrs.blue);
hold off
ylabel(label.met_reac.atp{num_flux.lp_fe}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 


% J_lp,cit
fig_index = fig_index + b; 
subplot(a, b, fig_index(2)); 
hold on; 
plot(D,   WT_met_reac_steady_1gL.flux(:,num_flux.lp_cit),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.flux(:,num_flux.lp_cit),   'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.lp_cit})    
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_lp
subplot(a, b, fig_index(3)); 
hold on; 
plot(D,   WT_met_reac_steady_1gL.prot(:,num_flux.lp_cit),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.prot(:,num_flux.lp_cit),   'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.lp_cit})
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D,   WT_met_reac_steady_1gL.substrate(:,num_flux.lp_cit),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.substrate(:,num_flux.lp_cit),   'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.lp_cit}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 


% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D,   WT_met_reac_steady_1gL.sig(:,num_flux.lp_cit),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.sig(:,num_flux.lp_cit),   'Color', plt_clrs.blue);
ylabel(label.met_reac.sig{num_flux.lp_cit}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

%
% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on 
plot(D,   WT_met_reac_steady_1gL.atp(:,num_flux.lp_cit),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.atp(:,num_flux.lp_cit),   'Color', plt_clrs.blue);
hold off
ylabel(label.met_reac.atp{num_flux.lp_cit}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 


fig_index = fig_index + b; 

%--------------------------------------------------------------------------

% J_lo 
subplot(a, b, fig_index(2)); 
hold on; 
plot(D,   WT_met_reac_steady_1gL.flux(:,num_flux.lo),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.flux(:,num_flux.lo),   'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.lo})    
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_lo
subplot(a, b, fig_index(3)); 
hold on; 
plot(D,   WT_met_reac_steady_1gL.prot(:,num_flux.lo),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.prot(:,num_flux.lo),   'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.lo})
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D,   WT_met_reac_steady_1gL.substrate(:,num_flux.lo),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.substrate(:,num_flux.lo),   'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.lo}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 


% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D,   WT_met_reac_steady_1gL.sig(:,num_flux.lo),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.sig(:,num_flux.lo),   'Color', plt_clrs.blue);
ylabel(label.met_reac.sig{num_flux.lo}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

%
% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on 
plot(D,   WT_met_reac_steady_1gL.atp(:,num_flux.lo),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.atp(:,num_flux.lo),   'Color', plt_clrs.blue);
hold off
ylabel(label.met_reac.atp{num_flux.lo}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 


% J_lp,cit/J_lp,fe
subplot(a, b, fig_index(1)); 
hold on; 
plot(D,   WT_met_reac_steady_1gL.flux(:,num_flux.lp_cit)./WT_met_reac_steady_1gL.flux(:,num_flux.lp_fe),...
'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.flux(:,num_flux.lp_cit)./WT_met_reac_steady_10gL.flux(:,num_flux.lp_fe),...
'Color', plt_clrs.blue);
ylabel('Jlp fe/cit')    
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

set(gcf,'Position' ,[10 140 912 658])




sgtitle('real metabolites')

%% virtual metabolites - precursor 

figure; 
a = 6+2;
b = 7;
fig_index = 1:6;

% precursor 

subplot(a, b, fig_index(1));
hold on; 
plot(D, WT_y_steady_1gL(:,num_y.pc), 'Color', plt_clrs.green);
plot(D, WT_y_steady_10gL(:,num_y.pc), 'Color', plt_clrs.blue);
ylabel(label.met(num_met.pc))    
xlabel(x_label)
ymax1 = max(WT_y_steady_1gL(:,num_y.pc)); 
ymax2 = max(WT_y_steady_10gL(:,num_y.pc)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on;
hold off; 

% J_gy
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.flux(:,num_flux.gy), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.gy), 'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.gy})    
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.flux(:,num_flux.gy)); 
ymax2 = max(WT_met_reac_steady_10gL.flux(:,num_flux.gy)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_gy
subplot(a, b, fig_index(3)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.prot(:,num_flux.gy), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.prot(:,num_flux.gy), 'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.gy})
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.prot(:,num_flux.gy)); 
ymax2 = max(WT_met_reac_steady_10gL.prot(:,num_flux.gy)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Gl/(km+Gl) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.substrate(:,num_flux.gy), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.substrate(:,num_flux.gy), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.gy}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.substrate(:,num_flux.gy)); 
ymax2 = max(WT_met_reac_steady_10gL.substrate(:,num_flux.gy)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.sig(:,num_flux.gy), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.sig(:,num_flux.gy), 'Color', plt_clrs.blue);
ylabel(label.met_reac.sig{num_flux.gy}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.sig(:,num_flux.gy)); 
ymax2 = max(WT_met_reac_steady_10gL.sig(:,num_flux.gy)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.atp(:,num_flux.gy), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.atp(:,num_flux.gy), 'Color', plt_clrs.blue);
ylabel(label.met_reac.atp{num_flux.gy}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.atp(:,num_flux.gy)); 
ymax2 = max(WT_met_reac_steady_10gL.atp(:,num_flux.gy)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

fig_index = fig_index + b; 
%--------------------------------------------------------------------------

% J_gn
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.flux(:,num_flux.gn), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.gn), 'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.gn})    
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.flux(:,num_flux.gn)); 
ymax2 = max(WT_met_reac_steady_10gL.flux(:,num_flux.gn)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_gn
subplot(a, b, fig_index(3)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.prot(:,num_flux.gn), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.prot(:,num_flux.gn), 'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.gn})
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.prot(:,num_flux.gn)); 
ymax2 = max(WT_met_reac_steady_10gL.prot(:,num_flux.gn)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Eh/(km+Eh) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.substrate(:,num_flux.gn), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.substrate(:,num_flux.gn), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.gn}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.substrate(:,num_flux.gn)); 
ymax2 = max(WT_met_reac_steady_10gL.substrate(:,num_flux.gn)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.sig(:,num_flux.gn), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.sig(:,num_flux.gn), 'Color', plt_clrs.blue);
ylabel(label.met_reac.sig{num_flux.gn}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.sig(:,num_flux.gn)); 
ymax2 = max(WT_met_reac_steady_10gL.sig(:,num_flux.gn)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold on; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.atp(:,num_flux.gn), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.atp(:,num_flux.gn), 'Color', plt_clrs.blue);
ylabel(label.met_reac.atp{num_flux.gn}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.atp(:,num_flux.gn)); 
ymax2 = max(WT_met_reac_steady_10gL.atp(:,num_flux.gn)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

fig_index = fig_index + b; 
%--------------------------------------------------------------------------

% J_mt
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.flux(:,num_flux.mt), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.mt), 'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.mt})    
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.flux(:,num_flux.mt)); 
ymax2 = max(WT_met_reac_steady_10gL.flux(:,num_flux.mt)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_mt
subplot(a, b, fig_index(3)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.prot(:,num_flux.mt), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.prot(:,num_flux.mt), 'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.mt})
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.prot(:,num_flux.mt)); 
ymax2 = max(WT_met_reac_steady_10gL.prot(:,num_flux.mt)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.substrate(:,num_flux.mt), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.substrate(:,num_flux.mt), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.mt}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.substrate(:,num_flux.mt)); 
ymax2 = max(WT_met_reac_steady_10gL.substrate(:,num_flux.mt)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.sig(:,num_flux.mt), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.sig(:,num_flux.mt), 'Color', plt_clrs.blue);
ylabel(label.met_reac.sig{num_flux.mt}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.sig(:,num_flux.mt)); 
ymax2 = max(WT_met_reac_steady_10gL.sig(:,num_flux.mt)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.atp(:,num_flux.mt), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.atp(:,num_flux.mt), 'Color', plt_clrs.blue);
ylabel(label.met_reac.atp{num_flux.mt}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.atp(:,num_flux.mt)); 
ymax2 = max(WT_met_reac_steady_10gL.atp(:,num_flux.mt)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

fig_index = fig_index + b; 
%--------------------------------------------------------------------------

% J_fe 
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.flux(:,num_flux.fe), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.fe})    
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.flux(:,num_flux.fe)); 
ymax2 = max(WT_met_reac_steady_10gL.flux(:,num_flux.fe)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_fe
subplot(a, b, fig_index(3)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.prot(:,num_flux.fe), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.prot(:,num_flux.fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.fe})
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.prot(:,num_flux.fe)); 
ymax2 = max(WT_met_reac_steady_10gL.prot(:,num_flux.fe)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.substrate(:,num_flux.fe), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.substrate(:,num_flux.fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.fe}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.substrate(:,num_flux.fe)); 
ymax2 = max(WT_met_reac_steady_10gL.substrate(:,num_flux.fe)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

%{
% signaling contribution 
subplot(a, b, fig_index(5)); 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.sig(:,num_flux.fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.sig{num_flux.fe}) 
xlabel(x_label)
ylim([0 1.05.*max(bart_WT_batch_min_met_reac.sig(:,num_flux.fe))])
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.atp(:,num_flux.fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.atp{num_flux.fe}) 
xlabel(x_label)
ylim([0 1.05.*max(bart_WT_batch_min_met_reac.atp(:,num_flux.fe))])
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
%}
fig_index = fig_index + b; 

%--------------------------------------------------------------------------

%{
% precursor
subplot(a, b, fig_index(1));
plot(bart_WT_batch_min_t, y(:,3), 'Color', plt_clrs.blue);
ylabel(label.met(3))    
xlabel(x_label)
axis square; 
box on;
%}

% J_as
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, 2*WT_met_reac_steady_1gL.flux(:,num_flux.as), 'Color', plt_clrs.green);
plot(D, 2*WT_met_reac_steady_10gL.flux(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.as})    
xlabel(x_label) 
ymax1 = max(2*WT_met_reac_steady_1gL.flux(:,num_flux.as)); 
ymax2 = max(2*WT_met_reac_steady_10gL.flux(:,num_flux.as)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_as
subplot(a, b, fig_index(3)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.prot(:,num_flux.as), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.prot(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.as})
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.prot(:,num_flux.as)); 
ymax2 = max(WT_met_reac_steady_10gL.prot(:,num_flux.as)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.substrate(:,num_flux.as), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.substrate(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.as}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.substrate(:,num_flux.as)); 
ymax2 = max(WT_met_reac_steady_10gL.substrate(:,num_flux.as)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.sig(:,num_flux.as), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.sig(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.sig{num_flux.as}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.sig(:,num_flux.as)); 
ymax2 = max(WT_met_reac_steady_10gL.sig(:,num_flux.as)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on;
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.atp(:,num_flux.as), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.atp(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.atp{num_flux.as}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.atp(:,num_flux.as)); 
ymax2 = max(WT_met_reac_steady_10gL.atp(:,num_flux.as)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% fig_index = fig_index + b; 
% J_po
subplot(a, b, fig_index(1)); 
hold on; 
plot(D,   WT_other_met_reac_steady_1gL.flux.po,   'Color', plt_clrs.green);
plot(D,   WT_other_met_reac_steady_10gL.flux.po,   'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_other_met_reac.flux.po, 'Color', plt_clrs.red);
ylabel('J_{po}')    
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold on; 

% -------------------------------------------------------------------------

fig_index = fig_index + b; 
% J_lp,fe
subplot(a, b, fig_index(2)); 
hold on; 
plot(D,   WT_met_reac_steady_1gL.flux(:,num_flux.lp_fe),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.flux(:,num_flux.lp_fe),   'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.flux(:,num_flux.lp_fe), 'Color', plt_clrs.red);
ylabel(label.met_reac.flux{num_flux.lp_fe})    
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_lp
subplot(a, b, fig_index(3)); 
hold on; 
plot(D,   WT_met_reac_steady_1gL.prot(:,num_flux.lp_fe),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.prot(:,num_flux.lp_fe),   'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.lp_fe), 'Color', plt_clrs.red);
ylabel(label.met_reac.prot{num_flux.lp_fe})
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D,   WT_met_reac_steady_1gL.substrate(:,num_flux.lp_fe),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.substrate(:,num_flux.lp_fe),   'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.lp_fe), 'Color', plt_clrs.red);
ylabel(label.met_reac.substrate{num_flux.lp_fe}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 


% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D,   WT_met_reac_steady_1gL.sig(:,num_flux.lp_fe),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.sig(:,num_flux.lp_fe),   'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.lp_fe), 'Color', plt_clrs.red);
ylabel(label.met_reac.sig{num_flux.lp_fe}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

%
% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on 
plot(D,   WT_met_reac_steady_1gL.atp(:,num_flux.lp_fe),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.atp(:,num_flux.lp_fe),   'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.lp_fe), 'Color', plt_clrs.red);
hold off
ylabel(label.met_reac.atp{num_flux.lp_fe}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 


% J_lp,cit
fig_index = fig_index + b; 
subplot(a, b, fig_index(2)); 
hold on; 
plot(D,   WT_met_reac_steady_1gL.flux(:,num_flux.lp_cit),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.flux(:,num_flux.lp_cit),   'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.flux(:,num_flux.lp_cit), 'Color', plt_clrs.red);
ylabel(label.met_reac.flux{num_flux.lp_cit})    
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_lp
subplot(a, b, fig_index(3)); 
hold on; 
plot(D,   WT_met_reac_steady_1gL.prot(:,num_flux.lp_cit),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.prot(:,num_flux.lp_cit),   'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.lp_cit), 'Color', plt_clrs.red);
ylabel(label.met_reac.prot{num_flux.lp_cit})
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D,   WT_met_reac_steady_1gL.substrate(:,num_flux.lp_cit),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.substrate(:,num_flux.lp_cit),   'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.lp_cit), 'Color', plt_clrs.red);
ylabel(label.met_reac.substrate{num_flux.lp_cit}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 


% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D,   WT_met_reac_steady_1gL.sig(:,num_flux.lp_cit),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.sig(:,num_flux.lp_cit),   'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.lp_cit), 'Color', plt_clrs.red);
ylabel(label.met_reac.sig{num_flux.lp_cit}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

%
% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on 
plot(D,   WT_met_reac_steady_1gL.atp(:,num_flux.lp_cit),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.atp(:,num_flux.lp_cit),   'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.lp_cit), 'Color', plt_clrs.red);
hold off
ylabel(label.met_reac.atp{num_flux.lp_cit}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 


fig_index = fig_index + b; 

%--------------------------------------------------------------------------

% J_lo 
subplot(a, b, fig_index(2)); 
hold on; 
plot(D,   WT_met_reac_steady_1gL.flux(:,num_flux.lo),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.flux(:,num_flux.lo),   'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.flux(:,num_flux.lo), 'Color', plt_clrs.red);
ylabel(label.met_reac.flux{num_flux.lo})    
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_lo
subplot(a, b, fig_index(3)); 
hold on; 
plot(D,   WT_met_reac_steady_1gL.prot(:,num_flux.lo),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.prot(:,num_flux.lo),   'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.lo), 'Color', plt_clrs.red);
ylabel(label.met_reac.prot{num_flux.lo})
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D,   WT_met_reac_steady_1gL.substrate(:,num_flux.lo),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.substrate(:,num_flux.lo),   'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.lo), 'Color', plt_clrs.red);
ylabel(label.met_reac.substrate{num_flux.lo}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 


% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D,   WT_met_reac_steady_1gL.sig(:,num_flux.lo),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.sig(:,num_flux.lo),   'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.lo), 'Color', plt_clrs.red);
ylabel(label.met_reac.sig{num_flux.lo}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

%
% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on 
plot(D,   WT_met_reac_steady_1gL.atp(:,num_flux.lo),   'Color', plt_clrs.green);
plot(D,   WT_met_reac_steady_10gL.atp(:,num_flux.lo),   'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.lo), 'Color', plt_clrs.red);
hold off
ylabel(label.met_reac.atp{num_flux.lo}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 

sgtitle('precursor')
set(gcf,'Position' ,[440 189 914 609])



%% virtual metabolites - ATP

figure; 
a = 6+1;
b = 6;
fig_index = 1:6;

%--------------------------------------------------------------------------
%                                   ATP
%--------------------------------------------------------------------------

% ATP
subplot(a, b, fig_index(1));
hold on; 
plot(D, WT_y_steady_1gL(:,num_y.ae), 'Color', plt_clrs.green);
plot(D, WT_y_steady_10gL(:,num_y.ae), 'Color', plt_clrs.blue);
ylabel(label.met{num_met.ae})    
xlabel(x_label)
ymax1 = max(WT_y_steady_1gL(:,num_y.ae)); 
ymax2 = max(WT_y_steady_10gL(:,num_y.ae)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on;
hold off; 
%}

% qgy * J_gy 
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, (par.q_gy).*(WT_met_reac_steady_1gL.flux(:,num_flux.gy)), 'Color', plt_clrs.green);
plot(D, (par.q_gy).*(WT_met_reac_steady_10gL.flux(:,num_flux.gy)), 'Color', plt_clrs.blue);
ylabel(label.atp_flux(1))    
xlabel(x_label)
ymax1 = max((par.q_gy).*WT_met_reac_steady_1gL.flux(:,num_flux.gy)); 
ymax2 = max((par.q_gy).*WT_met_reac_steady_10gL.flux(:,num_flux.gy)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_gy
subplot(a, b, fig_index(3)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.prot(:,num_flux.gy), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.prot(:,num_flux.gy), 'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.gy})
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.prot(:,num_flux.gy)); 
ymax2 = max(WT_met_reac_steady_10gL.prot(:,num_flux.gy)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Gl/(km+Gl) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.substrate(:,num_flux.gy), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.substrate(:,num_flux.gy), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.gy}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.substrate(:,num_flux.gy)); 
ymax2 = max(WT_met_reac_steady_10gL.substrate(:,num_flux.gy)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.sig(:,num_flux.gy), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.sig(:,num_flux.gy), 'Color', plt_clrs.blue);
ylabel(label.met_reac.sig{num_flux.gy}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.sig(:,num_flux.gy)); 
ymax2 = max(WT_met_reac_steady_10gL.sig(:,num_flux.gy)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.atp(:,num_flux.gy), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.atp(:,num_flux.gy), 'Color', plt_clrs.blue);
ylabel(label.met_reac.atp{num_flux.gy}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.atp(:,num_flux.gy)); 
ymax2 = max(WT_met_reac_steady_10gL.atp(:,num_flux.gy)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

fig_index = fig_index + b; 

%--------------------------------------------------------------------------

% qmt * J_mt
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, (par.q_mt).*(WT_met_reac_steady_1gL.flux(:,num_flux.mt)), 'Color', plt_clrs.green);
plot(D, (par.q_mt).*(WT_met_reac_steady_10gL.flux(:,num_flux.mt)), 'Color', plt_clrs.blue);
ylabel(label.atp_flux(2))    
xlabel(x_label)
ymax1 = max((par.q_mt).*WT_met_reac_steady_1gL.flux(:,num_flux.mt)); 
ymax2 = max((par.q_mt).*WT_met_reac_steady_10gL.flux(:,num_flux.mt)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_mt
subplot(a, b, fig_index(3)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.prot(:,num_flux.mt), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.prot(:,num_flux.mt), 'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.mt})
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.prot(:,num_flux.mt)); 
ymax2 = max(WT_met_reac_steady_10gL.prot(:,num_flux.mt)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.substrate(:,num_flux.mt), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.substrate(:,num_flux.mt), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.mt}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.substrate(:,num_flux.mt)); 
ymax2 = max(WT_met_reac_steady_10gL.substrate(:,num_flux.mt)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.sig(:,num_flux.mt), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.sig(:,num_flux.mt), 'Color', plt_clrs.blue);
ylabel(label.met_reac.sig{num_flux.mt}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.sig(:,num_flux.mt)); 
ymax2 = max(WT_met_reac_steady_10gL.sig(:,num_flux.mt)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.atp(:,num_flux.mt), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.atp(:,num_flux.mt), 'Color', plt_clrs.blue);
ylabel(label.met_reac.atp{num_flux.mt}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.atp(:,num_flux.mt)); 
ymax2 = max(WT_met_reac_steady_10gL.atp(:,num_flux.mt)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

fig_index = fig_index + b; 

%--------------------------------------------------------------------------

% q_gn * J_gn 
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, (par.q_gn).*(WT_met_reac_steady_1gL.flux(:,num_flux.gn)), 'Color', plt_clrs.green);
plot(D, (par.q_gn).*(WT_met_reac_steady_10gL.flux(:,num_flux.gn)), 'Color', plt_clrs.blue);
ylabel(label.atp_flux(3))    
xlabel(x_label)
ymax1 = max((par.q_gn).*WT_met_reac_steady_1gL.flux(:,num_flux.gn)); 
ymax2 = max((par.q_gn).*WT_met_reac_steady_10gL.flux(:,num_flux.gn)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_gn
subplot(a, b, fig_index(3)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.prot(:,num_flux.gn), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.prot(:,num_flux.gn), 'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.gn})
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.prot(:,num_flux.gn)); 
ymax2 = max(WT_met_reac_steady_10gL.prot(:,num_flux.gn)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Eh/(km+Eh) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.substrate(:,num_flux.gn), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.substrate(:,num_flux.gn), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.gn}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.substrate(:,num_flux.gn)); 
ymax2 = max(WT_met_reac_steady_10gL.substrate(:,num_flux.gn)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.sig(:,num_flux.gn), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.sig(:,num_flux.gn), 'Color', plt_clrs.blue);
ylabel(label.met_reac.sig{num_flux.gn}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.sig(:,num_flux.gn)); 
ymax2 = max(WT_met_reac_steady_10gL.sig(:,num_flux.gn)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

%{
% ATP contribution 
subplot(a, b, fig_index(6)); 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.atp(:,num_flux.gn), 'Color', plt_clrs.blue);
ylabel(label.met_reac.atp{num_flux.gn}) 
xlabel(x_label)
ylim([0 1.05.*max(bart_WT_batch_min_met_reac.atp(:,num_flux.gn))])
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
%}

fig_index = fig_index + b; 
%--------------------------------------------------------------------------

% qas * J_as
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, (par.q_as).*(WT_met_reac_steady_1gL.flux(:,num_flux.as)), 'Color', plt_clrs.green);
plot(D, (par.q_as).*(WT_met_reac_steady_10gL.flux(:,num_flux.as)), 'Color', plt_clrs.blue);
ylabel(label.atp_flux(4))    
xlabel(x_label)
ymax1 = max((par.q_as).*WT_met_reac_steady_1gL.flux(:,num_flux.as)); 
ymax2 = max((par.q_as).*WT_met_reac_steady_10gL.flux(:,num_flux.as)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_as
subplot(a, b, fig_index(3)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.prot(:,num_flux.as), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.prot(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.as})
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.prot(:,num_flux.as)); 
ymax2 = max(WT_met_reac_steady_10gL.prot(:,num_flux.as)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.substrate(:,num_flux.as), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.substrate(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.as}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.substrate(:,num_flux.as)); 
ymax2 = max(WT_met_reac_steady_10gL.substrate(:,num_flux.as)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.sig(:,num_flux.as), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.sig(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.sig{num_flux.as}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.sig(:,num_flux.as)); 
ymax2 = max(WT_met_reac_steady_10gL.sig(:,num_flux.as)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.atp(:,num_flux.as), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.atp(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.atp{num_flux.as}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.atp(:,num_flux.as)); 
ymax2 = max(WT_met_reac_steady_10gL.atp(:,num_flux.as)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

fig_index = fig_index + b; 
%--------------------------------------------------------------------------

% qp * protein synthesis rate
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, par.q_p.*WT_prot_syn_steady_1gL.prot_syn_rate(:), 'Color', plt_clrs.green);
plot(D, par.q_p.*WT_prot_syn_steady_10gL.prot_syn_rate(:), 'Color', plt_clrs.blue);
ylabel(label.atp_flux(5)) 
xlabel(x_label)
ymax1 = max(par.q_p.*WT_prot_syn_steady_1gL.prot_syn_rate(:)); 
ymax2 = max(par.q_p.*WT_prot_syn_steady_10gL.prot_syn_rate(:)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% free bart_WT_batch_min_ribosomes
subplot(a, b, fig_index(3)); 
hold on; 
plot(D, WT_y_steady_1gL(:,num_y.r0), 'Color', plt_clrs.green);
plot(D, WT_y_steady_10gL(:,num_y.r0), 'Color', plt_clrs.blue);
ylabel(label.rib_cells{1}) 
xlabel(x_label)
ymax1 = max(WT_y_steady_1gL(:,num_y.r0)); 
ymax2 = max(WT_y_steady_10gL(:,num_y.r0)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% tc
subplot(a, b, fig_index(4)); 
hold on; 
plot(D, WT_prot_syn_steady_1gL.tc(:), 'Color', plt_clrs.green);
plot(D, WT_prot_syn_steady_10gL.tc(:), 'Color', plt_clrs.blue);
ylabel(label.prot_syn.tc) 
xlabel(x_label)
ymax1 = max(WT_prot_syn_steady_1gL.tc(:)); 
ymax2 = max(WT_prot_syn_steady_10gL.tc(:)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% snf1 contribution
subplot(a, b, fig_index(5)); 
hold on; 
plot(D, WT_prot_syn_steady_1gL.eIF_a_s(:), 'Color', plt_clrs.green);
plot(D, WT_prot_syn_steady_10gL.eIF_a_s(:), 'Color', plt_clrs.blue);
ylabel(label.prot_syn.eIF_a_s) 
xlabel(x_label)
ymax1 = max(WT_prot_syn_steady_1gL.eIF_a_s(:)); 
ymax2 = max(WT_prot_syn_steady_10gL.eIF_a_s(:)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% tor contribution
subplot(a, b, fig_index(6)); 
hold on; 
plot(D, WT_prot_syn_steady_1gL.eIF_a_tau(:), 'Color', plt_clrs.green);
plot(D, WT_prot_syn_steady_10gL.eIF_a_tau(:), 'Color', plt_clrs.blue);
ylabel(label.prot_syn.eIF_a_tau) 
xlabel(x_label)
ymax1 = max(WT_prot_syn_steady_1gL.eIF_a_tau(:)); 
ymax2 = max(WT_prot_syn_steady_10gL.eIF_a_tau(:)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

%--------------------------------------------------------------------------

% J_eo
subplot(a, b, fig_index(1)); 
hold on; 
plot(D, WT_other_met_reac_steady_1gL.flux.eo, 'Color', plt_clrs.green);
plot(D, WT_other_met_reac_steady_10gL.flux.eo, 'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_other_met_reac.flux.eo, 'Color', plt_clrs.red);
ylabel(label.atp_flux(6)) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 



fig_index = fig_index + b; 


%--------------------------------------------------------------------------

% par.q_lo * J_lo

subplot(a, b, fig_index(2)); 
hold on; 
plot(D, (par.q_lo).*(WT_met_reac_steady_1gL.flux(:,num_flux.lo)), 'Color', plt_clrs.green);
plot(D, (par.q_lo).*(WT_met_reac_steady_10gL.flux(:,num_flux.lo)), 'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, (par.q_lo).*(murphy_WT_batch_YPD_met_reac.flux(:,num_flux.lo)), 'Color', plt_clrs.red);
ylabel(label.atp_flux(7))    
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_lo
subplot(a, b, fig_index(3)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.prot(:,num_flux.lo), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.prot(:,num_flux.lo), 'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.lo), 'Color', plt_clrs.red);
ylabel(label.met_reac.prot{num_flux.lo})
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.substrate(:,num_flux.lo), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.substrate(:,num_flux.lo), 'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.lo), 'Color', plt_clrs.red);
ylabel(label.met_reac.substrate{num_flux.lo}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.sig(:,num_flux.lo), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.sig(:,num_flux.lo), 'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.lo), 'Color', plt_clrs.red);
ylabel(label.met_reac.sig{num_flux.lo}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.atp(:,num_flux.lo), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.atp(:,num_flux.lo), 'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.lo), 'Color', plt_clrs.red);
ylabel(label.met_reac.atp{num_flux.lo}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

fig_index = fig_index + b; 


% par.q_lp * J_lp
subplot(a, b, fig_index(2)); 
hold on; 
plot(D,   (par.q_lp).*(WT_met_reac_steady_1gL.flux(:,num_flux.lp_fe)   + WT_met_reac_steady_1gL.flux(:,num_flux.lp_cit)), 'Color', plt_clrs.green);
plot(D,   (par.q_lp).*(WT_met_reac_steady_10gL.flux(:,num_flux.lp_fe)   + WT_met_reac_steady_10gL.flux(:,num_flux.lp_cit)), 'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, (par.q_lp).*(murphy_WT_batch_YPD_met_reac.flux(:,num_flux.lp_fe) + murphy_WT_batch_YPD_met_reac.flux(:,num_flux.lp_cit)), 'Color', plt_clrs.red);
ylabel(label.atp_flux(8))    
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_lp_cit
subplot(a, b, fig_index(3)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.prot(:,num_flux.lp_cit), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.prot(:,num_flux.lp_cit), 'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.lp_cit), 'Color', plt_clrs.red);
ylabel(label.met_reac.prot{num_flux.lp_cit})
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.substrate(:,num_flux.lp_cit), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.substrate(:,num_flux.lp_cit), 'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.lp_cit), 'Color', plt_clrs.red);
ylabel(label.met_reac.substrate{num_flux.lp_cit}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.sig(:,num_flux.lp_cit), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.sig(:,num_flux.lp_cit), 'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.lp_cit), 'Color', plt_clrs.red);
ylabel(label.met_reac.sig{num_flux.lp_cit}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.atp(:,num_flux.lp_cit), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.atp(:,num_flux.lp_cit), 'Color', plt_clrs.blue);
%plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.lp_cit), 'Color', plt_clrs.red);
ylabel(label.met_reac.atp{num_flux.as}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 
sgtitle('ATP')
set(gcf,'Position' ,[70 67 885 738])


%% virtual metabolites - amino acids 

figure; 
a = 5;
b = 6;
fig_index = 1:6;

%--------------------------------------------------------------------------
%                       intracellular amino acids
%--------------------------------------------------------------------------

% intracellular amino acids 
subplot(a, b, fig_index(1));
hold on; 
plot(D, WT_y_steady_1gL(:,num_y.aa_in), 'Color', plt_clrs.green);
plot(D, WT_y_steady_10gL(:,num_y.aa_in), 'Color', plt_clrs.blue);
ylabel(label.met{num_met.aa_in})    
xlabel(x_label)
ymax1 = max(WT_y_steady_1gL(:,num_y.aa_in)); 
ymax2 = max(WT_y_steady_10gL(:,num_y.aa_in)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on;
hold off; 

% J_as
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.flux(:,num_flux.as), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.as})    
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.flux(:,num_flux.as)); 
ymax2 = max(WT_met_reac_steady_10gL.flux(:,num_flux.as)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_as
subplot(a, b, fig_index(3)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.prot(:,num_flux.as), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.prot(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.as})
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.prot(:,num_flux.as)); 
ymax2 = max(WT_met_reac_steady_10gL.prot(:,num_flux.as)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.substrate(:,num_flux.as), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.substrate(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.as}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.substrate(:,num_flux.as)); 
ymax2 = max(WT_met_reac_steady_10gL.substrate(:,num_flux.as)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.sig(:,num_flux.as), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.sig(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.sig{num_flux.as}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.sig(:,num_flux.as)); 
ymax2 = max(WT_met_reac_steady_10gL.sig(:,num_flux.as)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.atp(:,num_flux.as), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.atp(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.atp{num_flux.as}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.atp(:,num_flux.as)); 
ymax2 = max(WT_met_reac_steady_10gL.atp(:,num_flux.as)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

fig_index = fig_index + b; 

%--------------------------------------------------------------------------

% protein synthesis rate
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, WT_prot_syn_steady_1gL.prot_syn_rate(:), 'Color', plt_clrs.green);
plot(D, WT_prot_syn_steady_10gL.prot_syn_rate(:), 'Color', plt_clrs.blue);
ylabel(label.prot_syn.prot_syn_rate) 
xlabel(x_label)
ymax1 = max(WT_prot_syn_steady_1gL.prot_syn_rate(:)); 
ymax2 = max(WT_prot_syn_steady_10gL.prot_syn_rate(:)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% Rat
subplot(a, b, fig_index(3)); 
hold on; 
plot(D, WT_rib_steady_1gL.rat, 'Color', plt_clrs.green);
plot(D, WT_rib_steady_10gL.rat, 'Color', plt_clrs.blue);
ylabel('R_{at}') 
xlabel(x_label)
ymax1 = max(WT_rib_steady_1gL.rat); 
ymax2 = max(WT_rib_steady_10gL.rat); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

%{
% tc
subplot(a, b, fig_index(4)); 
plot(bart_WT_batch_min_t, bart_WT_batch_min_prot_syn.tc(:), 'Color', plt_clrs.blue);
ylabel(label.bart_WT_batch_min_prot_syn.tc) 
xlabel(x_label)
ylim([0 1.05.*max(bart_WT_batch_min_prot_syn.tc(:))])
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 

% snf1 contribution
subplot(a, b, fig_index(5)); 
plot(bart_WT_batch_min_t, bart_WT_batch_min_prot_syn.eIF_a_s(:), 'Color', plt_clrs.blue);
ylabel(label.bart_WT_batch_min_prot_syn.eIF_a_s) 
xlabel(x_label)
ylim([0 1.05.*max(bart_WT_batch_min_prot_syn.eIF_a_s(:))])
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 

% tor contribution
subplot(a, b, fig_index(6)); 
plot(bart_WT_batch_min_t, bart_WT_batch_min_prot_syn.eIF_a_tau(:), 'Color', plt_clrs.blue);
ylabel(label.bart_WT_batch_min_prot_syn.eIF_a_tau) 
xlabel(x_label)
ylim([0 1.05.*max(bart_WT_batch_min_prot_syn.eIF_a_tau(:))])
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
%}

fig_index = fig_index + b; 
%--------------------------------------------------------------------------
%                       extracellular amino acids
%--------------------------------------------------------------------------

% extracellular amino acids 
subplot(a, b, fig_index(1));
hold on; 
plot(D, WT_y_steady_1gL(:,num_y.aaex), 'Color', plt_clrs.green);
plot(D, WT_y_steady_10gL(:,num_y.aaex), 'Color', plt_clrs.blue);
ylabel(label.met{1})    
xlabel(x_label)
ymax1 = max(WT_y_steady_1gL(:,num_y.aaex)); 
ymax2 = max(WT_y_steady_10gL(:,num_y.aaex)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on;
hold off; 

% J_at
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.flux(:,num_flux.at), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.at), 'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.at})    
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.flux(:,num_flux.at)); 
ymax2 = max(WT_met_reac_steady_10gL.flux(:,num_flux.at)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_at
subplot(a, b, fig_index(3)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.prot(:,num_flux.at), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.prot(:,num_flux.at), 'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.at})
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.prot(:,num_flux.at)); 
ymax2 = max(WT_met_reac_steady_10gL.prot(:,num_flux.at)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Aa_ex/(km+Aa_ex) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.substrate(:,num_flux.at), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.substrate(:,num_flux.at), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.at}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.substrate(:,num_flux.at)); 
ymax2 = max(WT_met_reac_steady_10gL.substrate(:,num_flux.at)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.sig(:,num_flux.at), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.sig(:,num_flux.at), 'Color', plt_clrs.blue);
ylabel(label.met_reac.sig{num_flux.at}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.sig(:,num_flux.at)); 
ymax2 = max(WT_met_reac_steady_10gL.sig(:,num_flux.at)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.atp(:,num_flux.at), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.atp(:,num_flux.at), 'Color', plt_clrs.blue);
ylabel(label.met_reac.atp{num_flux.at}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.atp(:,num_flux.at)); 
ymax2 = max(WT_met_reac_steady_10gL.atp(:,num_flux.at)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

%--------------------------------------------------------------------------
%                                     NH4
%--------------------------------------------------------------------------


fig_index = fig_index + b; 
% extracellular NH4 
subplot(a, b, fig_index(1));
hold on; 
plot(D, WT_y_steady_1gL(:,num_y.nh4), 'Color', plt_clrs.green);
plot(D, WT_y_steady_10gL(:,num_y.nh4), 'Color', plt_clrs.blue);
ylabel('NH4')    
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on;
hold off; 


% J_nh4
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, par.n_nh* WT_met_reac_steady_1gL.flux(:,num_flux.as), 'Color', plt_clrs.green);
plot(D, par.n_nh* WT_met_reac_steady_10gL.flux(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel('J_{NH4}')    
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_as
subplot(a, b, fig_index(3)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.prot(:,num_flux.as), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.prot(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.as})
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Aa_ex/(km+Aa_ex) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.substrate(:,num_flux.as), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.substrate(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.as}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.sig(:,num_flux.as), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.sig(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.sig{num_flux.as}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.atp(:,num_flux.as), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.atp(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.atp{num_flux.as}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% -------------------------------------------------------------------------


fig_index = fig_index + b; 
% Jno, NH4 to others (others = other than Jas)
subplot(a, b, fig_index(1));
hold on; 
plot(D, WT_other_met_reac_steady_1gL.flux.no, 'Color', plt_clrs.green);
plot(D, WT_other_met_reac_steady_10gL.flux.no, 'Color', plt_clrs.blue);
ylabel('J_{no}')    
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on;
hold off; 
sgtitle('amino acids')

set(gcf,'Position' ,[384 209 1017 589])

%% ------------------------------------------------------------------------
%                     other growth related intermediate
%--------------------------------------------------------------------------
% other intermediate module_gene_expression_update_ribosome
figure; 
for i = 1: numel(label.others)
subplot(ceil(numel(label.others)/4),4,i)
hold on; 
plot(D, WT_prot_syn_steady_1gL.rib.others(:,i), 'Color', plt_clrs.green);
plot(D, WT_prot_syn_steady_10gL.rib.others(:,i), 'Color', plt_clrs.blue);
hold off; 
ylabel(label.others{i}, 'Interpreter', 'none')    
xlabel(x_label)
axis square; 
box on;
ylim([0 max(ylim)*1.05])
end 
sgtitle('other growth intermediates')
set(gcf,'Position' , [682 67 605 731])

%% partition

figure; 
a = 4;  b = 4;
subplot(a, b, 1+2*b); 
hold on; 
plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.gy), 'Color', plt_clrs.blue);
plot(kumar_hg.dr, kumar_hg.Jgy, 'o', 'Color', plt_clrs.green)
ylabel('Jgy')    
xlabel(x_label) 
ylim([0 max(ylim)*1.05])
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% J_as
subplot(a, b, 1); 
hold on; 
plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.as), 'Color', plt_clrs.blue);
%plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.mt), 'Color', plt_clrs.red);
ylabel('Jas')    
xlabel(x_label) 
ylim([0 max(ylim)*1.05])
%legend('Jmt')
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% J_mt
subplot(a, b, 1+b); 
hold on; 
%plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.as), 'Color', plt_clrs.blue);
plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.mt), 'Color', plt_clrs.red);
ylabel('Jmt')    
xlabel(x_label) 
ylim([0 max(ylim)*1.05])
%legend('Jmt')
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

subplot(a, b, 2+b); 
hold on; 
%plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.as), 'Color', plt_clrs.blue);
plot(kumar_hg.dr, kumar_hg.o2_uptake, 'o-', 'Color', plt_clrs.red);
plot(kumar_lg.dr, kumar_lg.o2_uptake, 'o-', 'Color', plt_clrs.blue);
ylabel('O2 uptake')    
xlabel(x_label) 
ylim([0 max(ylim)*1.05])
xlim([0 0.42])
%legend('Jmt')
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% subplot(a, b, 1+b); 
% hold on; 
% plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.as), 'Color', plt_clrs.blue);
% plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.mt), 'Color', plt_clrs.red);
% ylabel('Jas,mt')    
% xlabel(x_label) 
% ylim([0 max(ylim)*1.05])
% legend('Jas','Jmt')
% %xline(t_gl_depl)
% %xline(t_eh_depl)
% axis square; 
% box on; 
% hold off; 

subplot(a, b, 1+3*b); 
hold on; 
plot(D, 4*WT_met_reac_steady_10gL.flux(:,num_flux.as)'./(2*WT_met_reac_steady_10gL.flux(:,num_flux.gy)'), 'Color', plt_clrs.blue);
plot(D, (4*WT_met_reac_steady_10gL.flux(:,num_flux.as)'+WT_met_reac_steady_10gL.flux(:,num_flux.mt)')./(2*WT_met_reac_steady_10gL.flux(:,num_flux.gy)'), 'Color', plt_clrs.red);
plot(D, (4*WT_met_reac_steady_10gL.flux(:,num_flux.as)'+WT_met_reac_steady_10gL.flux(:,num_flux.mt)' + WT_met_reac_steady_10gL.flux(:,num_flux.fe)')./(2*WT_met_reac_steady_10gL.flux(:,num_flux.gy)'), 'Color', plt_clrs.green);

ylabel('ratio')    
xlabel(x_label) 
ylim([0 max(ylim)*1.05])
legend('4*Jas/2*Jgy','4*Jas+Jmt/2*Jgy','4*Jas+Jmt +Jfe/2*Jgy')
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

subplot(a, b, 3); 
hold on; 
plot(D, (par.q_gy*WT_met_reac_steady_10gL.flux(:,num_flux.gy)')./(par.q_gy*WT_met_reac_steady_10gL.flux(:,num_flux.gy)'+ par.q_mt*WT_met_reac_steady_10gL.flux(:,num_flux.mt)'), 'Color', plt_clrs.blue);
plot(D, (par.q_mt*WT_met_reac_steady_10gL.flux(:,num_flux.mt)')./(par.q_gy*WT_met_reac_steady_10gL.flux(:,num_flux.gy)'+ par.q_mt*WT_met_reac_steady_10gL.flux(:,num_flux.mt)'), 'Color', plt_clrs.red);
%plot(kumar_hg.dr, kumar_hg.Jgy, 'o', 'Color', plt_clrs.green)
ylabel('atp prod')   
legend('Jgy% ','Jmt% ')
xlabel(x_label) 
ylim([0 max(ylim)*1.05])
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

subplot(a, b, 3+b); 
hold on; 
plot(D, ((par.q_p + par.q_as + par.q_eo/par.pro_den).*WT_met_reac_steady_10gL.flux(:,num_flux.as)')./(par.q_gy*WT_met_reac_steady_10gL.flux(:,num_flux.gy)'+ par.q_mt*WT_met_reac_steady_10gL.flux(:,num_flux.mt)'), 'Color', plt_clrs.blue);
%plot(kumar_hg.dr, kumar_hg.Jgy, 'o', 'Color', plt_clrs.green)
ylabel('ratio')   
legend('(q_p  + (q_eo/pro_den) + q_as) *Jas vs Jatp prod')
xlabel(x_label) 
ylim([0 max(ylim)*1.05])
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

%gl_in_1gL  = 1.0  * gl_gPerL_to_uM; 
gl_in_10gL = 10.0 * gl_gPerL_to_uM; 
%gl_in_08gL = 0.8  * gl_gPerL_to_uM; 
%biomass_yield_1gL  = (WT_y_steady_1gL(:,end)  .* (1/gdw_to_cell))./((gl_in_1gL*(ones(size(WT_y_steady_1gL(:,num_y.gl),1),1))) - ((WT_y_steady_1gL(:,num_y.gl) .* 1/gl_gPerL_to_uM)));
biomass_yield_10gL = (WT_y_steady_10gL(:,end) .* (1/gdw_to_cell))./((gl_in_10gL*(ones(size(WT_y_steady_10gL(:,num_y.gl),1),1))) - ((WT_y_steady_10gL(:,num_y.gl) .* 1/gl_gPerL_to_uM)));
%biomass_yield_08gL = (WT_y_steady_08gL(:,end) .* (1/gdw_to_cell))./((gl_in_08gL*(ones(size(WT_y_steady_08gL(:,num_y.gl),1),1))) - ((WT_y_steady_08gL(:,num_y.gl) .* 1/gl_gPerL_to_uM)));

% biomass yield/Jmt
subplot(a, b, 4+2*b); 
hold on; 
%plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.as), 'Color', plt_clrs.blue);
plot(WT_met_reac_steady_10gL.flux(:,num_flux.mt), biomass_yield_10gL./WT_met_reac_steady_10gL.flux(:,num_flux.mt), 'Color', plt_clrs.red);
ylabel('biomass yield / J_{mt}')    
xlabel('J_{mt}') 
%ylim([0 max(ylim)*1.05])
%legend('Jmt')
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

%%
%{
2*WT_met_reac_steady_10gL.flux(:,num_flux.gy)'
2*WT_met_reac_steady_10gL.flux(:,num_flux.as)' + WT_met_reac_steady_10gL.flux(:,num_flux.mt)' + 2*WT_met_reac_steady_10gL.flux(:,num_flux.as)'
4*WT_met_reac_steady_10gL.flux(:,num_flux.as)'./WT_met_reac_steady_10gL.flux(:,num_flux.mt)'
WT_met_reac_steady_10gL.flux(:,num_flux.mt)'./WT_met_reac_steady_10gL.flux(:,num_flux.as)'
% 2*((par.q_p + par.q_as + par.q_eo/par.pro_den)- par.q_gy)/(par.q_gy +2*par.q_mt) % Jmt /Jas 
4*WT_met_reac_steady_10gL.flux(:,num_flux.as)'./(2*WT_met_reac_steady_10gL.flux(:,num_flux.gy)')

(par.q_gy*WT_met_reac_steady_10gL.flux(:,num_flux.gy)' + par.q_mt*WT_met_reac_steady_10gL.flux(:,num_flux.mt)')./WT_met_reac_steady_10gL.flux(:,num_flux.as)'
(par.q_p + par.q_as + par.q_eo/par.pro_den)*WT_met_reac_steady_10gL.flux(:,num_flux.as)'

2*WT_met_reac_steady_10gL.flux(:,num_flux.sp)' + WT_met_reac_steady_10gL.flux(:,num_flux.sd)'
WT_met_reac_steady_10gL.flux(:,num_flux.fe)'

2*WT_met_reac_steady_10gL.flux(:,num_flux.as)' + WT_met_reac_steady_10gL.flux(:,num_flux.mt)' + WT_met_reac_steady_10gL.flux(:,num_flux.fe)'
%}

%     kumar high gl Jgy 
%     0.2653    0.5591    1.9427    3.5632
%     new par v2 update:
%     0.2863    0.6623    2.6581    3.5067    
%     old par: 
%      0.3747    0.8213    1.8874    3.7680 --> lead to low cell number

%     kumar high gl Jgy 
% 0.5976    0.9971    3.2883    4.5551 --> CO2
% 0.1992    0.3324    1.0961    1.5184 --> Jmt
% 2* Jgy - Jmt = 0.2863*2 - 0.1992  = 93350 (Jas)
% 93350/0.12 = 7.7792e+05 (should be cell density)

% (approx.) Jgy equation before Dcr: 
% (0.12* par.pro_den)*((par.q_p + par.q_as + par.q_eo/par.pro_den) + 4*par.q_mt)/(par.q_gy + 2*par.q_mt)



%% virtual metabolites - storage carbon 

figure; 
a = 2;
b = 6;
fig_index = 1:6;

%--------------------------------------------------------------------------
%                       storage carbon
%--------------------------------------------------------------------------

% storage carbon
subplot(a, b, fig_index(1));
hold on; 
plot(D, WT_y_steady_1gL(:,num_y.sc), 'Color', plt_clrs.green);
plot(D, WT_y_steady_10gL(:,num_y.sc), 'Color', plt_clrs.blue);
ylabel(label.met{num_met.sc})    
xlabel(x_label)
ymax1 = max(WT_y_steady_1gL(:,num_y.sc)); 
ymax2 = max(WT_y_steady_10gL(:,num_y.sc)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on;
hold off; 

% J_sp
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.flux(:,num_flux.sp), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.sp), 'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.sp})    
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.flux(:,num_flux.sp)); 
ymax2 = max(WT_met_reac_steady_10gL.flux(:,num_flux.sp)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_sp
subplot(a, b, fig_index(3)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.prot(:,num_flux.sp), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.prot(:,num_flux.sp), 'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.sp})
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.prot(:,num_flux.sp)); 
ymax2 = max(WT_met_reac_steady_10gL.prot(:,num_flux.sp)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.substrate(:,num_flux.sp), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.substrate(:,num_flux.sp), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.sp}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.substrate(:,num_flux.sp)); 
ymax2 = max(WT_met_reac_steady_10gL.substrate(:,num_flux.sp)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.sig(:,num_flux.sp), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.sig(:,num_flux.sp), 'Color', plt_clrs.blue);
ylabel(label.met_reac.sig{num_flux.sp}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.sig(:,num_flux.sp)); 
ymax2 = max(WT_met_reac_steady_10gL.sig(:,num_flux.sp)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.atp(:,num_flux.sp), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.atp(:,num_flux.sp), 'Color', plt_clrs.blue);
ylabel(label.met_reac.atp{num_flux.sp}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.atp(:,num_flux.sp)); 
ymax2 = max(WT_met_reac_steady_10gL.atp(:,num_flux.sp)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

fig_index = fig_index + b; 

%--------------------------------------------------------------------------

% J_sd
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.flux(:,num_flux.sd), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.sd), 'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.sd})    
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.flux(:,num_flux.sd)); 
ymax2 = max(WT_met_reac_steady_10gL.flux(:,num_flux.sd)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_sd
subplot(a, b, fig_index(3)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.prot(:,num_flux.sd), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.prot(:,num_flux.sd), 'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.sd})
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.prot(:,num_flux.sd)); 
ymax2 = max(WT_met_reac_steady_10gL.prot(:,num_flux.sd)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.substrate(:,num_flux.sd), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.substrate(:,num_flux.sd), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.sd}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.substrate(:,num_flux.sd)); 
ymax2 = max(WT_met_reac_steady_10gL.substrate(:,num_flux.sd)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.sig(:,num_flux.sd), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.sig(:,num_flux.sd), 'Color', plt_clrs.blue);
ylabel(label.met_reac.sig{num_flux.sd}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.sig(:,num_flux.sd)); 
ymax2 = max(WT_met_reac_steady_10gL.sig(:,num_flux.sd)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.atp(:,num_flux.sd), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.atp(:,num_flux.sd), 'Color', plt_clrs.blue);
ylabel(label.met_reac.atp{num_flux.sd}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.atp(:,num_flux.sd)); 
ymax2 = max(WT_met_reac_steady_10gL.atp(:,num_flux.sd)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 


%fig_index = fig_index + b; 


end 
