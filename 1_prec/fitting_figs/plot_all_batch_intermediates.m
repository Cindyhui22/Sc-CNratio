function plot_all_batch_intermediates(bart_WT_batch_min_t, zamp_WT_batch_min_t, murphy_WT_batch_YPD_t, ...
bart_WT_batch_min_y, zamp_WT_batch_min_y, murphy_WT_batch_YPD_y, ...
bart_WT_batch_min_met_reac, zamp_WT_batch_min_met_reac, murphy_WT_batch_YPD_met_reac, ...
bart_WT_batch_min_other_met_reac, zamp_WT_batch_min_other_met_reac, murphy_WT_batch_YPD_other_met_reac, ...
bart_WT_batch_min_prot_syn, zamp_WT_batch_min_prot_syn, murphy_WT_batch_YPD_prot_syn, ...
bart_WT_batch_min_rib, zamp_WT_batch_min_rib, murphy_WT_batch_YPD_rib, ...
label, plt_clrs, par, num_y, num_flux,num_prot,  x_label)

% set(0, 'DefaultLineLineWidth',  1.5);
% set(0, 'DefaultLineMarkerSize', 3);
plot_num;

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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_y(:,num_y.gl),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_y(:,num_y.gl),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y(:,num_y.gl), 'Color', plt_clrs.red);
ylabel(label.met{2})    
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_y(:,num_y.gl)); 
ymax2 = max(zamp_WT_batch_min_y(:,num_y.gl)); 
ymax3 = max(murphy_WT_batch_YPD_y(:,num_y.gl)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on;
hold off; 

% J_gy
subplot(a, b, fig_index(2)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.flux(:,num_flux.gy),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.flux(:,num_flux.gy),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.flux(:,num_flux.gy), 'Color', plt_clrs.red);
ylabel(label.met_reac.flux{num_flux.gy})    
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.flux(:,num_flux.gy)); 
ymax2 = max(zamp_WT_batch_min_met_reac.flux(:,num_flux.gy)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.flux(:,num_flux.gy)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_gy
subplot(a, b, fig_index(3)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.prot(:,num_flux.gy),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.prot(:,num_flux.gy),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.gy), 'Color', plt_clrs.red);
ylabel(label.met_reac.prot{num_flux.gy})
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.prot(:,num_flux.gy)); 
ymax2 = max(zamp_WT_batch_min_met_reac.prot(:,num_flux.gy)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.prot(:,num_flux.gy)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Gl/(km+Gl) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.substrate(:,num_flux.gy),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.substrate(:,num_flux.gy),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.gy), 'Color', plt_clrs.red);
ylabel(label.met_reac.substrate{num_flux.gy}); 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.substrate(:,num_flux.gy)); 
ymax2 = max(zamp_WT_batch_min_met_reac.substrate(:,num_flux.gy)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.gy)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.sig(:,num_flux.gy),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.sig(:,num_flux.gy),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.gy), 'Color', plt_clrs.red);
ylabel(label.met_reac.sig{num_flux.gy}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.sig(:,num_flux.gy)); 
ymax2 = max(zamp_WT_batch_min_met_reac.sig(:,num_flux.gy)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.sig(:,num_flux.gy)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.atp(:,num_flux.gy),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.atp(:,num_flux.gy),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.gy), 'Color', plt_clrs.red);
ylabel(label.met_reac.atp{num_flux.gy}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.atp(:,num_flux.gy)); 
ymax2 = max(zamp_WT_batch_min_met_reac.atp(:,num_flux.gy)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.atp(:,num_flux.gy)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_y(:,num_y.eh),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_y(:,num_y.eh),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y(:,num_y.eh), 'Color', plt_clrs.red);
ylabel(label.met{num_met.eh})    
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_y(:,num_y.eh)); 
ymax2 = max(zamp_WT_batch_min_y(:,num_y.eh)); 
ymax3 = max(murphy_WT_batch_YPD_y(:,num_y.eh)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on;
hold off; 

% J_fe 
subplot(a, b, fig_index(2)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.flux(:,num_flux.fe),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.flux(:,num_flux.fe),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.flux(:,num_flux.fe), 'Color', plt_clrs.red);
ylabel(label.met_reac.flux{num_flux.fe})    
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.flux(:,num_flux.fe)); 
ymax2 = max(zamp_WT_batch_min_met_reac.flux(:,num_flux.fe)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.flux(:,num_flux.fe)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_fe
subplot(a, b, fig_index(3)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.prot(:,num_flux.fe),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.prot(:,num_flux.fe),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.fe), 'Color', plt_clrs.red);
ylabel(label.met_reac.prot{num_flux.fe})
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.prot(:,num_flux.fe)); 
ymax2 = max(zamp_WT_batch_min_met_reac.prot(:,num_flux.fe)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.prot(:,num_flux.fe)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.substrate(:,num_flux.fe),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.substrate(:,num_flux.fe),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.fe), 'Color', plt_clrs.red);
ylabel(label.met_reac.substrate{num_flux.fe}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.substrate(:,num_flux.fe)); 
ymax2 = max(zamp_WT_batch_min_met_reac.substrate(:,num_flux.fe)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.fe)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.flux(:,num_flux.gn),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.flux(:,num_flux.gn),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.flux(:,num_flux.gn), 'Color', plt_clrs.red);
ylabel('J_{gn}')    
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.flux(:,num_flux.gn)); 
ymax2 = max(zamp_WT_batch_min_met_reac.flux(:,num_flux.gn)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.flux(:,num_flux.gn)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_gn
subplot(a, b, fig_index(3)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.prot(:,num_flux.gn),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.prot(:,num_flux.gn),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.gn), 'Color', plt_clrs.red);
ylabel(label.met_reac.prot{num_flux.gn})
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.prot(:,num_flux.gn)); 
ymax2 = max(zamp_WT_batch_min_met_reac.prot(:,num_flux.gn)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.prot(:,num_flux.gn)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Eh/(km+Eh) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.substrate(:,num_flux.gn),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.substrate(:,num_flux.gn),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.gn), 'Color', plt_clrs.red);
ylabel(label.met_reac.substrate{num_flux.gn}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.substrate(:,num_flux.gn)); 
ymax2 = max(zamp_WT_batch_min_met_reac.substrate(:,num_flux.gn)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.gn)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.sig(:,num_flux.gn),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.sig(:,num_flux.gn),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.gn), 'Color', plt_clrs.red);
ylabel(label.met_reac.sig{num_flux.gn}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.sig(:,num_flux.gn)); 
ymax2 = max(zamp_WT_batch_min_met_reac.sig(:,num_flux.gn)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.sig(:,num_flux.gn)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.atp(:,num_flux.gn),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.atp(:,num_flux.gn),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.gn), 'Color', plt_clrs.red);
hold off
ylabel(label.met_reac.atp{num_flux.gn}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.atp(:,num_flux.gn)); 
ymax2 = max(zamp_WT_batch_min_met_reac.atp(:,num_flux.gn)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.atp(:,num_flux.gn)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
%}

fig_index = fig_index + b; 


%--------------------------------------------------------------------------
%                              Lipid
%--------------------------------------------------------------------------
% lipid
subplot(a, b, fig_index(1));
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_y(:,num_y.lp_e),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_y(:,num_y.lp_e),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y(:,num_y.lp_e), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.flux(:,num_flux.lp_fe),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.flux(:,num_flux.lp_fe),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.flux(:,num_flux.lp_fe), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.prot(:,num_flux.lp_fe),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.prot(:,num_flux.lp_fe),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.lp_fe), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.substrate(:,num_flux.lp_fe),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.substrate(:,num_flux.lp_fe),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.lp_fe), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.sig(:,num_flux.lp_fe),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.sig(:,num_flux.lp_fe),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.lp_fe), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.atp(:,num_flux.lp_fe),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.atp(:,num_flux.lp_fe),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.lp_fe), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.flux(:,num_flux.lp_cit),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.flux(:,num_flux.lp_cit),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.flux(:,num_flux.lp_cit), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.prot(:,num_flux.lp_cit),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.prot(:,num_flux.lp_cit),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.lp_cit), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.substrate(:,num_flux.lp_cit),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.substrate(:,num_flux.lp_cit),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.lp_cit), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.sig(:,num_flux.lp_cit),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.sig(:,num_flux.lp_cit),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.lp_cit), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.atp(:,num_flux.lp_cit),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.atp(:,num_flux.lp_cit),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.lp_cit), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.flux(:,num_flux.lo),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.flux(:,num_flux.lo),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.flux(:,num_flux.lo), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.prot(:,num_flux.lo),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.prot(:,num_flux.lo),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.lo), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.substrate(:,num_flux.lo),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.substrate(:,num_flux.lo),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.lo), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.sig(:,num_flux.lo),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.sig(:,num_flux.lo),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.lo), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.atp(:,num_flux.lo),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.atp(:,num_flux.lo),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.lo), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.flux(:,num_flux.lp_cit)./bart_WT_batch_min_met_reac.flux(:,num_flux.lp_fe),...
'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.flux(:,num_flux.lp_cit)./zamp_WT_batch_min_met_reac.flux(:,num_flux.lp_fe),...
'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.flux(:,num_flux.lp_cit)./murphy_WT_batch_YPD_met_reac.flux(:,num_flux.lp_fe),...
'Color', plt_clrs.red);
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
a = 5+2;
b = 6;
fig_index = 1:6;

% precursor 

subplot(a, b, fig_index(1));
hold on;
plot(bart_WT_batch_min_t,   bart_WT_batch_min_y(:,num_y.pc),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_y(:,num_y.pc),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y(:,num_y.pc), 'Color', plt_clrs.red);
hold off; 
ylabel(label.met(num_met.pc))    
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_y(:,num_y.pc)); 
ymax2 = max(zamp_WT_batch_min_y(:,num_y.pc)); 
ymax3 = max(murphy_WT_batch_YPD_y(:,num_y.pc)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on;
hold off; 

% J_gy
subplot(a, b, fig_index(2)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.flux(:,num_flux.gy),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.flux(:,num_flux.gy),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.flux(:,num_flux.gy), 'Color', plt_clrs.red);
ylabel(label.met_reac.flux{num_flux.gy})    
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.flux(:,num_flux.gy)); 
ymax2 = max(zamp_WT_batch_min_met_reac.flux(:,num_flux.gy)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.flux(:,num_flux.gy)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_gy
subplot(a, b, fig_index(3)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.prot(:,num_flux.gy),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.prot(:,num_flux.gy),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.gy), 'Color', plt_clrs.red);
ylabel(label.met_reac.prot{num_flux.gy})
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.prot(:,num_flux.gy)); 
ymax2 = max(zamp_WT_batch_min_met_reac.prot(:,num_flux.gy)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.prot(:,num_flux.gy)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Gl/(km+Gl) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.substrate(:,num_flux.gy),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.substrate(:,num_flux.gy),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.gy), 'Color', plt_clrs.red);
ylabel(label.met_reac.substrate{num_flux.gy}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.substrate(:,num_flux.gy)); 
ymax2 = max(zamp_WT_batch_min_met_reac.substrate(:,num_flux.gy)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.gy)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.sig(:,num_flux.gy),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.sig(:,num_flux.gy),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.gy), 'Color', plt_clrs.red);
ylabel(label.met_reac.sig{num_flux.gy}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.sig(:,num_flux.gy)); 
ymax2 = max(zamp_WT_batch_min_met_reac.sig(:,num_flux.gy)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.sig(:,num_flux.gy)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.atp(:,num_flux.gy),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.atp(:,num_flux.gy),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.gy), 'Color', plt_clrs.red);
ylabel(label.met_reac.atp{num_flux.gy}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.atp(:,num_flux.gy)); 
ymax2 = max(zamp_WT_batch_min_met_reac.atp(:,num_flux.gy)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.atp(:,num_flux.gy)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.flux(:,num_flux.mt),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.flux(:,num_flux.mt),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.flux(:,num_flux.mt), 'Color', plt_clrs.red);
ylabel(label.met_reac.flux{num_flux.mt})    
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.flux(:,num_flux.mt)); 
ymax2 = max(zamp_WT_batch_min_met_reac.flux(:,num_flux.mt)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.flux(:,num_flux.mt)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_mt
subplot(a, b, fig_index(3)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.prot(:,num_flux.mt),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.prot(:,num_flux.mt),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.mt), 'Color', plt_clrs.red);
ylabel(label.met_reac.prot{num_flux.mt})
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.prot(:,num_flux.mt)); 
ymax2 = max(zamp_WT_batch_min_met_reac.prot(:,num_flux.mt)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.prot(:,num_flux.mt)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.substrate(:,num_flux.mt),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.substrate(:,num_flux.mt),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.mt), 'Color', plt_clrs.red);
ylabel(label.met_reac.substrate{num_flux.mt}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.substrate(:,num_flux.mt)); 
ymax2 = max(zamp_WT_batch_min_met_reac.substrate(:,num_flux.mt)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.mt)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.sig(:,num_flux.mt),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.sig(:,num_flux.mt),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.mt), 'Color', plt_clrs.red);
ylabel(label.met_reac.sig{num_flux.mt}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.sig(:,num_flux.mt)); 
ymax2 = max(zamp_WT_batch_min_met_reac.sig(:,num_flux.mt)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.sig(:,num_flux.mt)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.atp(:,num_flux.mt),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.atp(:,num_flux.mt),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.mt), 'Color', plt_clrs.red);
ylabel(label.met_reac.atp{num_flux.mt}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.atp(:,num_flux.mt)); 
ymax2 = max(zamp_WT_batch_min_met_reac.atp(:,num_flux.mt)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.atp(:,num_flux.mt)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

fig_index = fig_index + b; 
%--------------------------------------------------------------------------

% ethanol yield
subplot(a, b, fig_index(1));
hold on; 
plot(bart_WT_batch_min_t,   (bart_WT_batch_min_y(:,num_y.eh) - bart_WT_batch_min_y(1,num_y.eh))./(bart_WT_batch_min_y(1,num_y.gl) - bart_WT_batch_min_y(:,num_y.gl)),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   (zamp_WT_batch_min_y(:,num_y.eh) - zamp_WT_batch_min_y(1,num_y.eh))./(zamp_WT_batch_min_y(1,num_y.gl) - zamp_WT_batch_min_y(:,num_y.gl)),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, (murphy_WT_batch_YPD_y(:,num_y.eh) - murphy_WT_batch_YPD_y(1,num_y.eh))./(murphy_WT_batch_YPD_y(1,num_y.gl) - murphy_WT_batch_YPD_y(:,num_y.gl)), 'Color', plt_clrs.red);
ylabel('ethanol yield')    
xlabel(x_label)
%{
ymax1 = max(bart_WT_batch_min_met_reac.flux(:,num_flux.fe)); 
ymax2 = max(zamp_WT_batch_min_met_reac.flux(:,num_flux.fe)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.flux(:,num_flux.fe)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%}
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 


% J_fe 
subplot(a, b, fig_index(2)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.flux(:,num_flux.fe),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.flux(:,num_flux.fe),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.flux(:,num_flux.fe), 'Color', plt_clrs.red);
ylabel(label.met_reac.flux{num_flux.fe})    
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.flux(:,num_flux.fe)); 
ymax2 = max(zamp_WT_batch_min_met_reac.flux(:,num_flux.fe)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.flux(:,num_flux.fe)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_fe
subplot(a, b, fig_index(3)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.prot(:,num_flux.fe),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.prot(:,num_flux.fe),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.fe), 'Color', plt_clrs.red);
ylabel(label.met_reac.prot{num_flux.fe})
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.prot(:,num_flux.fe)); 
ymax2 = max(zamp_WT_batch_min_met_reac.prot(:,num_flux.fe)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.prot(:,num_flux.fe)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.substrate(:,num_flux.fe),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.substrate(:,num_flux.fe),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.fe), 'Color', plt_clrs.red);
ylabel(label.met_reac.substrate{num_flux.fe}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.substrate(:,num_flux.fe)); 
ymax2 = max(zamp_WT_batch_min_met_reac.substrate(:,num_flux.fe)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.fe)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.flux(:,num_flux.as),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.flux(:,num_flux.as),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.flux(:,num_flux.as), 'Color', plt_clrs.red);
ylabel(label.met_reac.flux{num_flux.as})    
xlabel(x_label) 
ymax1 = max(bart_WT_batch_min_met_reac.flux(:,num_flux.as)); 
ymax2 = max(zamp_WT_batch_min_met_reac.flux(:,num_flux.as)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.flux(:,num_flux.as)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_as
subplot(a, b, fig_index(3)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.prot(:,num_flux.as),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.prot(:,num_flux.as),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.as), 'Color', plt_clrs.red);
ylabel(label.met_reac.prot{num_flux.as})
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.prot(:,num_flux.as)); 
ymax2 = max(zamp_WT_batch_min_met_reac.prot(:,num_flux.as)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.prot(:,num_flux.as)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.substrate(:,num_flux.as),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.substrate(:,num_flux.as),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.as), 'Color', plt_clrs.red);
ylabel(label.met_reac.substrate{num_flux.as}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.substrate(:,num_flux.as)); 
ymax2 = max(zamp_WT_batch_min_met_reac.substrate(:,num_flux.as)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.as)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.sig(:,num_flux.as),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.sig(:,num_flux.as),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.as), 'Color', plt_clrs.red);
ylabel(label.met_reac.sig{num_flux.as}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.sig(:,num_flux.as)); 
ymax2 = max(zamp_WT_batch_min_met_reac.sig(:,num_flux.as)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.sig(:,num_flux.as)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on;
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.atp(:,num_flux.as),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.atp(:,num_flux.as),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.as), 'Color', plt_clrs.red);
ylabel(label.met_reac.atp{num_flux.as}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.atp(:,num_flux.as)); 
ymax2 = max(zamp_WT_batch_min_met_reac.atp(:,num_flux.as)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.atp(:,num_flux.as)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 



%--------------------------------------------------------------------------
% fig_index = fig_index + b; 
% J_po
subplot(a, b, fig_index(1)); 
hold on; 
plot(bart_WT_batch_min_t,   bart_WT_batch_min_other_met_reac.flux.po,   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_other_met_reac.flux.po,   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_other_met_reac.flux.po, 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.flux(:,num_flux.lp_fe),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.flux(:,num_flux.lp_fe),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.flux(:,num_flux.lp_fe), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.prot(:,num_flux.lp_fe),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.prot(:,num_flux.lp_fe),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.lp_fe), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.substrate(:,num_flux.lp_fe),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.substrate(:,num_flux.lp_fe),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.lp_fe), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.sig(:,num_flux.lp_fe),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.sig(:,num_flux.lp_fe),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.lp_fe), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.atp(:,num_flux.lp_fe),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.atp(:,num_flux.lp_fe),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.lp_fe), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.flux(:,num_flux.lp_cit),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.flux(:,num_flux.lp_cit),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.flux(:,num_flux.lp_cit), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.prot(:,num_flux.lp_cit),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.prot(:,num_flux.lp_cit),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.lp_cit), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.substrate(:,num_flux.lp_cit),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.substrate(:,num_flux.lp_cit),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.lp_cit), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.sig(:,num_flux.lp_cit),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.sig(:,num_flux.lp_cit),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.lp_cit), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.atp(:,num_flux.lp_cit),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.atp(:,num_flux.lp_cit),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.lp_cit), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.flux(:,num_flux.lo),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.flux(:,num_flux.lo),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.flux(:,num_flux.lo), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.prot(:,num_flux.lo),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.prot(:,num_flux.lo),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.lo), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.substrate(:,num_flux.lo),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.substrate(:,num_flux.lo),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.lo), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.sig(:,num_flux.lo),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.sig(:,num_flux.lo),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.lo), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   bart_WT_batch_min_met_reac.atp(:,num_flux.lo),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   zamp_WT_batch_min_met_reac.atp(:,num_flux.lo),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.lo), 'Color', plt_clrs.red);
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

%}


%% virtual metabolites - ATP

figure; 
a = 7;
b = 6;
fig_index = 1:6;

%--------------------------------------------------------------------------
%                                   ATP
%--------------------------------------------------------------------------

% ATP
subplot(a, b, fig_index(1));
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_y(:,num_y.ae), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_y(:,num_y.ae), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y(:,num_y.ae), 'Color', plt_clrs.red);
hold off;
ylabel(label.met{num_met.ae})    
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_y(:,num_y.ae)); 
ymax2 = max(zamp_WT_batch_min_y(:,num_y.ae)); 
ymax3 = max(murphy_WT_batch_YPD_y(:,num_y.ae)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on;

%}

% qgy * J_gy 
subplot(a, b, fig_index(2)); 
hold on; 
plot(bart_WT_batch_min_t,   (par.q_gy).*(bart_WT_batch_min_met_reac.flux(:,num_flux.gy)),   'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   (par.q_gy).*(zamp_WT_batch_min_met_reac.flux(:,num_flux.gy)),   'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, (par.q_gy).*(murphy_WT_batch_YPD_met_reac.flux(:,num_flux.gy)), 'Color', plt_clrs.red);
ylabel(label.atp_flux(1))    
xlabel(x_label)
ymax1 = max((par.q_gy).*bart_WT_batch_min_met_reac.flux(:,num_flux.gy)); 
ymax2 = max((par.q_gy).*zamp_WT_batch_min_met_reac.flux(:,num_flux.gy)); 
ymax3 = max((par.q_gy).*murphy_WT_batch_YPD_met_reac.flux(:,num_flux.gy)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_gy
subplot(a, b, fig_index(3)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.prot(:,num_flux.gy), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.prot(:,num_flux.gy), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.gy), 'Color', plt_clrs.red);
ylabel(label.met_reac.prot{num_flux.gy})
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.prot(:,num_flux.gy)); 
ymax2 = max(zamp_WT_batch_min_met_reac.prot(:,num_flux.gy)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.prot(:,num_flux.gy)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Gl/(km+Gl) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.substrate(:,num_flux.gy), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.substrate(:,num_flux.gy), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.gy), 'Color', plt_clrs.red);
ylabel(label.met_reac.substrate{num_flux.gy}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.substrate(:,num_flux.gy)); 
ymax2 = max(zamp_WT_batch_min_met_reac.substrate(:,num_flux.gy)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.gy)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.sig(:,num_flux.gy), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.sig(:,num_flux.gy), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.gy), 'Color', plt_clrs.red);
ylabel(label.met_reac.sig{num_flux.gy}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.sig(:,num_flux.gy)); 
ymax2 = max(zamp_WT_batch_min_met_reac.sig(:,num_flux.gy)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.sig(:,num_flux.gy)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.atp(:,num_flux.gy), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.atp(:,num_flux.gy), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.gy), 'Color', plt_clrs.red);
ylabel(label.met_reac.atp{num_flux.gy}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.atp(:,num_flux.gy)); 
ymax2 = max(zamp_WT_batch_min_met_reac.atp(:,num_flux.gy)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.atp(:,num_flux.gy)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
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
plot(bart_WT_batch_min_t, (par.q_mt).*(bart_WT_batch_min_met_reac.flux(:,num_flux.mt)), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, (par.q_mt).*(zamp_WT_batch_min_met_reac.flux(:,num_flux.mt)), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, (par.q_mt).*(murphy_WT_batch_YPD_met_reac.flux(:,num_flux.mt)), 'Color', plt_clrs.red);
ylabel(label.atp_flux(2))    
xlabel(x_label)
ymax1 = max((par.q_mt).*bart_WT_batch_min_met_reac.flux(:,num_flux.mt)); 
ymax2 = max((par.q_mt).*zamp_WT_batch_min_met_reac.flux(:,num_flux.mt)); 
ymax3 = max((par.q_mt).*murphy_WT_batch_YPD_met_reac.flux(:,num_flux.mt)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_mt
subplot(a, b, fig_index(3)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.prot(:,num_flux.mt), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.prot(:,num_flux.mt), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.mt), 'Color', plt_clrs.red);
ylabel(label.met_reac.prot{num_flux.mt})
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.prot(:,num_flux.mt)); 
ymax2 = max(zamp_WT_batch_min_met_reac.prot(:,num_flux.mt)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.prot(:,num_flux.mt)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.substrate(:,num_flux.mt), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.substrate(:,num_flux.mt), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.mt), 'Color', plt_clrs.red);
ylabel(label.met_reac.substrate{num_flux.mt}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.substrate(:,num_flux.mt)); 
ymax2 = max(zamp_WT_batch_min_met_reac.substrate(:,num_flux.mt)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.mt)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.sig(:,num_flux.mt), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.sig(:,num_flux.mt), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.mt), 'Color', plt_clrs.red);
ylabel(label.met_reac.sig{num_flux.mt}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.sig(:,num_flux.mt)); 
ymax2 = max(zamp_WT_batch_min_met_reac.sig(:,num_flux.mt)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.sig(:,num_flux.mt)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.atp(:,num_flux.mt), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.atp(:,num_flux.mt), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.mt), 'Color', plt_clrs.red);
ylabel(label.met_reac.atp{num_flux.mt}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.atp(:,num_flux.mt)); 
ymax2 = max(zamp_WT_batch_min_met_reac.atp(:,num_flux.mt)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.atp(:,num_flux.mt)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
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
plot(bart_WT_batch_min_t, (par.q_gn).*(bart_WT_batch_min_met_reac.flux(:,num_flux.gn)), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, (par.q_gn).*(zamp_WT_batch_min_met_reac.flux(:,num_flux.gn)), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, (par.q_gn).*(murphy_WT_batch_YPD_met_reac.flux(:,num_flux.gn)), 'Color', plt_clrs.red);
ylabel(label.atp_flux(3))    
xlabel(x_label)
ymax1 = max((par.q_gn).*bart_WT_batch_min_met_reac.flux(:,num_flux.gn)); 
ymax2 = max((par.q_gn).*zamp_WT_batch_min_met_reac.flux(:,num_flux.gn)); 
ymax3 = max((par.q_gn).*murphy_WT_batch_YPD_met_reac.flux(:,num_flux.gn)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_gn
subplot(a, b, fig_index(3)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.prot(:,num_flux.gn), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.prot(:,num_flux.gn), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.gn), 'Color', plt_clrs.red);
ylabel(label.met_reac.prot{num_flux.gn})
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.prot(:,num_flux.gn)); 
ymax2 = max(zamp_WT_batch_min_met_reac.prot(:,num_flux.gn)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.prot(:,num_flux.gn)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Eh/(km+Eh) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.substrate(:,num_flux.gn), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.substrate(:,num_flux.gn), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.gn), 'Color', plt_clrs.red);
ylabel(label.met_reac.substrate{num_flux.gn}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.substrate(:,num_flux.gn)); 
ymax2 = max(zamp_WT_batch_min_met_reac.substrate(:,num_flux.gn)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.gn)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.sig(:,num_flux.gn), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.sig(:,num_flux.gn), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.gn), 'Color', plt_clrs.red);
ylabel(label.met_reac.sig{num_flux.gn}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.sig(:,num_flux.gn)); 
ymax2 = max(zamp_WT_batch_min_met_reac.sig(:,num_flux.gn)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.sig(:,num_flux.gn)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
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
plot(bart_WT_batch_min_t, (par.q_as).*(bart_WT_batch_min_met_reac.flux(:,num_flux.as)), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, (par.q_as).*(zamp_WT_batch_min_met_reac.flux(:,num_flux.as)), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, (par.q_as).*(murphy_WT_batch_YPD_met_reac.flux(:,num_flux.as)), 'Color', plt_clrs.red);
ylabel(label.atp_flux(4))    
xlabel(x_label)
ymax1 = max(-(par.q_as).*bart_WT_batch_min_met_reac.flux(:,num_flux.as)); 
ymax2 = max(-(par.q_as).*zamp_WT_batch_min_met_reac.flux(:,num_flux.as)); 
ymax3 = max(-(par.q_as).*murphy_WT_batch_YPD_met_reac.flux(:,num_flux.as)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_as
subplot(a, b, fig_index(3)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.prot(:,num_flux.as), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.prot(:,num_flux.as), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.as), 'Color', plt_clrs.red);
ylabel(label.met_reac.prot{num_flux.as})
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.prot(:,num_flux.as)); 
ymax2 = max(zamp_WT_batch_min_met_reac.prot(:,num_flux.as)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.prot(:,num_flux.as)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.substrate(:,num_flux.as), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.substrate(:,num_flux.as), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.as), 'Color', plt_clrs.red);
ylabel(label.met_reac.substrate{num_flux.as}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.substrate(:,num_flux.as)); 
ymax2 = max(zamp_WT_batch_min_met_reac.substrate(:,num_flux.as)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.as)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.sig(:,num_flux.as), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.sig(:,num_flux.as), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.as), 'Color', plt_clrs.red);
ylabel(label.met_reac.sig{num_flux.as}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.sig(:,num_flux.as)); 
ymax2 = max(zamp_WT_batch_min_met_reac.sig(:,num_flux.as)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.sig(:,num_flux.as)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.atp(:,num_flux.as), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.atp(:,num_flux.as), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.as), 'Color', plt_clrs.red);
ylabel(label.met_reac.atp{num_flux.as}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.atp(:,num_flux.as)); 
ymax2 = max(zamp_WT_batch_min_met_reac.atp(:,num_flux.as)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.atp(:,num_flux.as)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
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
plot(bart_WT_batch_min_t, par.q_p.*bart_WT_batch_min_prot_syn.prot_syn_rate(:), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, par.q_p.*zamp_WT_batch_min_prot_syn.prot_syn_rate(:), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, par.q_p.*murphy_WT_batch_YPD_prot_syn.prot_syn_rate(:), 'Color', plt_clrs.red);
ylabel(label.atp_flux(5)) 
xlabel(x_label)
ymax1 = max(par.q_p.*bart_WT_batch_min_prot_syn.prot_syn_rate(:)); 
ymax2 = max(par.q_p.*zamp_WT_batch_min_prot_syn.prot_syn_rate(:)); 
ymax3 = max(par.q_p.*murphy_WT_batch_YPD_prot_syn.prot_syn_rate(:)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% free bart_WT_batch_min_ribosomes
subplot(a, b, fig_index(3)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_y(:,num_y.r0), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_y(:,num_y.r0), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y(:,num_y.r0), 'Color', plt_clrs.red);
ylabel(label.rib_cells{1}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_y(:,num_y.r0)); 
ymax2 = max(zamp_WT_batch_min_y(:,num_y.r0)); 
ymax3 = max(murphy_WT_batch_YPD_y(:,num_y.r0)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% tc
subplot(a, b, fig_index(4)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_prot_syn.tc(:), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_prot_syn.tc(:), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_prot_syn.tc(:), 'Color', plt_clrs.red);
ylabel(label.prot_syn.tc) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_prot_syn.tc(:)); 
ymax2 = max(zamp_WT_batch_min_prot_syn.tc(:)); 
ymax3 = max(murphy_WT_batch_YPD_prot_syn.tc(:)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% snf1 contribution
subplot(a, b, fig_index(5)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_prot_syn.eIF_a_s(:), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_prot_syn.eIF_a_s(:), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_prot_syn.eIF_a_s(:), 'Color', plt_clrs.red);
ylabel(label.prot_syn.eIF_a_s) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_prot_syn.eIF_a_s(:)); 
ymax2 = max(zamp_WT_batch_min_prot_syn.eIF_a_s(:)); 
ymax3 = max(murphy_WT_batch_YPD_prot_syn.eIF_a_s(:)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% tor contribution
subplot(a, b, fig_index(6)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_prot_syn.eIF_a_tau(:), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_prot_syn.eIF_a_tau(:), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_prot_syn.eIF_a_tau(:), 'Color', plt_clrs.red);
ylabel(label.prot_syn.eIF_a_tau) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_prot_syn.eIF_a_tau(:)); 
ymax2 = max(zamp_WT_batch_min_prot_syn.eIF_a_tau(:)); 
ymax3 = max(murphy_WT_batch_YPD_prot_syn.eIF_a_tau(:)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
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
plot(bart_WT_batch_min_t, bart_WT_batch_min_other_met_reac.flux.eo, 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_other_met_reac.flux.eo, 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_other_met_reac.flux.eo, 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t, (par.q_lo).*(bart_WT_batch_min_met_reac.flux(:,num_flux.lo)), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, (par.q_lo).*(zamp_WT_batch_min_met_reac.flux(:,num_flux.lo)), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, (par.q_lo).*(murphy_WT_batch_YPD_met_reac.flux(:,num_flux.lo)), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.prot(:,num_flux.lo), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.prot(:,num_flux.lo), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.lo), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.substrate(:,num_flux.lo), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.substrate(:,num_flux.lo), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.lo), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.sig(:,num_flux.lo), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.sig(:,num_flux.lo), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.lo), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.atp(:,num_flux.lo), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.atp(:,num_flux.lo), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.lo), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t,   (par.q_lp).*(bart_WT_batch_min_met_reac.flux(:,num_flux.lp_fe)   + bart_WT_batch_min_met_reac.flux(:,num_flux.lp_cit)), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t,   (par.q_lp).*(zamp_WT_batch_min_met_reac.flux(:,num_flux.lp_fe)   + zamp_WT_batch_min_met_reac.flux(:,num_flux.lp_cit)), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, (par.q_lp).*(murphy_WT_batch_YPD_met_reac.flux(:,num_flux.lp_fe) + murphy_WT_batch_YPD_met_reac.flux(:,num_flux.lp_cit)), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.prot(:,num_flux.lp_cit), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.prot(:,num_flux.lp_cit), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.lp_cit), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.substrate(:,num_flux.lp_cit), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.substrate(:,num_flux.lp_cit), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.lp_cit), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.sig(:,num_flux.lp_cit), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.sig(:,num_flux.lp_cit), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.lp_cit), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.atp(:,num_flux.lp_cit), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.atp(:,num_flux.lp_cit), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.lp_cit), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t, bart_WT_batch_min_y(:,num_y.aa_in), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_y(:,num_y.aa_in), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y(:,num_y.aa_in), 'Color', plt_clrs.red);
ylabel(label.met{num_met.aa_in})    
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_y(:,num_y.aa_in)); 
ymax2 = max(zamp_WT_batch_min_y(:,num_y.aa_in)); 
ymax3 = max(murphy_WT_batch_YPD_y(:,num_y.aa_in)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on;
hold off; 

% J_as
subplot(a, b, fig_index(2)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.flux(:,num_flux.as), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.flux(:,num_flux.as), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.flux(:,num_flux.as), 'Color', plt_clrs.red);
ylabel(label.met_reac.flux{num_flux.as})    
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.flux(:,num_flux.as)); 
ymax2 = max(zamp_WT_batch_min_met_reac.flux(:,num_flux.as)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.flux(:,num_flux.as)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_as
subplot(a, b, fig_index(3)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.prot(:,num_flux.as), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.prot(:,num_flux.as), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.as), 'Color', plt_clrs.red);
ylabel(label.met_reac.prot{num_flux.as})
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.prot(:,num_flux.as)); 
ymax2 = max(zamp_WT_batch_min_met_reac.prot(:,num_flux.as)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.prot(:,num_flux.as)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.substrate(:,num_flux.as), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.substrate(:,num_flux.as), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.as), 'Color', plt_clrs.red);
ylabel(label.met_reac.substrate{num_flux.as}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.substrate(:,num_flux.as)); 
ymax2 = max(zamp_WT_batch_min_met_reac.substrate(:,num_flux.as)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.as)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.sig(:,num_flux.as), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.sig(:,num_flux.as), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.as), 'Color', plt_clrs.red);
ylabel(label.met_reac.sig{num_flux.as}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.sig(:,num_flux.as)); 
ymax2 = max(zamp_WT_batch_min_met_reac.sig(:,num_flux.as)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.sig(:,num_flux.as)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.atp(:,num_flux.as), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.atp(:,num_flux.as), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.as), 'Color', plt_clrs.red);
ylabel(label.met_reac.atp{num_flux.as}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.atp(:,num_flux.as)); 
ymax2 = max(zamp_WT_batch_min_met_reac.atp(:,num_flux.as)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.atp(:,num_flux.as)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
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
plot(bart_WT_batch_min_t, bart_WT_batch_min_prot_syn.prot_syn_rate(:), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_prot_syn.prot_syn_rate(:), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_prot_syn.prot_syn_rate(:), 'Color', plt_clrs.red);
ylabel(label.prot_syn.prot_syn_rate) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_prot_syn.prot_syn_rate(:)); 
ymax2 = max(zamp_WT_batch_min_prot_syn.prot_syn_rate(:)); 
ymax3 = max(murphy_WT_batch_YPD_prot_syn.prot_syn_rate(:)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% Rat
subplot(a, b, fig_index(3)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_rib.rat, 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_rib.rat, 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_rib.rat, 'Color', plt_clrs.red);
ylabel('R_{at}') 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_rib.rat); 
ymax2 = max(zamp_WT_batch_min_rib.rat); 
ymax3 = max(murphy_WT_batch_YPD_rib.rat); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
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
plot(bart_WT_batch_min_t, bart_WT_batch_min_y(:,num_y.aaex), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_y(:,num_y.aaex), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y(:,num_y.aaex), 'Color', plt_clrs.red);
ylabel(label.met{1})    
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_y(:,num_y.aaex)); 
ymax2 = max(zamp_WT_batch_min_y(:,num_y.aaex)); 
ymax3 = max(murphy_WT_batch_YPD_y(:,num_y.aaex)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on;
hold off; 

% J_at
subplot(a, b, fig_index(2)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.flux(:,num_flux.at), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.flux(:,num_flux.at), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.flux(:,num_flux.at), 'Color', plt_clrs.red);
ylabel(label.met_reac.flux{num_flux.at})    
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.flux(:,num_flux.at)); 
ymax2 = max(zamp_WT_batch_min_met_reac.flux(:,num_flux.at)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.flux(:,num_flux.at)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_at
subplot(a, b, fig_index(3)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.prot(:,num_flux.at), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.prot(:,num_flux.at), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.at), 'Color', plt_clrs.red);
ylabel(label.met_reac.prot{num_flux.at})
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.prot(:,num_flux.at)); 
ymax2 = max(zamp_WT_batch_min_met_reac.prot(:,num_flux.at)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.prot(:,num_flux.at)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Aa_ex/(km+Aa_ex) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.substrate(:,num_flux.at), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.substrate(:,num_flux.at), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.at), 'Color', plt_clrs.red);
ylabel(label.met_reac.substrate{num_flux.at}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.substrate(:,num_flux.at)); 
ymax2 = max(zamp_WT_batch_min_met_reac.substrate(:,num_flux.at)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.at)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.sig(:,num_flux.at), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.sig(:,num_flux.at), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.at), 'Color', plt_clrs.red);
ylabel(label.met_reac.sig{num_flux.at}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.sig(:,num_flux.at)); 
ymax2 = max(zamp_WT_batch_min_met_reac.sig(:,num_flux.at)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.sig(:,num_flux.at)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.atp(:,num_flux.at), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.atp(:,num_flux.at), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.at), 'Color', plt_clrs.red);
ylabel(label.met_reac.atp{num_flux.at}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.atp(:,num_flux.at)); 
ymax2 = max(zamp_WT_batch_min_met_reac.atp(:,num_flux.at)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.atp(:,num_flux.at)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
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
plot(bart_WT_batch_min_t, bart_WT_batch_min_y(:,num_y.nh4), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_y(:,num_y.nh4), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y(:,num_y.nh4), 'Color', plt_clrs.red);
ylabel('NH4')    
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_y(:,num_y.nh4)); 
ymax2 = max(zamp_WT_batch_min_y(:,num_y.nh4)); 
ymax3 = max(murphy_WT_batch_YPD_y(:,num_y.nh4)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on;
hold off; 


% J_nh4
subplot(a, b, fig_index(2)); 
hold on; 
plot(bart_WT_batch_min_t, par.n_nh* bart_WT_batch_min_met_reac.flux(:,num_flux.as), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, par.n_nh* zamp_WT_batch_min_met_reac.flux(:,num_flux.as), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, par.n_nh* murphy_WT_batch_YPD_met_reac.flux(:,num_flux.as), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.prot(:,num_flux.as), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.prot(:,num_flux.as), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.as), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.substrate(:,num_flux.as), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.substrate(:,num_flux.as), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.as), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.sig(:,num_flux.as), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.sig(:,num_flux.as), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.as), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.atp(:,num_flux.as), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.atp(:,num_flux.as), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.as), 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t, bart_WT_batch_min_other_met_reac.flux.no, 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_other_met_reac.flux.no, 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_other_met_reac.flux.no, 'Color', plt_clrs.red);
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
plot(bart_WT_batch_min_t, bart_WT_batch_min_prot_syn.rib.others(:,i), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_prot_syn.rib.others(:,i), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_prot_syn.rib.others(:,i), 'Color', plt_clrs.red);
hold off; 
ylabel(label.others{i}, 'Interpreter', 'none')    
xlabel(x_label)
axis square; 
box on;
ylim([0 max(ylim)*1.05])
end 

set(gcf,'Position' , [682 67 605 731])
% subplot(5,4,18)
% hold on; 
% plot(bart_WT_batch_min_t, bart_WT_batch_min_rib.ras./bart_WT_batch_min_y(:,num_y.r0), 'Color', plt_clrs.green);
% plot(zamp_WT_batch_min_t, zamp_WT_batch_min_rib.ras./zamp_WT_batch_min_y(:,num_y.r0), 'Color', plt_clrs.blue);
% plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_rib.ras./murphy_WT_batch_YPD_y(:,num_y.r0), 'Color', plt_clrs.red);
% hold off;
% ylabel('R_{as}/R_{0}', 'Interpreter', 'none')   
% xlabel(x_label)
% axis square; 
% box on;
% ylim([0 max(ylim)*1.05])

% subplot(5,4,19)
% hold on; 
% plot(bart_WT_batch_min_t, bart_WT_batch_min_rib.rf./bart_WT_batch_min_y(:,num_y.r0), 'Color', plt_clrs.green);
% plot(zamp_WT_batch_min_t, zamp_WT_batch_min_rib.rf./zamp_WT_batch_min_y(:,num_y.r0), 'Color', plt_clrs.blue);
% plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_rib.rf./murphy_WT_batch_YPD_y(:,num_y.r0), 'Color', plt_clrs.red);
% hold off;
% ylabel('R_{f}/R_{0}', 'Interpreter', 'none')    
% xlabel(x_label)
% axis square; 
% box on;
% ylim([0 max(ylim)*1.05])
% 
% subplot(5,4,20)
% hold on; 
% plot(bart_WT_batch_min_t, bart_WT_batch_min_rib.rat./bart_WT_batch_min_y(:,num_y.r0), 'Color', plt_clrs.green);
% plot(zamp_WT_batch_min_t, zamp_WT_batch_min_rib.rat./zamp_WT_batch_min_y(:,num_y.r0), 'Color', plt_clrs.blue);
% plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_rib.rat./murphy_WT_batch_YPD_y(:,num_y.r0), 'Color', plt_clrs.red);
% hold off;
% ylabel('R_{at}/R_{0}', 'Interpreter', 'none')    
% xlabel(x_label)
% axis square; 
% box on;
% ylim([0 max(ylim)*1.05])

sgtitle('other growth intermediates')




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
plot(bart_WT_batch_min_t, bart_WT_batch_min_y(:,num_y.sc), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_y(:,num_y.sc), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y(:,num_y.sc), 'Color', plt_clrs.red);
ylabel(label.met{num_met.sc})    
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_y(:,num_y.sc)); 
ymax2 = max(zamp_WT_batch_min_y(:,num_y.sc)); 
ymax3 = max(murphy_WT_batch_YPD_y(:,num_y.sc)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on;
hold off; 

% J_sp
subplot(a, b, fig_index(2)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.flux(:,num_flux.sp), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.flux(:,num_flux.sp), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.flux(:,num_flux.sp), 'Color', plt_clrs.red);
ylabel(label.met_reac.flux{num_flux.sp})    
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.flux(:,num_flux.sp)); 
ymax2 = max(zamp_WT_batch_min_met_reac.flux(:,num_flux.sp)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.flux(:,num_flux.sp)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_sp
subplot(a, b, fig_index(3)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.prot(:,num_flux.sp), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.prot(:,num_flux.sp), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.sp), 'Color', plt_clrs.red);
ylabel(label.met_reac.prot{num_flux.sp})
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.prot(:,num_flux.sp)); 
ymax2 = max(zamp_WT_batch_min_met_reac.prot(:,num_flux.sp)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.prot(:,num_flux.sp)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.substrate(:,num_flux.sp), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.substrate(:,num_flux.sp), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.sp), 'Color', plt_clrs.red);
ylabel(label.met_reac.substrate{num_flux.sp}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.substrate(:,num_flux.sp)); 
ymax2 = max(zamp_WT_batch_min_met_reac.substrate(:,num_flux.sp)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.sp)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.sig(:,num_flux.sp), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.sig(:,num_flux.sp), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.sp), 'Color', plt_clrs.red);
ylabel(label.met_reac.sig{num_flux.sp}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.sig(:,num_flux.sp)); 
ymax2 = max(zamp_WT_batch_min_met_reac.sig(:,num_flux.sp)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.sig(:,num_flux.sp)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.atp(:,num_flux.sp), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.atp(:,num_flux.sp), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.sp), 'Color', plt_clrs.red);
ylabel(label.met_reac.atp{num_flux.sp}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.atp(:,num_flux.sp)); 
ymax2 = max(zamp_WT_batch_min_met_reac.atp(:,num_flux.sp)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.atp(:,num_flux.sp)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

fig_index = fig_index + b; 


% J_sd
subplot(a, b, fig_index(2)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.flux(:,num_flux.sd), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.flux(:,num_flux.sd), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.flux(:,num_flux.sd), 'Color', plt_clrs.red);
ylabel(label.met_reac.flux{num_flux.sd})    
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.flux(:,num_flux.sd)); 
ymax2 = max(zamp_WT_batch_min_met_reac.flux(:,num_flux.sd)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.flux(:,num_flux.sd)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_sp
subplot(a, b, fig_index(3)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.prot(:,num_flux.sd), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.prot(:,num_flux.sd), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.prot(:,num_flux.sd), 'Color', plt_clrs.red);
ylabel(label.met_reac.prot{num_flux.sd})
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.prot(:,num_flux.sd)); 
ymax2 = max(zamp_WT_batch_min_met_reac.prot(:,num_flux.sd)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.prot(:,num_flux.sd)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.substrate(:,num_flux.sd), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.substrate(:,num_flux.sd), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.sd), 'Color', plt_clrs.red);
ylabel(label.met_reac.substrate{num_flux.sd}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.substrate(:,num_flux.sd)); 
ymax2 = max(zamp_WT_batch_min_met_reac.substrate(:,num_flux.sd)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.substrate(:,num_flux.sd)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.sig(:,num_flux.sd), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.sig(:,num_flux.sd), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.sig(:,num_flux.sd), 'Color', plt_clrs.red);
ylabel(label.met_reac.sig{num_flux.sd}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.sig(:,num_flux.sd)); 
ymax2 = max(zamp_WT_batch_min_met_reac.sig(:,num_flux.sd)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.sig(:,num_flux.sd)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.atp(:,num_flux.sd), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.atp(:,num_flux.sd), 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_met_reac.atp(:,num_flux.sd), 'Color', plt_clrs.red);
ylabel(label.met_reac.atp{num_flux.sd}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.atp(:,num_flux.sd)); 
ymax2 = max(zamp_WT_batch_min_met_reac.atp(:,num_flux.sd)); 
ymax3 = max(murphy_WT_batch_YPD_met_reac.atp(:,num_flux.sd)); 
if max([ymax1, ymax2, ymax3]) > 0 
ylim([0 1.05*max([ymax1, ymax2, ymax3])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 


% --------------------------------------------------------------------------



bart_total_protein_con   = sum(bart_WT_batch_min_y(:,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);
zamp_total_protein_con   = sum(zamp_WT_batch_min_y(:,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);
murphy_total_protein_con = sum(murphy_WT_batch_YPD_y(:,num_y.r:num_y.lo).* par.l(num_prot.r:num_prot.lo)',2);


% proteins frac (sim)
bart_WT_batch_min_y(:,[num_y.r: num_y.sd]) = (bart_WT_batch_min_y(:,[num_y.r: num_y.sd]).*par.l')./ bart_total_protein_con;
zamp_WT_batch_min_y(:,[num_y.r: num_y.sd]) = (zamp_WT_batch_min_y(:,[num_y.r: num_y.sd]).*par.l')./ zamp_total_protein_con;
murphy_WT_batch_YPD_y(:,[num_y.r: num_y.sd]) = (murphy_WT_batch_YPD_y(:,[num_y.r: num_y.sd]).*par.l')./ murphy_total_protein_con;


%% protein fraction 
figure; 
a = 3; b = 5; 
prot_ = fieldnames(num_prot); 
for i = 1: max(struct2array(num_prot))
subplot(a, b, i)
hold on;
%plot(hackett.dr, hackett.r, 'o', 'MarkerFaceColor', color_sc.hackett)
plot(bart_WT_batch_min_t, bart_WT_batch_min_y(:,num_y.(prot_{i}))  , '-', 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_y(:,num_y.(prot_{i})) , '-', 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y(:,num_y.(prot_{i})) , '-', 'Color', plt_clrs.red);
ylabel(label.prot{num_prot.(prot_{i})})
xlabel(x_label);
if max(ylim) > 0 
ylim([0 1.1.*max(ylim)])
end 
%xlim([0 50])
%xline(d_crit)
axis square; 
box on;
hold off;
end 
sgtitle('protein fraction')

set(gcf, 'Position', [1441 161 698 637])
%% protein concentration

figure; 
a = 4; b = 4; 
prot_ = fieldnames(num_prot); 
for i = 1: max(struct2array(num_prot))
subplot(a, b, i)
hold on;
%plot(hackett.dr, hackett.r, 'o', 'MarkerFaceColor', color_sc.hackett)
plot(bart_WT_batch_min_t, bart_WT_batch_min_y(:,num_y.(prot_{i})).*bart_total_protein_con  , '-', 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_y(:,num_y.(prot_{i})).*zamp_total_protein_con , '-', 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_y(:,num_y.(prot_{i})).*murphy_total_protein_con  , '-', 'Color', plt_clrs.red);
ylabel(label.prot{num_prot.(prot_{i})})
xlabel(x_label);
if max(ylim) > 0 
ylim([0 1.1.*max(ylim)])
end 
%xlim([0 50])
%xline(d_crit)
axis square; 
box on;
hold off;
end 

sgtitle('protein concentration')
set(gcf, 'Position', [1441 161 698 637])

%% protein alpha


figure; 
a = 4; b = 4; 
prot_ = fieldnames(num_prot); 
for i = 1: max(struct2array(num_prot))
subplot(a, b, i)
hold on;
%plot(hackett.dr, hackett.r, 'o', 'MarkerFaceColor', color_sc.hackett)
plot(bart_WT_batch_min_t, bart_WT_batch_min_prot_syn.alpha(:,num_prot.(prot_{i}))  , '-', 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_prot_syn.alpha(:,num_prot.(prot_{i})) , '-', 'Color', plt_clrs.blue);
plot(murphy_WT_batch_YPD_t, murphy_WT_batch_YPD_prot_syn.alpha(:,num_prot.(prot_{i}))  , '-', 'Color', plt_clrs.red);
ylabel(label.prot{num_prot.(prot_{i})})
xlabel(x_label);
if max(ylim) > 0 
ylim([0 1.1.*max(ylim)])
end 
%xlim([0 50])
%xline(d_crit)
axis square; 
box on;
hold off;
end 

sgtitle('protein concentration')
set(gcf, 'Position', [1441 161 698 637])

end 
