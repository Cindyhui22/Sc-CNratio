function plot_all_batch_intermediates_tingaksed(bart_WT_batch_min_t, zamp_WT_batch_min_t, murphy_WT_batch_YPD_t, ...
bart_WT_batch_min_y, zamp_WT_batch_min_y, murphy_WT_batch_YPD_y, ...
bart_WT_batch_min_met_reac, zamp_WT_batch_min_met_reac, murphy_WT_batch_YPD_met_reac, ...
bart_WT_batch_min_other_met_reac, zamp_WT_batch_min_other_met_reac, murphy_WT_batch_YPD_other_met_reac, ...
bart_WT_batch_min_prot_syn, zamp_WT_batch_min_prot_syn, murphy_WT_batch_YPD_prot_syn, ...
bart_WT_batch_min_rib, zamp_WT_batch_min_rib, murphy_WT_batch_YPD_rib, ...
label, plt_colors, par, num_y, num_flux, x_label)



% tRNA contributions
%tc_contri = bart_WT_batch_min_prot_syn.tc./par.K_bart_WT_batch_min_rib_tc; 
%tu_contri = bart_WT_batch_min_prot_syn.tu./par.K_bart_WT_batch_min_rib_tu;
%rat_tc_contri = tc_contri./(1 + tc_contri + tu_contri);
%ras_tu_contri = tu_contri./(1 + tc_contri + tu_contri);
set(0, 'DefaultLineLineWidth',  1.5);


paper_figure_color;
x_lim = [0 20]; 
y_lim_flux = [0 *10^4 4 *10^6]; 

%% virtual metabolites - precursor 
act_no_sub = 1; 
figure; 
a = 5+3;
b = 7;
fig_index = 0:6;
fig_index(2) = 1; 
fig_index(3) = 2; 
fig_index(1) = 4; 
fig_index(4) = 3; 
fig_index(5) = 5; 
fig_index(6) = 7; 
fig_index(7) = 6;
% precursor 


% J_gy
subplot(a, b, fig_index(2)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.flux(:,num_flux.gy), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.flux(:,num_flux.gy), 'Color', color_zamp);
ylabel(label.met_reac.flux{num_flux.gy})    
xlabel(x_label)
ylim(y_lim_flux)
xlim(x_lim)
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_gy
subplot(a, b, fig_index(3)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.prot(:,num_flux.gy), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.prot(:,num_flux.gy), 'Color', color_zamp);
ylabel(label.met_reac.prot{num_flux.gy})
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.prot(:,num_flux.gy)); 
ymax2 = max(zamp_WT_batch_min_met_reac.prot(:,num_flux.gy)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.1*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

activity = bart_WT_batch_min_met_reac.substrate(:,num_flux.gy).*bart_WT_batch_min_met_reac.snf1(:,num_flux.gy).* bart_WT_batch_min_met_reac.atp(:,num_flux.gy);
if act_no_sub == 1
activity = 1.*bart_WT_batch_min_met_reac.sig(:,num_flux.gy).* bart_WT_batch_min_met_reac.atp(:,num_flux.gy);
activity2 = 1.*zamp_WT_batch_min_met_reac.sig(:,num_flux.gy).* zamp_WT_batch_min_met_reac.atp(:,num_flux.gy);
end 
subplot(a, b, fig_index(1)); 
hold on; 
plot(bart_WT_batch_min_t, activity, 'Color', color_bart);
plot(zamp_WT_batch_min_t, activity2, 'Color', color_zamp);
ylabel('activity'); ylim([0 1.1*max(ylim)])
xlabel(x_label)
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 

% k*Gl/(km+Gl) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.substrate(:,num_flux.gy), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.substrate(:,num_flux.gy), 'Color', color_zamp);
ylabel(label.met_reac.substrate{num_flux.gy}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.substrate(:,num_flux.gy)); 
ymax2 = max(zamp_WT_batch_min_met_reac.substrate(:,num_flux.gy)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.1*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.snf1(:,num_flux.gy), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.snf1(:,num_flux.gy), 'Color', color_zamp);
ylabel(label.met_reac.snf1{num_flux.gy}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.snf1(:,num_flux.gy)); 
ymax2 = max(zamp_WT_batch_min_met_reac.snf1(:,num_flux.gy)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.1*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 


subplot(a, b, fig_index(7)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.pc(:,num_flux.gy), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.pc(:,num_flux.gy), 'Color', color_zamp);
ylabel(label.met_reac.pc{num_flux.gy}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.pc(:,num_flux.gy)); 
ymax2 = max(zamp_WT_batch_min_met_reac.pc(:,num_flux.gy)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.1*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 


% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.atp(:,num_flux.gy), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.atp(:,num_flux.gy), 'Color', color_zamp);
ylabel(label.met_reac.atp{num_flux.gy}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.atp(:,num_flux.gy)); 
ymax2 = max(zamp_WT_batch_min_met_reac.atp(:,num_flux.gy)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.1*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 


fig_index = fig_index + b; 
%--------------------------------------------------------------------------

% J_as
subplot(a, b, fig_index(2)); 
hold on; 
plot(bart_WT_batch_min_t, 2*bart_WT_batch_min_met_reac.flux(:,num_flux.as), 'Color', color_bart);
plot(zamp_WT_batch_min_t, 2*zamp_WT_batch_min_met_reac.flux(:,num_flux.as), 'Color', color_zamp);
ylabel(label.met_reac.flux{num_flux.as})    
xlabel(x_label) 
ylim(y_lim_flux)
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_as
subplot(a, b, fig_index(3)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.prot(:,num_flux.as), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.prot(:,num_flux.as), 'Color', color_zamp);
ylabel(label.met_reac.prot{num_flux.as})
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.prot(:,num_flux.as)); 
ymax2 = max(zamp_WT_batch_min_met_reac.prot(:,num_flux.as)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.1*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

activity = bart_WT_batch_min_met_reac.substrate(:,num_flux.as).*bart_WT_batch_min_met_reac.snf1(:,num_flux.as).* bart_WT_batch_min_met_reac.atp(:,num_flux.as);
if act_no_sub == 1
activity = 1.*bart_WT_batch_min_met_reac.sig(:,num_flux.as).* bart_WT_batch_min_met_reac.atp(:,num_flux.as);
activity2 = 1.*zamp_WT_batch_min_met_reac.sig(:,num_flux.as).* zamp_WT_batch_min_met_reac.atp(:,num_flux.as);
end 
subplot(a, b, fig_index(1)); 
hold on; 
plot(bart_WT_batch_min_t, activity, 'Color', color_bart);
plot(zamp_WT_batch_min_t, activity2, 'Color', color_zamp);
ylabel('activity'); ylim([0 1.1*max(ylim)])
xlabel(x_label)
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.substrate(:,num_flux.as), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.substrate(:,num_flux.as), 'Color', color_zamp);
ylabel(label.met_reac.substrate{num_flux.as}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.substrate(:,num_flux.as)); 
ymax2 = max(zamp_WT_batch_min_met_reac.substrate(:,num_flux.as)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.1*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(7)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.pc(:,num_flux.as), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.pc(:,num_flux.as), 'Color', color_zamp);
ylabel(label.met_reac.pc{num_flux.as}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.snf1(:,num_flux.as)); 
ymax2 = max(zamp_WT_batch_min_met_reac.snf1(:,num_flux.as)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.1*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on;
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.atp(:,num_flux.as), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.atp(:,num_flux.as), 'Color', color_zamp);
ylabel(label.met_reac.atp{num_flux.as}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.atp(:,num_flux.as)); 
ymax2 = max(zamp_WT_batch_min_met_reac.atp(:,num_flux.as)); 
% if max([ymax1, ymax2]) > 0 
%     ylim([0 1.1*max([ymax1, ymax2])])
% end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

fig_index = fig_index + b; 

% J_mt
subplot(a, b, fig_index(2)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.flux(:,num_flux.mt), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.flux(:,num_flux.mt), 'Color', color_zamp);
ylabel(label.met_reac.flux{num_flux.mt})    
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.flux(:,num_flux.mt)); 
ymax2 = max(zamp_WT_batch_min_met_reac.flux(:,num_flux.mt)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.1*max([ymax1, ymax2])])
end
ylim(y_lim_flux)
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_mt
subplot(a, b, fig_index(3)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.prot(:,num_flux.mt), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.prot(:,num_flux.mt), 'Color', color_zamp);
ylabel(label.met_reac.prot{num_flux.mt})
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.prot(:,num_flux.mt)); 
ymax2 = max(zamp_WT_batch_min_met_reac.prot(:,num_flux.mt)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.1*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 


activity = bart_WT_batch_min_met_reac.substrate(:,num_flux.mt).*bart_WT_batch_min_met_reac.snf1(:,num_flux.mt).* bart_WT_batch_min_met_reac.atp(:,num_flux.mt);
if act_no_sub == 1
activity = 1*bart_WT_batch_min_met_reac.sig(:,num_flux.mt).* bart_WT_batch_min_met_reac.atp(:,num_flux.mt);
activity2 = 1.*zamp_WT_batch_min_met_reac.sig(:,num_flux.mt).* zamp_WT_batch_min_met_reac.atp(:,num_flux.mt);
end
subplot(a, b, fig_index(1)); 
hold on; 
plot(bart_WT_batch_min_t, activity, 'Color', color_bart);
plot(zamp_WT_batch_min_t, activity2, 'Color', color_zamp);
ylabel('activity'); ylim([0 1.1*max(ylim)])
xlabel(x_label)
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 

% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.substrate(:,num_flux.mt), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.substrate(:,num_flux.mt), 'Color', color_zamp);
ylabel(label.met_reac.substrate{num_flux.mt}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.substrate(:,num_flux.mt)); 
ymax2 = max(zamp_WT_batch_min_met_reac.substrate(:,num_flux.mt)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.1*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(7)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.pc(:,num_flux.mt), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.pc(:,num_flux.mt), 'Color', color_zamp);
ylabel(label.met_reac.pc{num_flux.mt}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.snf1(:,num_flux.mt)); 
ymax2 = max(zamp_WT_batch_min_met_reac.snf1(:,num_flux.mt)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.1*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.atp(:,num_flux.mt), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.atp(:,num_flux.mt), 'Color', color_zamp);
ylabel(label.met_reac.atp{num_flux.mt}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.atp(:,num_flux.mt)); 
ymax2 = max(zamp_WT_batch_min_met_reac.atp(:,num_flux.mt)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.1*max([ymax1, ymax2])])
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
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.flux(:,num_flux.fe), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.flux(:,num_flux.fe), 'Color', color_zamp);
ylabel(label.met_reac.flux{num_flux.fe})    
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.flux(:,num_flux.fe)); 
ymax2 = max(zamp_WT_batch_min_met_reac.flux(:,num_flux.fe)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.1*max([ymax1, ymax2])])
end
ylim(y_lim_flux)
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_fe
subplot(a, b, fig_index(3)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.prot(:,num_flux.fe), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.prot(:,num_flux.fe), 'Color', color_zamp);
ylabel(label.met_reac.prot{num_flux.fe})
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.prot(:,num_flux.fe)); 
ymax2 = max(zamp_WT_batch_min_met_reac.prot(:,num_flux.fe)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.1*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

activity = bart_WT_batch_min_met_reac.substrate(:,num_flux.fe).*bart_WT_batch_min_met_reac.snf1(:,num_flux.fe).* bart_WT_batch_min_met_reac.atp(:,num_flux.fe);
if act_no_sub == 1
activity = 1.*bart_WT_batch_min_met_reac.sig(:,num_flux.fe).* bart_WT_batch_min_met_reac.atp(:,num_flux.fe);
activity2 = 1.*zamp_WT_batch_min_met_reac.sig(:,num_flux.fe).* zamp_WT_batch_min_met_reac.atp(:,num_flux.fe);

end 
subplot(a, b, fig_index(1)); 
hold on; 
plot(bart_WT_batch_min_t, activity, 'Color', color_bart);
plot(zamp_WT_batch_min_t, activity2, 'Color', color_zamp);
ylabel('activity'); ylim([0 1.1*max(ylim)])
xlabel(x_label)
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 


% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.substrate(:,num_flux.fe), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.substrate(:,num_flux.fe), 'Color', color_zamp);
ylabel(label.met_reac.substrate{num_flux.fe}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.substrate(:,num_flux.fe)); 
ymax2 = max(zamp_WT_batch_min_met_reac.substrate(:,num_flux.fe)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.1*max([ymax1, ymax2])])
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
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.flux(:,num_flux.gn), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.flux(:,num_flux.gn), 'Color', color_zamp);
ylabel(label.met_reac.flux{num_flux.gn})    
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.flux(:,num_flux.gn)); 
ymax2 = max(zamp_WT_batch_min_met_reac.flux(:,num_flux.gn)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.1*max([ymax1, ymax2])])
end
ylim(y_lim_flux)
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_gn
subplot(a, b, fig_index(3)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.prot(:,num_flux.gn), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.prot(:,num_flux.gn), 'Color', color_zamp);
ylabel(label.met_reac.prot{num_flux.gn})
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.prot(:,num_flux.gn)); 
ymax2 = max(zamp_WT_batch_min_met_reac.prot(:,num_flux.gn)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.1*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

activity = bart_WT_batch_min_met_reac.substrate(:,num_flux.gn).*bart_WT_batch_min_met_reac.snf1(:,num_flux.gn).* bart_WT_batch_min_met_reac.atp(:,num_flux.gn);
if act_no_sub == 1
activity = 1.*bart_WT_batch_min_met_reac.sig(:,num_flux.gn).* bart_WT_batch_min_met_reac.atp(:,num_flux.gn);
activity2 = 1.*zamp_WT_batch_min_met_reac.sig(:,num_flux.gn).* zamp_WT_batch_min_met_reac.atp(:,num_flux.gn);

end 
subplot(a, b, fig_index(1)); 
hold on; 
plot(bart_WT_batch_min_t, activity, 'Color', color_bart);
plot(zamp_WT_batch_min_t, activity2, 'Color', color_zamp);
ylabel('activity'); ylim([0 1.1*max(ylim)])
xlabel(x_label)
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 


% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.substrate(:,num_flux.gn), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.substrate(:,num_flux.gn), 'Color', color_zamp);
ylabel(label.met_reac.substrate{num_flux.gn}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.substrate(:,num_flux.gn)); 
ymax2 = max(zamp_WT_batch_min_met_reac.substrate(:,num_flux.gn)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.1*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 


% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.snf1(:,num_flux.gn), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.snf1(:,num_flux.gn), 'Color', color_zamp);
ylabel(label.met_reac.snf1{num_flux.gn}) 
xlabel(x_label)
ymax1 = max(bart_WT_batch_min_met_reac.snf1(:,num_flux.gn)); 
ymax2 = max(zamp_WT_batch_min_met_reac.snf1(:,num_flux.gn)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.1*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold on;

% ---------------------------------------------------------------------------------
fig_index = fig_index + b; 
% J_fe_cit
subplot(a, b, fig_index(2)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.flux(:,num_flux.lp_cit), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.flux(:,num_flux.lp_cit), 'Color', color_zamp);
ylabel(label.met_reac.flux{num_flux.lp_cit})    
xlabel(x_label)
ylim(y_lim_flux)
xlim(x_lim)
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_lp
subplot(a, b, fig_index(3)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.prot(:,num_flux.lp_cit), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.prot(:,num_flux.lp_cit), 'Color', color_zamp);
ylabel(label.met_reac.prot{num_flux.lp_cit})
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.1*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

activity = bart_WT_batch_min_met_reac.substrate(:,num_flux.lp_cit).*bart_WT_batch_min_met_reac.snf1(:,num_flux.lp_cit).* bart_WT_batch_min_met_reac.atp(:,num_flux.lp_cit);
if act_no_sub == 1
activity = 1.*bart_WT_batch_min_met_reac.sig(:,num_flux.lp_cit).* bart_WT_batch_min_met_reac.atp(:,num_flux.lp_cit);
activity2 = 1.*zamp_WT_batch_min_met_reac.sig(:,num_flux.lp_cit).* zamp_WT_batch_min_met_reac.atp(:,num_flux.lp_cit);
end 
subplot(a, b, fig_index(1)); 
hold on; 
plot(bart_WT_batch_min_t, activity, 'Color', color_bart);
plot(zamp_WT_batch_min_t, activity2, 'Color', color_zamp);
ylabel('activity'); 
if max(ylim) > 0 
ylim([0 1.1*max(ylim)])
end
xlabel(x_label)
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 

% k*Gl/(km+Gl) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.substrate(:,num_flux.lp_cit), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.substrate(:,num_flux.lp_cit), 'Color', color_zamp);
ylabel(label.met_reac.substrate{num_flux.lp_cit}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.1*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.snf1(:,num_flux.lp_cit), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.snf1(:,num_flux.lp_cit), 'Color', color_zamp);
ylabel(label.met_reac.snf1{num_flux.lp_cit}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.1*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 


subplot(a, b, fig_index(7)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.pc(:,num_flux.lp_cit), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.pc(:,num_flux.lp_cit), 'Color', color_zamp);
ylabel(label.met_reac.pc{num_flux.lp_cit}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.1*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 


% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.atp(:,num_flux.lp_cit), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.atp(:,num_flux.lp_cit), 'Color', color_zamp);
ylabel(label.met_reac.atp{num_flux.lp_cit}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.1*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 


% ---------------------------------------------------------------------------------
fig_index = fig_index + b; 
% J_lp_fe
subplot(a, b, fig_index(2)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.flux(:,num_flux.lp_fe), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.flux(:,num_flux.lp_fe), 'Color', color_zamp);
ylabel(label.met_reac.flux{num_flux.lp_fe})    
xlabel(x_label)
ylim(y_lim_flux)
xlim(x_lim)
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_lp
subplot(a, b, fig_index(3)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.prot(:,num_flux.lp_fe), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.prot(:,num_flux.lp_fe), 'Color', color_zamp);
ylabel(label.met_reac.prot{num_flux.lp_fe})
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.1*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

activity = bart_WT_batch_min_met_reac.substrate(:,num_flux.lp_fe).*bart_WT_batch_min_met_reac.snf1(:,num_flux.lp_fe).* bart_WT_batch_min_met_reac.atp(:,num_flux.lp_fe);
if act_no_sub == 1
activity = 1.*bart_WT_batch_min_met_reac.sig(:,num_flux.lp_fe).* bart_WT_batch_min_met_reac.atp(:,num_flux.lp_fe);
activity2 = 1.*zamp_WT_batch_min_met_reac.sig(:,num_flux.lp_fe).* zamp_WT_batch_min_met_reac.atp(:,num_flux.lp_fe);
end 
subplot(a, b, fig_index(1)); 
hold on; 
plot(bart_WT_batch_min_t, activity, 'Color', color_bart);
plot(zamp_WT_batch_min_t, activity2, 'Color', color_zamp);
ylabel('activity'); 
if max(ylim) > 0 
ylim([0 1.1*max(ylim)])
end
xlabel(x_label)
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 

% k*Gl/(km+Gl) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.substrate(:,num_flux.lp_fe), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.substrate(:,num_flux.lp_fe), 'Color', color_zamp);
ylabel(label.met_reac.substrate{num_flux.lp_fe}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.1*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.snf1(:,num_flux.lp_fe), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.snf1(:,num_flux.lp_fe), 'Color', color_zamp);
ylabel(label.met_reac.snf1{num_flux.lp_fe}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.1*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 


subplot(a, b, fig_index(7)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.pc(:,num_flux.lp_fe), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.pc(:,num_flux.lp_fe), 'Color', color_zamp);
ylabel(label.met_reac.pc{num_flux.lp_fe}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.1*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 


% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.atp(:,num_flux.lp_fe), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.atp(:,num_flux.lp_fe), 'Color', color_zamp);
ylabel(label.met_reac.atp{num_flux.lp_fe}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.1*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 


% ---------------------------------------------------------------------------------
fig_index = fig_index + b; 
% J_lo
subplot(a, b, fig_index(2)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.flux(:,num_flux.lo), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.flux(:,num_flux.lo), 'Color', color_zamp);
ylabel(label.met_reac.flux{num_flux.lo})    
xlabel(x_label)
ylim(y_lim_flux)
xlim(x_lim)
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% E_lo
subplot(a, b, fig_index(3)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.prot(:,num_flux.lo), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.prot(:,num_flux.lo), 'Color', color_zamp);
ylabel(label.met_reac.prot{num_flux.lo})
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.1*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

activity = bart_WT_batch_min_met_reac.substrate(:,num_flux.lo).*bart_WT_batch_min_met_reac.snf1(:,num_flux.lo).* bart_WT_batch_min_met_reac.atp(:,num_flux.lo);
if act_no_sub == 1
activity = 1.*bart_WT_batch_min_met_reac.sig(:,num_flux.lo).* bart_WT_batch_min_met_reac.atp(:,num_flux.lo);
activity2 = 1.*zamp_WT_batch_min_met_reac.sig(:,num_flux.lo).* zamp_WT_batch_min_met_reac.atp(:,num_flux.lo);
end 
subplot(a, b, fig_index(1)); 
hold on; 
plot(bart_WT_batch_min_t, activity, 'Color', color_bart);
plot(zamp_WT_batch_min_t, activity2, 'Color', color_zamp);
ylabel('activity'); 
if max(ylim) > 0 
ylim([0 1.1*max(ylim)])
end
xlabel(x_label)
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 

% k*Gl/(km+Gl) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.substrate(:,num_flux.lo), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.substrate(:,num_flux.lo), 'Color', color_zamp);
ylabel(label.met_reac.substrate{num_flux.lo}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.1*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.snf1(:,num_flux.lo), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.snf1(:,num_flux.lo), 'Color', color_zamp);
ylabel(label.met_reac.snf1{num_flux.lo}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.1*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 


subplot(a, b, fig_index(7)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.pc(:,num_flux.lo), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.pc(:,num_flux.lo), 'Color', color_zamp);
ylabel(label.met_reac.pc{num_flux.lo}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.1*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 


% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.atp(:,num_flux.lo), 'Color', color_bart);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.atp(:,num_flux.lo), 'Color', color_zamp);
ylabel(label.met_reac.atp{num_flux.lo}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.1*max(ylim)])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 


sgtitle('precursor')



%mmki(aa_in,  par.I_lp_aa,  par.n_lp_aa)
lp_aa_term_1gL = (par.I_lp_aa^par.n_lp_aa)./(par.I_lp_aa^par.n_lp_aa + bart_WT_batch_min_y(:,num_y.aa_in).^par.n_lp_aa);
lp_aa_term_10gL = (par.I_lp_aa^par.n_lp_aa)./(par.I_lp_aa^par.n_lp_aa + zamp_WT_batch_min_y(:,num_y.aa_in).^par.n_lp_aa);

%hill(nh4, par.K_as_nh, 1)
as_nh4_term_1gL = (bart_WT_batch_min_y(:,num_y.nh4).^1)./(par.K_as_nh^1+ bart_WT_batch_min_y(:,num_y.nh4).^1);
as_nh4_term_10gL = (zamp_WT_batch_min_y(:,num_y.nh4).^1)./(par.K_as_nh^1+ zamp_WT_batch_min_y(:,num_y.nh4).^1); 

%hill(pc,    par.K_as_p,        par.n_as_pc)
as_sub_term_1gL  = (bart_WT_batch_min_y(:,num_y.pc).^par.n_as_pc)./(par.K_as_p^par.n_as_pc+ bart_WT_batch_min_y(:,num_y.pc).^par.n_as_pc);
as_sub_term_10gL = (zamp_WT_batch_min_y(:,num_y.pc).^par.n_as_pc)./(par.K_as_p^par.n_as_pc+ zamp_WT_batch_min_y(:,num_y.pc).^par.n_as_pc);
%% addtional figures for newly added variables
figure; 
a = 4 ;b = 5; 

% J_as
subplot(a, b, 1); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.flux(:,num_flux.as), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.flux(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.as})    
xlabel(x_label) 
if max(ylim) > 0
ylim([0 max(ylim)])
end 
xlim(x_lim)
%xticks(x_ticks)
axis square; 
box on; 
hold off; 

% E_as
subplot(a, b, 2); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.prot(:,num_flux.as), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.prot(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.as})
xlabel(x_label)
if max(ylim) > 0
ylim([0 max(ylim)])
end 
xlim(x_lim)
%xticks(x_ticks)
axis square; 
box on; 
hold off; 




% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, 3); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.substrate(:,num_flux.as), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.substrate(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.as}) 
xlabel(x_label)
if max(ylim) > 0
ylim([0 max(ylim)])
end 
xlim(x_lim)
%xticks(x_ticks)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, 4); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.pc(:,num_flux.as), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.pc(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.pc{num_flux.as}) 
xlabel(x_label)
if max(ylim) > 0
ylim([0 max(ylim)])
end 
xlim(x_lim)
%xticks(x_ticks)
axis square; 
box on;
hold off; 

% ATP contribution 
subplot(a, b, 5); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.atp(:,num_flux.as), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.atp(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.atp{num_flux.as}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
%xticks(x_ticks)
axis square; 
box on; 
hold off; 



% Pc part in Jas substrate 
subplot(a, b, 3+b); 
hold on; 
plot(bart_WT_batch_min_t, as_sub_term_1gL, 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, as_sub_term_10gL, 'Color', plt_clrs.blue);
ylabel({'[P_{c}]/' , '(K+[P_{c}])'}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
%xticks(x_ticks)
axis square; 
box on;
hold off; 


% NH4 part in Jas substrate 
subplot(a, b, 4+b); 
hold on; 
plot(bart_WT_batch_min_t, as_nh4_term_1gL, 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, as_nh4_term_10gL, 'Color', plt_clrs.blue);
ylabel({'[NH_{4}]/' , '(K+[NH_{4}])'}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
%xticks(x_ticks)
axis square; 
box on;
hold off; 



% ------------------------------------------------------------------------
% Jlp,fe 

subplot(a, b, 1+2*b); 
hold on; 
plot(bart_WT_batch_min_t, 2*bart_WT_batch_min_met_reac.flux(:,num_flux.lp_fe), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, 2*zamp_WT_batch_min_met_reac.flux(:,num_flux.lp_fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.lp_fe})    
xlabel(x_label) 
if max(ylim) > 0
ylim([0 max(ylim)])
end 
xlim(x_lim)
%xticks(x_ticks)
axis square; 
box on; 
hold off; 

% E_lp
subplot(a, b, 2+2*b); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.prot(:,num_flux.lp_fe), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.prot(:,num_flux.lp_fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.lp_fe})
xlabel(x_label)
if max(ylim) > 0
ylim([0 max(ylim)])
end 
xlim(x_lim)
%xticks(x_ticks)
axis square; 
box on; 
hold off; 


% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, 3+2*b); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.substrate(:,num_flux.lp_fe), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.substrate(:,num_flux.lp_fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.lp_fe}) 
xlabel(x_label)
if max(ylim) > 0
ylim([0 max(ylim)])
end 
xlim(x_lim)
%xticks(x_ticks)
axis square; 
box on; 
hold off; 

% SNF1 contribution
subplot(a, b, 4+2*b); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.snf1(:,num_flux.lp_fe), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.snf1(:,num_flux.lp_fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.snf1{num_flux.lp_fe}) 
xlabel(x_label)
if max(ylim) > 0
ylim([0 max(ylim)])
end 
xlim(x_lim)
%xticks(x_ticks)
axis square; 
box on;
hold off; 



% ------------------------------------------------------------------------
% Jlp,cit 

subplot(a, b, 1+3*b); 
hold on; 
plot(bart_WT_batch_min_t, 2*bart_WT_batch_min_met_reac.flux(:,num_flux.lp_cit), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, 2*zamp_WT_batch_min_met_reac.flux(:,num_flux.lp_cit), 'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.lp_cit})    
xlabel(x_label) 
if max(ylim) > 0
ylim([0 max(ylim)])
end 
xlim(x_lim)
%xticks(x_ticks)
axis square; 
box on; 
hold off; 

% E_lp
subplot(a, b, 2+3*b); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.prot(:,num_flux.lp_cit), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.prot(:,num_flux.lp_cit), 'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.lp_cit})
xlabel(x_label)
if max(ylim) > 0
ylim([0 max(ylim)])
end 
xlim(x_lim)
%xticks(x_ticks)
axis square; 
box on; 
hold off; 


% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, 3+3*b); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.substrate(:,num_flux.lp_cit), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.substrate(:,num_flux.lp_cit), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.lp_cit}) 
xlabel(x_label)
if max(ylim) > 0
ylim([0 max(ylim)])
end 
xlim(x_lim)
%xticks(x_ticks)
axis square; 
box on; 
hold off; 

% SNF1 contribution
subplot(a, b, 4+3*b); 
hold on; 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.snf1(:,num_flux.lp_cit), 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, zamp_WT_batch_min_met_reac.snf1(:,num_flux.lp_cit), 'Color', plt_clrs.blue);
ylabel(label.met_reac.snf1{num_flux.lp_fe}) 
xlabel(x_label)
if max(ylim) > 0
ylim([0 max(ylim)])
end 
xlim(x_lim)
%xticks(x_ticks)
axis square; 
box on;
hold off; 


% SNF1 contribution
subplot(a, b, 5+3*b); 
hold on; 
plot(bart_WT_batch_min_t, lp_aa_term_1gL, 'Color', plt_clrs.green);
plot(zamp_WT_batch_min_t, lp_aa_term_10gL, 'Color', plt_clrs.blue);
ylabel('N term') 
xlabel(x_label)
if max(ylim) > 0
ylim([0 max(ylim)])
end 
xlim(x_lim)
%xticks(x_ticks)
axis square; 
box on;
hold off; 

end 
