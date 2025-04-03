function plot_all_chemostat_intermediates_tingaksed(D, ...
WT_y_steady_1gL, WT_y_steady_10gL, ...
WT_met_reac_steady_1gL, WT_met_reac_steady_10gL, ...
WT_other_met_reac_steady_1gL, WT_other_met_reac_steady_10gL, ...
WT_prot_syn_steady_1gL, WT_prot_syn_steady_10gL, ...
WT_rib_steady_1gL, WT_rib_steady_10gL, ...
label, plt_clrs, par, par2, num_y, num_flux, x_label)

% tRNA contributions
%tc_contri = bart_WT_batch_min_prot_syn.tc./par.K_bart_WT_batch_min_rib_tc; 
%tu_contri = bart_WT_batch_min_prot_syn.tu./par.K_bart_WT_batch_min_rib_tu;
%rat_tc_contri = tc_contri./(1 + tc_contri + tu_contri);
%ras_tu_contri = tu_contri./(1 + tc_contri + tu_contri);
set(0, 'DefaultLineLineWidth',  1.5);
y_lim_flux = [0 3*10^6]; 
x_lim      = [0   0.42]; 
x_ticks    = [0 0.2 0.4]; 

% %mmki(aa_in,  par.I_lp_aa,  par.n_lp_aa)
% lp_aa_term_1gL = (par.I_lp_aa^par.n_lp_aa)./(par.I_lp_aa^par.n_lp_aa + WT_y_steady_1gL(:,num_y.aa_in).^par.n_lp_aa);
% lp_aa_term_10gL = (par.I_lp_aa^par.n_lp_aa)./(par.I_lp_aa^par.n_lp_aa + WT_y_steady_10gL(:,num_y.aa_in).^par.n_lp_aa);
% 
% %hill(nh4, par.K_as_nh, 1)
% as_nh4_term_1gL = (WT_y_steady_1gL(:,num_y.nh4).^1)./(par.K_as_nh^1+ WT_y_steady_1gL(:,num_y.nh4).^1);
% as_nh4_term_10gL = (WT_y_steady_10gL(:,num_y.nh4).^1)./(par.K_as_nh^1+ WT_y_steady_10gL(:,num_y.nh4).^1); 
% 
% %hill(pc,    par.K_as_p,        par.n_as_pc)
% as_sub_term_1gL  = (WT_y_steady_1gL(:,num_y.pc).^par.n_as_pc)./(par.K_as_p^par.n_as_pc+ WT_y_steady_1gL(:,num_y.pc).^par.n_as_pc);
% as_sub_term_10gL = (WT_y_steady_10gL(:,num_y.pc).^par.n_as_pc)./(par.K_as_p^par.n_as_pc+ WT_y_steady_10gL(:,num_y.pc).^par.n_as_pc);


%mmki(aa_in,  par.I_lp_aa,  par.n_lp_aa)
lp_aa_term_1gL = (par.I_lp_aa^par.n_lp_aa)./(par.I_lp_aa^par.n_lp_aa + WT_y_steady_1gL(:,num_y.aa_in).^par.n_lp_aa);
lp_aa_term_10gL = (par2.I_lp_aa^par2.n_lp_aa)./(par2.I_lp_aa^par2.n_lp_aa + WT_y_steady_10gL(:,num_y.aa_in).^par2.n_lp_aa);

%hill(nh4, par.K_as_nh, 1)
as_nh4_term_1gL = (WT_y_steady_1gL(:,num_y.nh4).^1)./(par.K_as_nh^1+ WT_y_steady_1gL(:,num_y.nh4).^1);
as_nh4_term_10gL = (WT_y_steady_10gL(:,num_y.nh4).^1)./(par2.K_as_nh4^1+ WT_y_steady_10gL(:,num_y.nh4).^1); 

%hill(pc,    par.K_as_p,        par.n_as_pc)
as_sub_term_1gL  = (WT_y_steady_1gL(:,num_y.pc).^par.n_as_pc)./(par.K_as_p^par.n_as_pc+ WT_y_steady_1gL(:,num_y.pc).^par.n_as_pc);
as_sub_term_10gL = (WT_y_steady_10gL(:,num_y.pc).^par2.n_as_p)./(par2.K_as_p^par2.n_as_p+ WT_y_steady_10gL(:,num_y.pc).^par2.n_as_p);

%{
%% real metabolites 
figure; 
a = 3;
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
ylabel(label.met{2})    
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
plot(D, WT_met_reac_steady_1gL.snf1(:,num_flux.gy), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.snf1(:,num_flux.gy), 'Color', plt_clrs.blue);
ylabel(label.met_reac.snf1{num_flux.gy}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.snf1(:,num_flux.gy)); 
ymax2 = max(WT_met_reac_steady_10gL.snf1(:,num_flux.gy)); 
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
ylabel(label.met{num_y.eh})    
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



fig_index = fig_index + b; 

%--------------------------------------------------------------------------
% ethanol

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
plot(D, WT_met_reac_steady_1gL.snf1(:,num_flux.gn), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.snf1(:,num_flux.gn), 'Color', plt_clrs.blue);
ylabel(label.met_reac.snf1{num_flux.gn}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.snf1(:,num_flux.gn)); 
ymax2 = max(WT_met_reac_steady_10gL.snf1(:,num_flux.gn)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
%xline(t_gl_depl)
%xline(t_eh_depl)
axis square; 
box on; 
hold off; 

%}
%sgtitle('real metabolites')

%% virtual metabolites - precursor 
act_no_sub = 1; 
figure; 
a = 5+3;
b = 7;
fig_index = 0:6;
% fig_index(8) = 1; 
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
plot(D, WT_met_reac_steady_1gL.flux(:,num_flux.gy), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.gy), 'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.gy})    
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
ylim(y_lim_flux)
xlim(x_lim)
xticks(x_ticks)
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
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

activity = WT_met_reac_steady_1gL.substrate(:,num_flux.gy).*WT_met_reac_steady_1gL.snf1(:,num_flux.gy).* WT_met_reac_steady_10gL.atp(:,num_flux.gy);
if act_no_sub == 1
activity = 1.*WT_met_reac_steady_1gL.sig(:,num_flux.gy).* WT_met_reac_steady_10gL.atp(:,num_flux.gy);
activity2 = 1.*WT_met_reac_steady_10gL.sig(:,num_flux.gy).* WT_met_reac_steady_10gL.atp(:,num_flux.gy);
end 
subplot(a, b, fig_index(1)); 
hold on; 
plot(D, activity, 'Color', plt_clrs.green);
plot(D, activity2, 'Color', plt_clrs.blue);
ylabel('activity') 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 

% k*Gl/(km+Gl) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.substrate(:,num_flux.gy), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.substrate(:,num_flux.gy), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.gy}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.snf1(:,num_flux.gy), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.snf1(:,num_flux.gy), 'Color', plt_clrs.blue);
ylabel(label.met_reac.snf1{num_flux.gy}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 


subplot(a, b, fig_index(7)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.pc(:,num_flux.gy), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.pc(:,num_flux.gy), 'Color', plt_clrs.blue);
ylabel(label.met_reac.pc{num_flux.gy}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
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
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 


fig_index = fig_index + b; 
%--------------------------------------------------------------------------

% J_as
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, 2*WT_met_reac_steady_1gL.flux(:,num_flux.as), 'Color', plt_clrs.green);
plot(D, 2*WT_met_reac_steady_10gL.flux(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.as})    
xlabel(x_label) 
ymax1 = max(2*WT_met_reac_steady_1gL.flux(:,num_flux.as)); 
ymax2 = max(2*WT_met_reac_steady_10gL.flux(:,num_flux.as)); 
ylim(y_lim_flux)
xlim(x_lim)
xticks(x_ticks)
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
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

activity = WT_met_reac_steady_1gL.substrate(:,num_flux.as).*WT_met_reac_steady_1gL.snf1(:,num_flux.as).* WT_met_reac_steady_10gL.atp(:,num_flux.as);
if act_no_sub == 1
activity = 1.*WT_met_reac_steady_1gL.sig(:,num_flux.as).* WT_met_reac_steady_1gL.atp(:,num_flux.as);
activity2 = 1.*WT_met_reac_steady_10gL.sig(:,num_flux.as).* WT_met_reac_steady_10gL.atp(:,num_flux.as);

end 
subplot(a, b, fig_index(1)); 
hold on; 
plot(D, activity, 'Color', plt_clrs.green);
plot(D, activity2, 'Color', plt_clrs.blue);
ylabel('activity') 
xlabel(x_label)
ylim([0 1.1*max(ylim)])
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 

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
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(7)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.pc(:,num_flux.as), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.pc(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.pc{num_flux.as}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.snf1(:,num_flux.as)); 
ymax2 = max(WT_met_reac_steady_10gL.snf1(:,num_flux.as)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
xlim(x_lim)
xticks(x_ticks)
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
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

fig_index = fig_index + b; 

% J_mt
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.flux(:,num_flux.mt), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.mt), 'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.mt})    
xlabel(x_label)
ylim(y_lim_flux)
xlim(x_lim)
xticks(x_ticks)
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
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 


activity = WT_met_reac_steady_1gL.substrate(:,num_flux.mt).*WT_met_reac_steady_1gL.snf1(:,num_flux.mt).* WT_met_reac_steady_10gL.atp(:,num_flux.mt);
if act_no_sub == 1
activity = 1*WT_met_reac_steady_1gL.sig(:,num_flux.mt).* WT_met_reac_steady_1gL.atp(:,num_flux.mt);
activity2 = 1*WT_met_reac_steady_10gL.sig(:,num_flux.mt).* WT_met_reac_steady_10gL.atp(:,num_flux.mt);

end
subplot(a, b, fig_index(1)); 
hold on; 
plot(D, activity, 'Color', plt_clrs.green);
plot(D, activity2, 'Color', plt_clrs.blue);
ylabel('activity') 
xlabel(x_label)
ylim([0 1.1*max(ylim)])
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 

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
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(7)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.pc(:,num_flux.mt), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.pc(:,num_flux.mt), 'Color', plt_clrs.blue);
ylabel(label.met_reac.pc{num_flux.mt}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.snf1(:,num_flux.mt)); 
ymax2 = max(WT_met_reac_steady_10gL.snf1(:,num_flux.mt)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
xlim(x_lim)
xticks(x_ticks)
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
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

fig_index = fig_index + b; 
%--------------------------------------------------------------------------

% J_fe - J_gn 
% subplot(a, b, fig_index(8)); 
% hold on; 
% plot(D, WT_met_reac_steady_1gL.flux(:,num_flux.fe) - WT_met_reac_steady_1gL.flux(:,num_flux.gn), 'Color', plt_clrs.green);
% plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.fe) - WT_met_reac_steady_10gL.flux(:,num_flux.gn), 'Color', plt_clrs.blue);
% ylabel('J_{eh}')    
% xlabel(x_label)
% ylim(y_lim_flux)
% xlim(x_lim)
% xticks(x_ticks)
% axis square; 
% box on; 
% hold off; 

% J_fe 
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.flux(:,num_flux.fe), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.fe})    
xlabel(x_label)
ylim(y_lim_flux)
xlim(x_lim)
xticks(x_ticks)
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
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

activity = WT_met_reac_steady_1gL.substrate(:,num_flux.fe).*WT_met_reac_steady_1gL.snf1(:,num_flux.fe).* WT_met_reac_steady_10gL.atp(:,num_flux.fe);
if act_no_sub == 1
activity = 1.*WT_met_reac_steady_1gL.sig(:,num_flux.fe).* WT_met_reac_steady_1gL.atp(:,num_flux.fe);
activity2 = 1.*WT_met_reac_steady_10gL.sig(:,num_flux.fe).* WT_met_reac_steady_10gL.atp(:,num_flux.fe);

end 
subplot(a, b, fig_index(1)); 
hold on; 
plot(D, activity, 'Color', plt_clrs.green);
plot(D, activity2, 'Color', plt_clrs.blue);
ylabel('activity') 
xlabel(x_label)
ylim([0 1.1*max(ylim)])
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 

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
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

%{
% signaling contribution 
subplot(a, b, fig_index(5)); 
plot(bart_WT_batch_min_t, bart_WT_batch_min_met_reac.snf1(:,num_flux.fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.snf1{num_flux.fe}) 
xlabel(x_label)
ylim([0 1.05.*max(bart_WT_batch_min_met_reac.snf1(:,num_flux.fe))])
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

% J_gn
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.flux(:,num_flux.gn), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.gn), 'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.gn})    
xlabel(x_label)
ylim(y_lim_flux)
xlim(x_lim)
xticks(x_ticks)
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
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

activity = WT_met_reac_steady_1gL.substrate(:,num_flux.gn).*WT_met_reac_steady_1gL.snf1(:,num_flux.gn).* WT_met_reac_steady_10gL.atp(:,num_flux.gn);
if act_no_sub == 1
activity = 1.*WT_met_reac_steady_1gL.sig(:,num_flux.gn).* WT_met_reac_steady_1gL.atp(:,num_flux.gn);
activity2 = 1.*WT_met_reac_steady_10gL.sig(:,num_flux.gn).* WT_met_reac_steady_10gL.atp(:,num_flux.gn);

end 
subplot(a, b, fig_index(1)); 
hold on; 
plot(D, activity, 'Color', plt_clrs.green);
plot(D, activity2, 'Color', plt_clrs.blue);
ylabel('activity') 
xlabel(x_label)
ylim([0 1.1*max(ylim)])
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 

% k*Gl/(km+Gl) - substrate contribution 
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
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.snf1(:,num_flux.gn), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.snf1(:,num_flux.gn), 'Color', plt_clrs.blue);
ylabel(label.met_reac.snf1{num_flux.gn}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.snf1(:,num_flux.gn)); 
ymax2 = max(WT_met_reac_steady_10gL.snf1(:,num_flux.gn)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

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
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

fig_index = fig_index + b; 
%{
% precursor
subplot(a, b, fig_index(1));
plot(bart_WT_batch_min_t, y(:,3), 'Color', plt_clrs.blue);
ylabel(label.met(3))    
xlabel(x_label)
axis square; 
box on;
%}

%--------------------------------------------------------------------------
% Jlp,fe
% J_lp,fe
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.flux(:,num_flux.lp_fe), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.lp_fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.lp_fe})    
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
ylim(y_lim_flux)
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

% E_lp,fe
subplot(a, b, fig_index(3)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.prot(:,num_flux.lp_fe), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.prot(:,num_flux.lp_fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.lp_fe})
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

activity = WT_met_reac_steady_1gL.substrate(:,num_flux.lp_fe).*WT_met_reac_steady_1gL.snf1(:,num_flux.lp_fe).* WT_met_reac_steady_10gL.atp(:,num_flux.lp_fe);
if act_no_sub == 1
activity = 1.*WT_met_reac_steady_1gL.sig(:,num_flux.lp_fe).* WT_met_reac_steady_10gL.atp(:,num_flux.lp_fe);
activity2 = 1.*WT_met_reac_steady_10gL.sig(:,num_flux.lp_fe).* WT_met_reac_steady_10gL.atp(:,num_flux.lp_fe);
end 
subplot(a, b, fig_index(1)); 
hold on; 
plot(D, activity, 'Color', plt_clrs.green);
plot(D, activity2, 'Color', plt_clrs.blue);
ylabel('activity') 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 

% k*Gl/(km+Gl) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.substrate(:,num_flux.lp_fe), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.substrate(:,num_flux.lp_fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.lp_fe}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.snf1(:,num_flux.lp_fe), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.snf1(:,num_flux.lp_fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.snf1{num_flux.lp_fe}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 


subplot(a, b, fig_index(7)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.pc(:,num_flux.lp_fe), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.pc(:,num_flux.lp_fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.pc{num_flux.lp_fe}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 


% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.atp(:,num_flux.lp_fe), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.atp(:,num_flux.lp_fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.atp{num_flux.lp_fe}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 


fig_index = fig_index + b; 

% -------------------------------------------------------------------------
% J_lp,cit
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.flux(:,num_flux.lp_cit), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.lp_cit), 'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.lp_cit})    
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
ylim(y_lim_flux)
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

% E_lp
subplot(a, b, fig_index(3)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.prot(:,num_flux.lp_cit), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.prot(:,num_flux.lp_cit), 'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.lp_cit})
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

activity = WT_met_reac_steady_1gL.substrate(:,num_flux.lp_cit).*WT_met_reac_steady_1gL.snf1(:,num_flux.lp_cit).* WT_met_reac_steady_10gL.atp(:,num_flux.lp_cit);
if act_no_sub == 1
activity = 1.*WT_met_reac_steady_1gL.sig(:,num_flux.lp_cit).* WT_met_reac_steady_10gL.atp(:,num_flux.lp_cit);
activity2 = 1.*WT_met_reac_steady_10gL.sig(:,num_flux.lp_cit).* WT_met_reac_steady_10gL.atp(:,num_flux.lp_cit);
end 
subplot(a, b, fig_index(1)); 
hold on; 
plot(D, activity, 'Color', plt_clrs.green);
plot(D, activity2, 'Color', plt_clrs.blue);
ylabel('activity') 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 

% k*Gl/(km+Gl) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.substrate(:,num_flux.lp_cit), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.substrate(:,num_flux.lp_cit), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.lp_cit}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.snf1(:,num_flux.lp_cit), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.snf1(:,num_flux.lp_cit), 'Color', plt_clrs.blue);
ylabel(label.met_reac.snf1{num_flux.lp_cit}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 


subplot(a, b, fig_index(7)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.pc(:,num_flux.lp_cit), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.pc(:,num_flux.lp_cit), 'Color', plt_clrs.blue);
ylabel(label.met_reac.pc{num_flux.lp_cit}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 


% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.atp(:,num_flux.lp_cit), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.atp(:,num_flux.lp_cit), 'Color', plt_clrs.blue);
ylabel(label.met_reac.atp{num_flux.lp_cit}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 


fig_index = fig_index + b; 

% -------------------------------------------------------------------------
% J_lo
subplot(a, b, fig_index(2)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.flux(:,num_flux.lo), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.flux(:,num_flux.lo), 'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.lo})    
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
ylim(y_lim_flux)
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

% E_lp
subplot(a, b, fig_index(3)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.prot(:,num_flux.lo), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.prot(:,num_flux.lo), 'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.lo})
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

activity = WT_met_reac_steady_1gL.substrate(:,num_flux.lo).*WT_met_reac_steady_1gL.snf1(:,num_flux.lo).* WT_met_reac_steady_10gL.atp(:,num_flux.lo);
if act_no_sub == 1
activity = 1.*WT_met_reac_steady_1gL.sig(:,num_flux.lo).* WT_met_reac_steady_10gL.atp(:,num_flux.lo);
activity2 = 1.*WT_met_reac_steady_10gL.sig(:,num_flux.lo).* WT_met_reac_steady_10gL.atp(:,num_flux.lo);
end 
subplot(a, b, fig_index(1)); 
hold on; 
plot(D, activity, 'Color', plt_clrs.green);
plot(D, activity2, 'Color', plt_clrs.blue);
ylabel('activity') 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 

% k*Gl/(km+Gl) - substrate contribution 
subplot(a, b, fig_index(4)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.substrate(:,num_flux.lo), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.substrate(:,num_flux.lo), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.lo}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, fig_index(5)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.snf1(:,num_flux.lo), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.snf1(:,num_flux.lo), 'Color', plt_clrs.blue);
ylabel(label.met_reac.snf1{num_flux.lo}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 


subplot(a, b, fig_index(7)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.pc(:,num_flux.lo), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.pc(:,num_flux.lo), 'Color', plt_clrs.blue);
ylabel(label.met_reac.pc{num_flux.lo}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 


% ATP contribution 
subplot(a, b, fig_index(6)); 
hold on; 
plot(D, WT_met_reac_steady_1gL.atp(:,num_flux.lo), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.atp(:,num_flux.lo), 'Color', plt_clrs.blue);
ylabel(label.met_reac.atp{num_flux.lo}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

sgtitle('precursor')


%% addtional figures for newly added variables
figure; 
a = 4 ;b = 5; 

% J_as
subplot(a, b, 1); 
hold on; 
plot(D, 2*WT_met_reac_steady_1gL.flux(:,num_flux.as), 'Color', plt_clrs.green);
plot(D, 2*WT_met_reac_steady_10gL.flux(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.as})    
xlabel(x_label) 
ymax1 = max(2*WT_met_reac_steady_1gL.flux(:,num_flux.as)); 
ymax2 = max(2*WT_met_reac_steady_10gL.flux(:,num_flux.as)); 
ylim(y_lim_flux)
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

% E_as
subplot(a, b, 2); 
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
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 




% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, 3); 
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
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

% signaling contribution 
subplot(a, b, 4); 
hold on; 
plot(D, WT_met_reac_steady_1gL.pc(:,num_flux.as), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.pc(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.pc{num_flux.as}) 
xlabel(x_label)
ymax1 = max(WT_met_reac_steady_1gL.snf1(:,num_flux.as)); 
ymax2 = max(WT_met_reac_steady_10gL.snf1(:,num_flux.as)); 
if max([ymax1, ymax2]) > 0 
ylim([0 1.05*max([ymax1, ymax2])])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on;
hold off; 

% ATP contribution 
subplot(a, b, 5); 
hold on; 
plot(D, WT_met_reac_steady_1gL.atp(:,num_flux.as), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.atp(:,num_flux.as), 'Color', plt_clrs.blue);
ylabel(label.met_reac.atp{num_flux.as}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 



% Pc part in Jas substrate 
subplot(a, b, 3+b); 
hold on; 
plot(D, as_sub_term_1gL, 'Color', plt_clrs.green);
plot(D, as_sub_term_10gL, 'Color', plt_clrs.blue);
ylabel({'[P_{c}]/' , '(K+[P_{c}])'}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on;
hold off; 


% NH4 part in Jas substrate 
subplot(a, b, 4+b); 
hold on; 
plot(D, as_nh4_term_1gL, 'Color', plt_clrs.green);
plot(D, as_nh4_term_10gL, 'Color', plt_clrs.blue);
ylabel({'[NH_{4}]/' , '(K+[NH_{4}])'}) 
xlabel(x_label)
if max(ylim) > 0 
ylim([0 1.05*max(ylim)])
end
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on;
hold off; 



% ------------------------------------------------------------------------
% Jlp,fe 

subplot(a, b, 1+2*b); 
hold on; 
plot(D, 2*WT_met_reac_steady_1gL.flux(:,num_flux.lp_fe), 'Color', plt_clrs.green);
plot(D, 2*WT_met_reac_steady_10gL.flux(:,num_flux.lp_fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.lp_fe})    
xlabel(x_label) 
if max(ylim) > 0
ylim([0 max(ylim)])
end 
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

% E_lp
subplot(a, b, 2+2*b); 
hold on; 
plot(D, WT_met_reac_steady_1gL.prot(:,num_flux.lp_fe), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.prot(:,num_flux.lp_fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.lp_fe})
xlabel(x_label)
if max(ylim) > 0
ylim([0 max(ylim)])
end 
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 


% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, 3+2*b); 
hold on; 
plot(D, WT_met_reac_steady_1gL.substrate(:,num_flux.lp_fe), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.substrate(:,num_flux.lp_fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.lp_fe}) 
xlabel(x_label)
if max(ylim) > 0
ylim([0 max(ylim)])
end 
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

% SNF1 contribution
subplot(a, b, 4+2*b); 
hold on; 
plot(D, WT_met_reac_steady_1gL.snf1(:,num_flux.lp_fe), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.snf1(:,num_flux.lp_fe), 'Color', plt_clrs.blue);
ylabel(label.met_reac.snf1{num_flux.lp_fe}) 
xlabel(x_label)
if max(ylim) > 0
ylim([0 max(ylim)])
end 
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on;
hold off; 



% ------------------------------------------------------------------------
% Jlp,cit 

subplot(a, b, 1+3*b); 
hold on; 
plot(D, 2*WT_met_reac_steady_1gL.flux(:,num_flux.lp_cit), 'Color', plt_clrs.green);
plot(D, 2*WT_met_reac_steady_10gL.flux(:,num_flux.lp_cit), 'Color', plt_clrs.blue);
ylabel(label.met_reac.flux{num_flux.lp_cit})    
xlabel(x_label) 
if max(ylim) > 0
ylim([0 max(ylim)])
end 
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

% E_lp
subplot(a, b, 2+3*b); 
hold on; 
plot(D, WT_met_reac_steady_1gL.prot(:,num_flux.lp_cit), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.prot(:,num_flux.lp_cit), 'Color', plt_clrs.blue);
ylabel(label.met_reac.prot{num_flux.lp_cit})
xlabel(x_label)
if max(ylim) > 0
ylim([0 max(ylim)])
end 
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 


% k*Pc/(km+Pc) - substrate contribution 
subplot(a, b, 3+3*b); 
hold on; 
plot(D, WT_met_reac_steady_1gL.substrate(:,num_flux.lp_cit), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.substrate(:,num_flux.lp_cit), 'Color', plt_clrs.blue);
ylabel(label.met_reac.substrate{num_flux.lp_cit}) 
xlabel(x_label)
if max(ylim) > 0
ylim([0 max(ylim)])
end 
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on; 
hold off; 

% SNF1 contribution
subplot(a, b, 4+3*b); 
hold on; 
plot(D, WT_met_reac_steady_1gL.snf1(:,num_flux.lp_cit), 'Color', plt_clrs.green);
plot(D, WT_met_reac_steady_10gL.snf1(:,num_flux.lp_cit), 'Color', plt_clrs.blue);
ylabel(label.met_reac.snf1{num_flux.lp_fe}) 
xlabel(x_label)
if max(ylim) > 0
ylim([0 max(ylim)])
end 
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on;
hold off; 


% SNF1 contribution
subplot(a, b, 5+3*b); 
hold on; 
plot(D, lp_aa_term_1gL, 'Color', plt_clrs.green);
plot(D, lp_aa_term_10gL, 'Color', plt_clrs.blue);
ylabel('N term') 
xlabel(x_label)
if max(ylim) > 0
ylim([0 max(ylim)])
end 
xlim(x_lim)
xticks(x_ticks)
axis square; 
box on;
hold off; 




end 
