% FA Hackett20161
% collected from C-limit and N-limit condition
% original unit: fraction of biomass
close all; clc;
fa = [0.006678443	0.007197644	0.007179225	0.005887776	0.004491896	0.004802612	0.003640444	0.005052841	0.004394802	0.003610893	0.01079081	0.00957392	0.007205246	0.00585251	0.004900328	0.008732372	0.006762352	0.00625494	0.004498592	0.00418617	0.006840384	0.005997257	0.006658933	0.005248563	0.007889887;
0.018560812	0.020237484	0.017532586	0.018984058	0.017110032	0.016148111	0.016473037	0.021493343	0.018142027	0.016391696	0.031990649	0.029861104	0.02058455	0.01862279	0.015339959	0.012734697	0.014090408	0.013382197	0.012501266	0.013898508	0.014582031	0.01703861	0.01645267	0.015061376	0.02115657;
0.002381052	0.002646855	0.005168672	0.001180821	0.000595741	0.001505975	0.000469034	0.001061315	0.001081977	0.000362134	0.004722919	0.003694329	0.002479982	0.001457323	0.000338824	0.002970217	0.001633099	0.001489654	0.000538667	0.000695772	0.000879597	0.001029502	0.002434681	0.000467144	0.00070916;
0.014837911	0.012965456	0.009812459	0.008463637	0.00662838	0.008356945	0.006557496	0.009023905	0.009078053	0.007889476	0.024391808	0.016711511	0.011231336	0.008095156	0.006572411	0.008306324	0.008045779	0.007064293	0.00629675	0.006582799	0.005016967	0.007161093	0.009677953	0.005559954	0.008204764;];

fa = fa.*100;
fa_tot = sum(fa);

figure();
subplot(2,1,1)
x_tot = [2, 8];
plot(fa_tot, fa(1,:), "o", 'Color',"#0072BD", 'MarkerFaceColor',"#0072BD",'MarkerSize',8); hold on %blue
plot(fa_tot, fa(2,:), "o", 'Color',"#77AC30", 'MarkerFaceColor',"#77AC30",'MarkerSize',8); % green
plot(fa_tot, fa(3,:), "o", 'Color',"#EDB120", 'MarkerFaceColor',"#EDB120",'MarkerSize',8); % orange
plot(fa_tot, fa(4,:), "o", 'Color',"#D95319", 'MarkerFaceColor',"#D95319",'MarkerSize',8); % red
plot(x_tot, 0.5 + 0.13.* (x_tot-2.4), "--", 'Color',"#0072BD",'LineWidth',2)
plot(x_tot, 1.3 + 0.40.* (x_tot-2.4), "--", 'Color',"#77AC30",'LineWidth',2)
plot(x_tot, 0.1 + 0.10.* (x_tot-2.4), "--", 'Color',"#EDB120",'LineWidth',2)
plot(x_tot, 0.5 + 0.37.* (x_tot-2.4), "--", 'Color',"#D95319",'LineWidth',2)
legend("C16:0", "C16:1", "C18:0", "C18:1", 'FontSize', 12, 'Location', 'southoutside')
ylim([0, 35]); 
xlim([3, 65]); hold off
xlabel("[FA_{tot}] (% of DW)")
ylabel("[FA_{i}] (% of DW)")
% axis square

subplot(2,1,2)
x0_tot = [2, 2.4];
x_tot = [2.4, 8];
plot(fa_tot, fa(1,:), "o", 'Color',"#0072BD", 'MarkerFaceColor',"#0072BD",'MarkerSize',8); hold on %blue
plot(fa_tot, fa(2,:), "o", 'Color',"#77AC30", 'MarkerFaceColor',"#77AC30",'MarkerSize',8); % green
plot(fa_tot, fa(3,:), "o", 'Color',"#EDB120", 'MarkerFaceColor',"#EDB120",'MarkerSize',8); % orange
plot(fa_tot, fa(4,:), "o", 'Color',"#D95319", 'MarkerFaceColor',"#D95319",'MarkerSize',8); % red
plot(x_tot, 0.5 + 0.13.* (x_tot-2.4), "--", 'Color',"#0072BD",'LineWidth',2)
plot(x_tot, 1.3 + 0.40.* (x_tot-2.4), "--", 'Color',"#77AC30",'LineWidth',2)
plot(x_tot, 0.1 + 0.10.* (x_tot-2.4), "--", 'Color',"#EDB120",'LineWidth',2)
plot(x_tot, 0.5 + 0.37.* (x_tot-2.4), "--", 'Color',"#D95319",'LineWidth',2)
plot(x0_tot, [0.5, 0.5] , "--", 'Color',"#0072BD",'LineWidth',2)
plot(x0_tot, [1.3, 1.3] , "--", 'Color',"#77AC30",'LineWidth',2)
plot(x0_tot, [0.1, 0.1] , "--", 'Color',"#EDB120",'LineWidth',2)
plot(x0_tot, [0.5, 0.5] , "--", 'Color',"#D95319",'LineWidth',2)
xline(2.4, "k--")
legend( ...
    "[FA_{16:0}] = 0.5 + 0.13 ([FA_{tot}]-2.4)", ...
    "[FA_{16:1}] = 1.3 + 0.40 ([FA_{tot}]-2.4)", ...
    "[FA_{18:0}] = 0.1 + 0.10 ([FA_{tot}]-2.4)", ...
    "[FA_{18:1}] = 0.5 + 0.37 ([FA_{tot}]-2.4)",...
    'FontSize', 12, 'Location', 'southoutside')
ylim([0, 4]); 
xlim([2, 8]); hold off
% xlim([3, 8]); hold off
xlabel("[FA_{tot}] (% of DW)")
ylabel("[FA_{i}] (% of DW)")
% axis square
fontsize(gcf, 14, 'points')
set(gcf, 'Position', [659 92 309 885])
