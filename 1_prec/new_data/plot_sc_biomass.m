close all; clear; clc; 
data_unit_conv;
plt_colors = plot_colors; 

%	protein	    glycogen	trehalose	total_sc	glucan/mannan	lipid	RNA		
Albers2007 = [47.2	3.1	    2.8	   5.9	   30.6	  4.9	5.1	; %	Diauxic shift
48.9	0.32	4	   4.32	   34.2	  4	     4.8	; %	C-starved
45.6	0.33	3.2	   3.53	   31.3	  3.3	4.3	; %	N and C
24.9	11.5	14.3	25.8	33	   10.1	2.2	; %	N starved glucsoe
35.4	4.1	     9.4	13.5	31.9	7.6	 2.9	;]; %	N starved ethanol

%	protein	carbohyd	lipid	RNA
Lange_2001 = [34.5	   44.9	    8	3.8	; %	D=0.022
38.5	40.9	7.8	4.6	; %	D=0.052
39.7	40.6	6.2	5.4	; %	D=0.087
40.7	38.7	6.1	6.6	; %	D=0.107
40.1	38.8	7	5.7	; %	D=0.126
42.2	35.9	7	5.9	; %	D=0.158
43.8	34.4	6.7	6.6	;]; %	D=0.211

%	protein	glycogen	trehalose	total_sc	lipid	RNA
Nissen_1997 = [45	    8.4	0.8	9.2	2.9	6.3	        ; %	D=0.1
50	    4.2	0.2	4.4	3	8.2	    ; %	D=0.2
55.5	0.6	0	0.6	3.8	10.1	; %	D=0.3
60.1	0	0	0	3.4	12.1	;]; %	D=0.4

D_Hackett    = [0.05 0.11 0.16 0.22 0.3];
biomass_density_Hackett_clim    = [0.204; 0.205;0.218;0.211; 0.187;];
% Carbohydrate	Protein	RNA	DNA	Polyphosphates	FattyAcids	glycogen	trehalose	total_sc
Hackett_clim = [42.79427713	43.10167818	5.160876547	0.464340759	4.708188862	3.77063852	6.517529816	11.85706521	18.37459502	; %	C_0.05
35.93818917	48.33600798	6.624288341	0.52956661	5.323059838	3.248888055	5.76970123	8.550747192	14.32044842	; %	C_0.11
36.30855581	46.18232506	6.974462457	0.465048489	5.738418926	4.331189253	6.452798984	5.678606228	12.13140521	; %	C_0.16
36.19544632	46.49293565	7.731532192	0.524790945	5.157973508	3.897321384	6.705372236	4.366589715	11.07196195	; %	C_0.22
34.03065389	47.45531735	8.800270034	0.586377107	5.758113529	3.369268084	6.880691847	1.369595508	8.250287355	;]; %	C_0.30]

biomass_density_Hackett_nlim    = [0.216;0.2;0.195;0.239; 0.182;];
% Carbohydrate	Protein	RNA	DNA	Polyphosphates	FattyAcids	glycogen	trehalose	total_sc
Hackett_nlim = [43.66484191	33.96302042	3.603473697	0.355324043	10.43438044	7.978959488	3.323467035	27.88910986	31.21257689	; %	N_0.05
40.49655727	38.46475273	4.870739075	0.504143075	9.071818206	6.591989642	5.725404807	13.3193802	19.04478501	; %	N_0.11
38.64183387	40.52640718	5.672433086	0.502298687	9.907876022	4.749151161	6.0173894	10.07867323	16.09606263	; %	N_0.16
36.2278625	43.90136579	7.011617104	0.395332625	8.416022966	4.047799023	7.117478175	2.442833981	9.560312156	; %	N_0.22
35.50744616	43.59564355	7.324368902	0.531537991	9.530437642	3.510565751	6.966058985	2.441169649	9.407228633	;]; %	N_0.30]

% ------------------------------------------------------------------------
% 1; RNA vs protein
% 2: lipid vs storage carbon
figure; 
a = 1; b = 2; 
subplot(a,b,1)
hold on 
plot(Albers2007(:,1), Albers2007(:,end),'o')
plot(Lange_2001(:,1), Lange_2001(:,end),'o')
plot(Nissen_1997(:,1), Nissen_1997(:,end),'o')
plot(Hackett_clim(:,2), Hackett_clim(:,3),'o')
plot(Hackett_nlim(:,2), Hackett_nlim(:,3),'o')
xlabel('Protein (%)')
ylabel('RNA (%)')
xlim([10 70])
ylim([0 14])
hold off
legend('Albers 2007', 'Lange 2001', 'Nissen 1997', 'Hackett clim', 'Hackett nlim')
legend box off
box on 
axis square

subplot(a,b,2)
hold on 
plot(Albers2007(:,4), Albers2007(:,end-1),'o')
%plot(Lange_2001(:,1), Lange_2001(:,end-1),'o')
plot(Nissen_1997(:,4), Nissen_1997(:,end-1),'o')
plot(Hackett_clim(:,end), Hackett_clim(:,end-3),'o')
plot(Hackett_nlim(:,end), Hackett_nlim(:,end-3),'o')

% plot(Albers2007(:,4), Albers2007(:,end-1),'o', 'Color', plt_colors.blue)
% %plot(Lange_2001(:,1), Lange_2001(:,end-1),'o')
% plot(Nissen_1997(:,4), Nissen_1997(:,end-1),'o', 'Color',plt_colors.green)
% plot(Hackett_clim(:,end), Hackett_clim(:,end-3),'o','Color', plt_colors.green)
% plot(Hackett_nlim(:,end), Hackett_nlim(:,end-3),'o','Color', plt_colors.blue)

xlabel('storage carbon (%)')
ylabel('lipid (%)')
xlim([0 30])
ylim([2 12])
hold off
legend('Albers 2007',  'Nissen 1997', 'Hackett clim', 'Hackett nlim')
legend box off
box on 
axis square
set(gcf, 'Position', [440 434 263 364])

% ------------------------------------------------------------------------
% protein, RNA vs D 
figure; 
a = 1; b = 2; 
subplot(a,b,1)
hold on 
plot(D_Hackett, Hackett_clim(:,2), 'o', 'Color',plt_colors.green)
plot(D_Hackett, Hackett_clim(:,3), 'o', 'Color',plt_colors.yellow)
hold off 
ylim([0 60])
xlim([0. 0.35])
xlabel('D (h^{-1})')
ylabel('%')
legend('protein', 'RNA')
legend box off
title('Clim chemostat')
box on 
axis square


subplot(a,b,2)
hold on 
plot(D_Hackett, Hackett_nlim(:,2),  'o', 'Color',plt_colors.green)
plot(D_Hackett, Hackett_nlim(:,3),  'o', 'Color',plt_colors.yellow)
hold off 
ylim([0 60])
xlim([0. 0.35])
xlabel('D (h^{-1})')
ylabel('%')
legend('protein', 'RNA')
legend box off
title('Nlim chemostat')
box on 
axis square

set(gcf, 'Position', [440 434 263 364])

% ------------------------------------------------------------------------
% storage carbon, lipid vs D 
figure; 
a = 1; b = 2; 
subplot(a,b,1)
hold on 
plot(D_Hackett, Hackett_clim(:,end), 'o')
plot(D_Hackett, Hackett_clim(:,end-3),'o')
hold off 
ylim([0 36])
xlim([0. 0.35])
xlabel('D (h^{-1})')
ylabel('%')
legend('storage carbon', 'lipid')
legend box off
title('Clim chemostat')
box on 
axis square


subplot(a,b,2)
hold on 
plot(D_Hackett, Hackett_nlim(:,end),  'o')
plot(D_Hackett, Hackett_nlim(:,end-3),'o')
hold off 
ylim([0 36])
xlim([0. 0.35])
xlabel('D (h^{-1})')
ylabel('%')
legend('storage carbon', 'lipid')
legend box off
title('Nlim chemostat')
box on 
axis square

set(gcf, 'Position', [440 434 263 364])





% ------------------------------------------------------------------------
% storage carbon, lipid vs D 
figure; 
a = 1; b = 2; 
subplot(a,b,1)
hold on 
plot(D_Hackett,biomass_density_Hackett_clim, 'ko')
hold off 
ylim([0 0.26])
xlim([0. 0.35])
xlabel('D (h^{-1})')
ylabel('gDW/ml')

title('Clim chemostat')
box on 
axis square


subplot(a,b,2)
hold on 
plot(D_Hackett, biomass_density_Hackett_nlim,  'ko')
hold off 
ylim([0 0.26])
xlim([0. 0.35])
xlabel('D (h^{-1})')
ylabel('gDW/ml')
title('Nlim chemostat')
box on 
axis square

set(gcf, 'Position', [440 434 263 364])
