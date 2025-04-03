%% plotting set up 
plt_clrs     = plot_colors;
xy_label_color  = 'k';

%xy_label_color  = 'k';
%%
% patch two region color
color1 = [17  17  17]/255; % lighter gray      
color2 = [1   1   1]/255;  % darker gray
alpha1 = 0.1;              
alpha2 = 0.15;  

% %% ----------------------------------------------------------------------------------------------------------------
% %                               growth opt colors 
% % ----------------------------------------------------------------------------------------------------------------- 
% 
% color_camp1  = [158,   1,  66]./255;
% color_camp2  = [213,  62,  79]./255;
% color_camp3  = [244, 109,  67]./255;
% color_camp4  = [253, 174,  97]./255;
% color_camp5  = [254, 224, 139]./255;
% color_camp6  = [230, 245, 152]./255;
% color_camp7  = [171, 221, 164]./255;
% color_camp8  = [102, 194, 165]./255;
% color_camp9  = [50,  136, 189]./255;
% color_camp10 = [30,   82, 113]./255;
% color_camp11 = [127,115,180]./255 ; 
% color_camp12 = [94,   79, 162]./255; 
% 
% camp_colors = [color_camp1;
% color_camp2;
% color_camp3; 
% color_camp4;
% color_camp5;
% color_camp6;
% color_camp7;
% color_camp8;
% color_camp9;
% color_camp10;];
% 
% color_wt     = 'k';
% 
% % background for figure 5,6 -----------------------------------------------
% color_low_camp  = '#F8F5F1'; % white
% color_high_camp = '#C7EAFB'; %'#EAE0D6';  % darker gray
% color_high_wt   = '#EAE0D6';%'#FDFAD9'; 
% 
% % cAMP colors -------------------------------------------------------------
% 
% % open-loop max color: 
% camp_gr_max   = plt_clrs.blue;               % gr. and cells
% camp_cell_max = plt_clrs.red;
% 
% 
% % GCN2 colors -------------------------------------------------------
% glu_gnc2_color = plt_clrs.lightgray;%plt_clrs.lightblue;
% wt_gnc2_color = 'k';% plt_clrs.lightblue;
% 
% eth_gnc2_color = plt_clrs.lightgray;%plt_clrs.lightblue;
% wt_eth_gnc2_color = 'k';% plt_clrs.lightblue;
% 
% glu_color_max = 'k';              % sugars 
% eth_color_max = plt_clrs.blue; 
% fru_color_max = plt_clrs.gray;   
% suc_color_max = plt_clrs.orange; 
% 
% %% Gr opt Figure 1 color 
% plt_clrs.orange      = '#fbb03a'; 
% 
% % orange --> white 
% orange_light_gradients = {'#fbb03a';
% '#ffbd5f';
% '#ffcb82';
% '#ffd9a3';
% '#ffe6c3';
% '#fff2e1';
% '#ffffff';};
% % orange --> black 
% orange_dark_gradients = {'#fbb03a';
% '#cd9033';
% '#a1722b';
% '#775423';
% '#4f391b';
% '#2a1f12';
% '#000000';};
% 
% % green --> white 
% green_light_gradients = {'#4cbd94';
% '#70c8a5';
% '#8fd4b7';
% '#acdfc8';
% '#c8eada';
% '#e3f4ed';
% '#ffffff';};
% % green --> black 
% green_dark_gradients = {'#4cbd94';
% '#419b7a';
% '#377a60';
% '#2c5a48';
% '#203d31';
% '#15211c';
% '#000000';};
% 
% % blue --> white 
% blue_light_gradients = {'#2078b4';
% '#568dc1';
% '#7ca3cd';
% '#9eb9da';
% '#bfd0e6';
% '#dfe7f3';
% '#ffffff';};
% 
% % blue --> black 
% blue_dark_gradients = {'#2078b4'
% '#226393'
% '#214f74'
% '#1e3c56'
% '#18293a'
% '#111820'
% '#000000'};
% 
% 
% % red --> white 
% red_light_gradients = {'#e31f26';
% '#ef5547';
% '#f97a69';
% '#ff9d8d';
% '#F9B7B7';
% '#ffe0d9';
% '#ffffff';};
% 
% % red --> black 
% red_dark_gradients = {'#e31f26';
% '#ba2121';
% '#93201c';
% '#6e1d18';
% '#4a1813';
% '#29110b';
% '#000000';};
% 
% % gray --> black 
% gray_light_gradients = {'#b8b8b8';
% '#c1c1c1';
% '#c9c9c9';
% '#d2d2d2';
% '#dbdbdb';
% '#e4e4e4';
% '#ededed';
% '#f6f6f6';
% '#ffffff';};
% 
% 
% gray_dark_gradients = {'#b8b8b8';
% '#9f9f9f';
% '#868686';
% '#6f6f6f';
% '#585858';
% '#424242';
% '#2d2d2d';
% '#1a1a1a';
% '#000000';};
% 
% plt_clrs.mgl0 =   plt_clrs.lightgray; gray_light_gradients{7};'#A97C50';
% plt_clrs.mgl6 =  plt_clrs.green; gray_light_gradients{3};
% plt_clrs.mgl8 =  plt_clrs.green;gray_dark_gradients{3};
% plt_clrs.mgl10 =  plt_clrs.green;gray_dark_gradients{7};
% 
% plt_clrs.lightorange = orange_light_gradients{3}; 
% plt_clrs.darkorange  = orange_dark_gradients{2}; 
% plt_clrs.lightgreen  = green_light_gradients{3};  
% plt_clrs.darkgreen   = green_dark_gradients{3}; 
% plt_clrs.lightblue   = blue_light_gradients{4}; 
% plt_clrs.darkblue    = blue_dark_gradients{3}; 
% plt_clrs.lightred    = red_light_gradients{5};
% plt_clrs.darkred     = red_dark_gradients{3}; 
% 
% camp_gr      =  plt_clrs.lightgray;  % gr. vs cAMP in glucose
% wt_camp_gr   =  'k';
% 
% glu_color    = plt_clrs.lightgray; 
% wt_camp_gl   =  'k';
% 
% eth_color    = plt_clrs.lightgreen; 
% wt_camp_eh   =  plt_clrs.darkgreen;
% 
% camp_cell    = plt_clrs.lightorange;
% wt_camp_cell =  plt_clrs.darkorange;
% 
% fru_color     = plt_clrs.lightblue ; 
% wt_camp_fr    = plt_clrs.darkblue; 
% 
% suc_color     = plt_clrs.lightred; 
% wt_camp_su    = plt_clrs.darkred;
% 
% plt_clrs.gl_colors = plt_clrs.gray;
% plt_clrs.eh_colors = plt_clrs.green;
% plt_clrs.cell_colors = plt_clrs.yellow;
% 
% wt_maker_size = 10;10; 
% wt_maker_shape = 'd'; 
% 
% %% ---------------------------- Chemostat -------------------------------------------------------------------------
% %                      Figure 2: Overflow (wild type)
% % ----------------------------------------------------------------------------------------------------------------- 
% % figure 2
% % color_kumar_hg = [102, 194, 165]./255;
% % color_hackett  = plt_clrs.blue; [252, 141,  98]./255;
% % color_kumar_lg = [141, 160, 203]./255;
% % color_boer     = plt_clrs.yellow; 
% % color_xia      = [231, 138, 195]./255;
% 
% %{
% % high glucose, green with different shape
% color_xia        = color_kumar_hg;
% color_johansson  = color_kumar_hg;
% color_postma     = color_kumar_hg;
% color_cortassa   = color_kumar_hg;
% shape_kumar_hg   = 'o'; 
% shape_xia        = '^'; 
% shape_johansson  = 'd'; 
% shape_postma     = 's'; 
% shape_cortassa   = 'v';
% 
% 
% % low glucose, orange with different shape
% color_kumar_lg = color_hackett;
% color_boer     = color_hackett;
% shape_hackett  = 'o'; 
% shape_kumar_lg = 'd'; 
% shape_boer     = 's';
% 
% % batch 
% color_solis = color_murphy; 
% shape_solis ='s';
% %}
% 
% %% ----------------------- Snf1 and Hxt mutant ---------------------------------------------------------------------
% %                     Figure 2: Overflow (mutant)
% % ------------------------------------------------------------------------------------------------------------------ 
% % color_hxt_wt   = [141, 160, 203]./255;
% % color_hxt_hxt1 = [252, 141,  98]./255;
% % color_hxt_hxt7 = [102, 194, 165]./255;
% % color_hxt_tm6  = [231, 138, 195]./255;
% 
% color_hxt_wt   = plt_clrs.blue;
% color_hxt_hxt1 = '#4894AB';
% color_hxt_hxt7 = '#71B0A2';
% color_hxt_tm6  = plt_clrs.green;
% 
% color_snf1_wt     = plt_clrs.blue; 
% color_snf1_mutant = plt_clrs.green;
% 
% color_hxt_bg    = '#C7EAFB';
% alpha_hxt_bg    = 0.3;
% color_snf1_bg   = '#FFF9AE';
% alpha_snf1_bg   = 0.3;
% %% ------------------------------  Batch --------------------------------------------------------------------------
% %                       Figure 3: Diauxic shift
% % ----------------------------------------------------------------------------------------------------------------- 
% 
% % figure 3
% %{
% color_bart   = [102, 194, 165]./255; % green
% color_murphy = [252, 141,  98]./255; % orange
% color_zamp   = [141, 160, 203]./255; % purple
% %} 
% color_bart   = [141, 160, 203]./255; % purple
% color_murphy = plt_clrs.blue;[252, 141,  98]./255; % orange
% color_zamp   = [102, 194, 165]./255; % green 
% 
% %% ------------------------------  diauxic lag  -------------------------------------------------------------------
% %                      Figure 5: diauxic lag time
% % ----------------------------------------------------------------------------------------------------------------- 
% %{
% color_1  = [158,   1,  66]./255;
% color_2  = [213,  62,  79]./255;
% color_3  = [244, 109,  67]./255;
% color_4  = [253, 174,  97]./255;
% color_5  = [254, 224, 139]./255;
% color_6  = [230, 245, 152]./255;
% color_7  = [171, 221, 164]./255;
% color_8  = [102, 194, 165]./255;
% color_9  = [50,  136, 189]./255;
% color_10 = [94,   79, 162]./255;
% 
% colors_map =   [color_1;
% color_3; 
% color_4;
% color_5;
% color_7;
% color_8;
% color_9;
% color_10;];
% %}
% % light green --> darl blue green
% color_1  = [247, 252,  240]./255;
% color_2  = [224, 243,  219]./255;
% color_3  = [204, 235,  197]./255;
% color_4  = [168, 221,  181]./255;
% color_5  = [123, 204, 196]./255;
% color_6  = [78, 197, 211]./255;
% color_7  = [43, 140, 190]./255;
% color_8  = [8, 104, 172]./255;
% color_9  = [8,  64, 129]./255;
% color_10 = [94,   79, 162]./255;
% 
% % light green --> darl blue green
% colors_map = ['#99CC99';
% '#81BB9E';
% '#69AAA4';
% '#509AA9';
% '#3889AF';
% '#2078B4';
% '#2078B4';];
% 
% color_hap_OE = '#99CC99'; plt_clrs.green;
% color_hap_wt = plt_clrs.blue;
% color_hap_mt = '#99CC99';plt_clrs.green;  
% 
% % color_3  = [204, 235,  197]./255;
% % color_4  = [168, 221,  181]./255;
% % color_5  = [123, 204, 196]./255;
% % color_6  = [78, 197, 211]./255;
% % color_7  = [43, 140, 190]./255;
% % color_8  = [8, 104, 172]./255;
% % color_9  = [8,  64, 129]./255;
% % color_10 = [94,   79, 162]./255;
% %% ------------------------------  Rich vs minimal ----------------------------------------------------------------
% %                      Figure 6: Response to external amino acid
% % ----------------------------------------------------------------------------------------------------------------- 
% 
% % figure 6
% color_complex = plt_clrs.green;  % complex 
% color_minimal = plt_clrs.blue;   % minimal 
% 
% 
% % SI figure 6
% pie_Jat = color_complex; %[252, 141,  98]./255;
% pie_Jas = color_minimal; %[141, 160, 203]./255 ;
% 
% color_snf1_1 = [102, 194, 165]./255;
% %color_tor   = [252, 141,  98]./255;
% color_tor = [141, 160, 203]./255;
