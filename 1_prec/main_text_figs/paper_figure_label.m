% 
% labelx      = -0.4;
% labely      = 1.15;
% fontsize    = 13;
% 
% set(0, 'DefaultLineLineWidth',  1);
% %set(0, 'DefaultLineLineWidth',  2.3);
% set(0, 'DefaultLineMarkerSize', 5.5);
% %set(0, 'DefaultLineMarkerEdgeColor', 'k');
% set(0, 'DefaultLineMarkerEdgeColor',"auto")
% set(0, 'DefaultLineMarkerFaceColor',"auto")
% set(0, 'DefaultAxesFontSize', 6)
% set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
% set(0,'DefaultTextColor', 'k') 
% %set(0,'DefaultTextColor', 'k') 
% %% ---------------------------- Chemostat -------------------------------------------------------------------------
% %                      Figure 2: Overflow (wild type)
% % ----------------------------------------------------------------------------------------------------------------- 
% % fig2.figure_output = ...
% % {'glucose', 'precursor', 'ethanol', 'aa_in' ; 
% % 'atp'     , 'lipid'    , 'NH4'    ,  'cells'    ; 
% % 'Egy'     ,  'Efe'     , 'Elp'    ,  'Elo' ,       ;
% % 'Egn'     , 'Emt'      , 'Eas'    ,  'REZ'   ; 
% % 'snf1'    , 'tor'      , 'Eat'    ,  ' '     }'; 
% fig2.figure_output = ...
% {'glucose', 'precursor'  , 'ethanol', 'aa_in' ; 
% 'atp'      , 'cells'      , 'NH4'    , 'JNH4'   ;
% 'snf1'     , 'tor'        , 'Egy'    , 'Efe'   ;
% 'Elp'      , 'Egn'        , 'Emt'    , 'Eas'   ; 
% 'Eat'      , 'Esp'        , 'Esd'    , 'REZ'     ; 
% 'RNA'      , 'protein'    , 'lipid'  ,  'carbo'}'; 
% 
% fig2.figure_id = ...
% {'a'      ,    'b'    ,  'c'      ,  'd'    ; 
% 'e'      ,    'f'    ,  'g'      ,  'h'    ;
% 'i'      ,    'j'    ,  'k'      ,  'l'    ;
% 'a'      ,    'b'    ,  'c'      , ' '     ; 
% 'a'      ,    'b'    ,  'c'      , ' '     ; }';  
% 
% 
% 
% fig2.y_label.glu  = {'glu. uptake (\muM/h)'};
% fig2.y_label.eth  = {'eth. prod. (\muM/h)'};
% fig2.y_label.jnh4 = {'J_{NH4} (\muM/h)'};
% fig2.y_label.cell = {'cell number'}; % 'number of cells'
% fig2.y_label.pc   = {'norm. precursor'};
% fig2.y_label.aa   = {'norm. amino acid_{in}'};
% fig2.y_label.ae   = {'norm ATP'};
% fig2.y_label.sc   = {'storage carb. (/muM)'}; % 'storage carb. (/muM)'
% fig2.y_label.snf1 = {'act. SNF1 fr.'};
% fig2.y_label.tor  = {'act. TORC1 fr.'};
% fig2.y_label.rez  = {'log_{2}(norm. REZ)'};
% fig2.y_label.gy   = {'log_{2}(norm. E_{gy})'};
% fig2.y_label.fe   = {'log_{2}(norm. E_{fe})'};
% fig2.y_label.gn   = {'log_{2}(norm. E_{gn})'};
% fig2.y_label.mt   = {'log_{2}(norm. E_{mt})'};
% fig2.y_label.as   = {'log_{2}(norm. E_{as})'};
% fig2.y_label.sp   = {'log_{2}(norm. E_{sp})'};
% fig2.y_label.sd   = {'log_{2}(norm. E_{sd})'};
% fig2.y_label.at   = {'log_{2}(norm. E_{at})'};
% fig2.y_label.lipid = {'lipid (\muM)'};
% fig2.y_label.nh4   = {'NH4 (\muM)'};
% fig2.y_label.lp    = {'log_{2}(norm. E_{lp})'};
% fig2.y_label.lo    = {'log_{2}(norm. E_{lo})'};
% 
% %% ----------------------- Snf1 and Hxt mutant ---------------------------------------------------------------------
% %                     Figure 2: Overflow (mutant)
% % ------------------------------------------------------------------------------------------------------------------ 
% 
% 
% fig2.subplt.hxt.glu_sim         = 't'; % hxt 
% fig2.y_label.hxt.glu_sim        = {'glucose (\muM)'}; 
% 
% fig2.subplt.hxt.eth_sim         = 't';
% fig2.y_label.hxt.eth_sim        = {'ethanol (\muM)'};
% 
% fig2.subplt.hxt.glu_exp         = 'u';
% fig2.y_label.hxt.glu_exp        = {'glucose (\muM)'};
% 
% fig2.subplt.hxt.eth_exp         = 'u';
% fig2.y_label.hxt.eth_exp        = {'ethanol (\muM)'};
% 
% fig2.subplt.snf1_mt.glu_sim     = 'v'; % snf1 glucose eth flux
% fig2.y_label.snf1_mt.glu_sim    = {'glu uptake';'eth prod.'}; 
% 
% fig2.subplt.snf1_mt.cell_sim    = 'v';
% fig2.y_label.snf1_mt.cell_sim   = {'cell number'};
% 
% fig2.subplt.snf1_mt.glu_exp     = 'w';
% fig2.y_label.snf1_mt.glu_exp    = {'glu uptake';'eth prod.'}; 
% 
% fig2.subplt.snf1_mt.cell_exp    = 'w';
% fig2.y_label.snf1_mt.cell_exp   = {'cell number'};
% 
% 
% y_label_flux_toward_pc = {'prec influx'}; 
% y_label_flux_away_pc   = {'prec outflux'}; 
% 
% y_label_flux_toward_ae = {'ATP influx'}; 
% y_label_flux_away_ae   = {'ATP outflux'};
% 
% %% ------------------------------  Batch --------------------------------------------------------------------------
% %                       Figure 3: Diauxic shift
% % -----------------------------------------------------------------------------------------------------------------
% 
% % 
% % fig3.figure_output = ...
% % {'glucose', 'precursor', 'ethanol', 'aa_in' ; 
% % 'atp'     , 'lipid'    , 'NH4'    ,  'cells'    ; 
% % 'Egy'     ,  'Efe'     , 'Elp'    ,  'Elo' ,       ;
% % 'Egn'     , 'Emt'      , 'Eas'    ,  'REZ'   ; 
% % 'snf1'    , 'tor'      , 'Eat'    ,  ' '     }'; 
% fig3.figure_output = ...
% {'glucose', 'precursor'  , 'ethanol', 'aa_in' ; 
% 'atp'      , 'cells'      , 'NH4'    , 'JNH4'   ;
% 'snf1'     , 'tor'        , 'Egy'    , 'Efe'   ;
% 'Elp'      , 'Egn'        , 'Emt'    , 'Eas'   ; 
% 'Eat'      , 'Esp'        , 'Esd'    , 'REZ'     ; 
% 'RNA'      , 'protein'    , 'lipid'  ,  'carbo'}'; 
% 
% 
% fig3.figure_id = ...
% {'a'      ,    'b'    ,  'c'      ,  'd'    ; 
% 'e'      ,    'f'    ,  'g'      ,  'h'    ;
% 'i'      ,    'j'    ,  'k'      ,  'l'    ;
% 'a'      ,    'b'    ,  'c'      , ' '     }';  
% 
% fig3.y_label.glu  = {'glucose (\muM)'};
% fig3.y_label.eth  = {'ethanol (\muM)'};
% fig3.y_label.jnh4 = {'J_{NH4} (\muM/h)'};
% fig3.y_label.cell = {'cell number'}; % 'number of cells'
% fig3.y_label.pc   = {'norm. precursor'};
% fig3.y_label.aa   = {'norm. amino acid_{in}'};
% fig3.y_label.ae   = {'norm. ATP'};
% fig3.y_label.sc   = {'storage carb. (/muM)'}; % 'storage carb. (/muM)'
% fig3.y_label.snf1 = {'act. SNF1 fr.'};
% fig3.y_label.tor  = {'act. TORC1 fr.'};
% fig3.y_label.rez  = {'protein fractions'};
% fig3.y_label.gy   = {'log_{2}(norm. E_{gy})'};
% fig3.y_label.fe   = {'log_{2}(norm. E_{fe})'};
% fig3.y_label.gn   = {'log_{2}(norm. E_{gn})'};
% fig3.y_label.mt   = {'log_{2}(norm. E_{mt})'};
% fig3.y_label.as   = {'log_{2}(norm. E_{as})'};
% fig3.y_label.sp   = {'log_{2}(norm. E_{sp})'};
% fig3.y_label.sd   = {'log_{2}(norm. E_{sd})'};
% fig3.y_label.at   = {'log_{2}(norm. E_{at})'};
% fig3.y_label.lipid = {'lipid (\muM)'};
% fig3.y_label.nh4   = {'NH4 (\muM)'};
% fig3.y_label.lp    = {'log_{2}(norm. E_{lp})'};
% fig3.y_label.lo    = {'log_{2}(norm. E_{lo})'};
% 
% %% ------------------------------  diauxic lag  -------------------------------------------------------------------
% %                      Figure 5: diauxic lag time
% % ----------------------------------------------------------------------------------------------------------------- 
% fig5.subplt.shift.cell_sim         = 'b';
% fig5.y_label.shift.cell_sim        = {'cells (fc)'};
% 
% fig5.subplt.shift.cell_exp         = 'c';
% fig5.y_label.shift.cell_exp        = {'cells (fc)'};
% 
% fig5.subplt.shift.lag_sim          = 'd';
% fig5.y_label.shift.lag_sim         = {'T_{L} (h)'};
% fig5.x_label.shift.lag_sim         = {'T_{2} (h)'};
% 
% fig5.subplt.shift.lag_exp          = 'e';
% fig5.y_label.shift.lag_exp         = {'T_{l} (h)'};
% 
% fig5.subplt.hap.lag_sim            = 'i';
% fig5.y_label.hap.lag_sim           = {'T_{L}  (h)'};
% 
% fig5.subplt.hap.lag_exp            = 'g';
% fig5.y_label.hap.lag_exp           = {'T_{L} (h)'};
% 
% fig5.subplt.hap.cell_sim           = 'g';
% fig5.y_label.hap.cell_sim          = {'cell number'};
% 
% fig5.subplt.hap.cell_exp           = 'h';
% fig5.y_label.hap.cell_exp          = {'cell number'};
% 
% 
% fig5.subplt.shift.cell_sim         = 'b';
% fig5.y_label.shift.cell_sim        = {'norm.'; 'cell number'};
% 
% fig5.subplt.shift.cell_exp         = 'c';
% fig5.y_label.shift.cell_exp        = {'norm.'; 'cell number'};
% 
% fig5.subplt.shift.lag_sim          = 'd';
% fig5.y_label.shift.lag_sim         = {'T_{L}(h)'};
% fig5.x_label.shift.lag_sim         = {'T_{1} (h)'};
% 
% fig5.subplt.shift.lag_exp          = 'e';
% fig5.y_label.shift.lag_exp         = {'T_{L} (h)'};
% 
% fig5.subplt.hap.lag_sim            = 'i';
% fig5.y_label.hap.lag_sim           = {'T_{L} (h)'};
% 
% fig5.subplt.hap.lag_exp            = 'g';
% fig5.y_label.hap.lag_exp           = {'T_{L} (h)'};
% 
% fig5.subplt.hap.cell_sim           = 'g';
% fig5.y_label.hap.cell_sim          = {'cell number'};
% 
% fig5.subplt.hap.cell_exp           = 'h';
% fig5.y_label.hap.cell_exp          = {'cell number'};
% 
% 
% 
% %% ------------------------------  Rich vs minimal ----------------------------------------------------------------
% %                      Figure 6: Response to external amino acid
% % ----------------------------------------------------------------------------------------------------------------- 
% fig6.subplt.glu_sim  = 'b';
% fig6.y_label.glu_sim  = {'glucose (\muM)'};
% 
% fig6.subplt.glu_exp  = 'c';
% fig6.y_label.glu_exp  = {'glucose (\muM)'};
% 
% fig6.subplt.eth_sim  = 'd';
% fig6.y_label.eth_sim  = {'ethanol (\muM)'};
% 
% fig6.subplt.eth_exp  = 'e';
% fig6.y_label.eth_exp  = {'ethanol (\muM)'};
% 
% fig6.subplt.flux_sim = 'f';
% fig6.y_label.flux_sim  = {'glu uptake,'; 'eth prod .(\muM/h)'};
% fig6.y_label.grate_sim = {'growth rate'; '(h^{-1})'};
% 
% fig6.subplt.flux_exp = 'g';
% fig6.y_label.flux_exp = {'glu uptake,'; 'eth prod .(\muM/h)'};
% fig6.y_label.grate_exp = {'growth rate'; '(h^{-1})'};
% 
% fig6.subplt.prot_sim = 'h';
% fig6.y_label.prot_sim = {'fractions'};
% 
% fig6.subplt.prot_exp = 'i';
% fig6.y_label.prot_exp = {'fractions'};
% 
% 
% 
% %%  -------------------  SI figures ------------
% %% lag SI 
% %
% fig5.subplt.shift.gn_vs_time       = 'a'; 
% fig5.y_label.shift.gn_vs_time      = {'E_{gn} (fc)'};
% 
% fig5.subplt.shift.gn_vs_lag        = 'b';
% fig5.y_label.shift.gn_vs_lag       = {'T_{L} (h)'};
% fig5.x_label.shift.gn_vs_lag       = {'E_{gn} distance'};
% %}
% 
% %
% fig5.subplt.hap.gn                 = 'c';
% fig5.y_label.hap.gn                = {'E_{gn}'};
% 
% fig5.subplt.hap.mt                 = 'd';
% fig5.y_label.hap.mt                = {'E_{mt}'};
% %}
% 
% %% figure 6 SI
% fig6_SI_all.subplt.glu   = 'a';
% fig6_SI_all.y_label.glu  = {'glucose (\muM)'};
% 
% fig6_SI_all.subplt.eth   = 'b';
% fig6_SI_all.y_label.eth  = {'ethanol (\muM)'};
% 
% fig6_SI_all.subplt.cell  = 'c';
% fig6_SI_all.y_label.cell = {'cell number'}; % 'number of cells'
% 
% fig6_SI_all.subplt.pc    = 'd';
% fig6_SI_all.y_label.pc   = {'precursor  (\muM)'};
% 
% fig6_SI_all.subplt.aa    = 'e';
% fig6_SI_all.y_label.aa   = {'amino acid_{in} (\muM)'};
% 
% fig6_SI_all.subplt.ae    = 'f';
% fig6_SI_all.y_label.ae   = {'atp  (\muM)'};
% 
% fig6_SI_all.subplt.sc    = 'g';
% fig6_SI_all.y_label.sc   = {'storage carb.'}; % 'storage carb. (/muM)'
% 
% fig6_SI_all.subplt.snf1  = 'h';
% fig6_SI_all.y_label.snf1 = {'active PKA/Snf1'};
% 
% fig6_SI_all.subplt.tor   = 'i';
% fig6_SI_all.y_label.tor  = {'active TORC1'};
% 
% fig6_SI_all.subplt.rez   = 'j';
% fig6_SI_all.y_label.rez  = {'protein sectors frac.'};
% 
% fig6_SI_all.subplt.gy    = 'k';
% fig6_SI_all.y_label.gy   = {'E_{gy} frac'};
% 
% fig6_SI_all.subplt.fe    = 'l';
% fig6_SI_all.y_label.fe   = {'E_{fe} frac'};
% 
% fig6_SI_all.subplt.gn    = 'm';
% fig6_SI_all.y_label.gn   = {'E_{gn} frac'};
% 
% fig6_SI_all.subplt.mt    = 'n';
% fig6_SI_all.y_label.mt   = {'E_{mt} frac'};
% 
% fig6_SI_all.subplt.as    = 'o';
% fig6_SI_all.y_label.as   = {'E_{as} frac'};
% 
% fig6_SI_all.subplt.sp    = 'p';
% fig6_SI_all.y_label.sp   = {'E_{sp} frac'};
% 
% fig6_SI_all.subplt.sd    = 'q';
% fig6_SI_all.y_label.sd   = {'E_{sd} frac'};
% 
% fig6_SI_all.subplt.at    = 'r';
% fig6_SI_all.y_label.at   = {'E_{at} frac'};
% 
% 
% %% figure X SI glucose pulse
% fig_SI_pulse.subplt.glu   = 'a';
% fig_SI_pulse.y_label.glu  = {'glucose (\muM)'};
% 
% fig_SI_pulse.subplt.eth   = 'b';
% fig_SI_pulse.y_label.eth  = {'ethanol (\muM)'};
% 
% fig_SI_pulse.subplt.cell  = 'c';
% fig_SI_pulse.y_label.cell = {'cell number'}; % 'number of cells'
% 
% 
% fig_SI_pulse.subplt.jgy    = 'd';
% fig_SI_pulse.y_label.jgy   = {'precursor  (fc)'};
% 
% fig_SI_pulse.subplt.jfe    = 'e';
% fig_SI_pulse.y_label.jfe   = {'amino acid_{in} (fc)'};
% 
% fig_SI_pulse.subplt.grate    = 'f';
% fig_SI_pulse.y_label.grate   = {'growth rate (h^{-1})'};
% 
% 
% fig_SI_pulse.subplt.pc    = 'g';
% fig_SI_pulse.y_label.pc   = {'precursor  (fc)'};
% 
% fig_SI_pulse.subplt.aa    = 'h';
% fig_SI_pulse.y_label.aa   = {'amino acid_{in} (fc)'};
% 
% fig_SI_pulse.subplt.ae    = 'i';
% fig_SI_pulse.y_label.ae   = {'atp  (fc)'};
% 
% fig_SI_pulse.subplt.sc    = 'j';
% fig_SI_pulse.y_label.sc   = {'storage carb.'}; % 'storage carb. (/muM)'
% 
% fig_SI_pulse.subplt.snf1  = 'k';
% fig_SI_pulse.y_label.snf1 = {'active PKA/Snf1'};
% 
% fig_SI_pulse.subplt.tor   = 'l';
% fig_SI_pulse.y_label.tor  = {'active TORC1'};
% 
% fig_SI_pulse.subplt.r    = 'm';
% fig_SI_pulse.y_label.r  = {'R log_{2}(fc)'};
% 
% fig_SI_pulse.subplt.z    = 'n';
% fig_SI_pulse.y_label.z   = {'Z log_{2}(fc)'};
% 
% fig_SI_pulse.subplt.sp    = 'o';
% fig_SI_pulse.y_label.sp   = {'E_{sp} log_{2}(fc)'};
% 
% fig_SI_pulse.subplt.gy    = 'p';
% fig_SI_pulse.y_label.gy   = {'E_{gy} log_{2}(fc)'};
% 
% fig_SI_pulse.subplt.fe    = 'q';
% fig_SI_pulse.y_label.fe   = {'E_{fe} log_{2}(fc)'};
% 
% fig_SI_pulse.subplt.sd    = 'r';
% fig_SI_pulse.y_label.sd   = {'E_{sd} log_{2}(fc)'};
% 
% fig_SI_pulse.subplt.gn    = 's';
% fig_SI_pulse.y_label.gn   = {'E_{gn} log_{2}(fc)'};
% 
% fig_SI_pulse.subplt.mt    = 't';
% fig_SI_pulse.y_label.mt   = {'E_{mt} log_{2}(fc)'};
% 
% fig_SI_pulse.subplt.as    = 'u';
% fig_SI_pulse.y_label.as   = {'E_{as} log_{2}(fc)'};
% 
% 
% 
% 
% fig_SI_pulse.subplt.Jgy_short_time     = 's';
% fig_SI_pulse.subplt.Jfe_short_time     = 't';
% fig_SI_pulse.subplt.pc_short_time      = 'u';
% fig_SI_pulse.subplt.snf1_short_time    = 'v';
% fig_SI_pulse.subplt.Emt_short_time     = 'w';
% fig_SI_pulse.subplt.R_short_time       = 'x';
% 
% 

%% ---------------------------- CABBI figures  -------------------------------------------------------------------------
%                      Figure Sc clim vs nlim
% ----------------------------------------------------------------------------------------------------------------- 


% figcabbi.figure_output = ...
%  {'glucose', 'precursor', 'ethanol', 'aa_in' , 'atp'  ; 
%  'cells'   , 'NH4'      , 'protein',  'lipid', 'snf1'  ;
%  'tor'     , 'Elp'      ,  'Egn'   , 'Emt'   , 'Eas',   ;
%  'Egy'     , 'Efe'      , 'Eat'    , 'Elo'   ,  'Esd'   ; 
%  'JNH4'    , 'REZ'      , 'RNA'    ,  'carbo',  'R'}'; 

%figcabbi.figure_output = ...
%  {'glucose', 'ethanol', 'cells'  ,'Egn'    , 'atp'  ; 
%  'protein' , 'lipid'  , 'REZ'    ,'Emt'    , 'snf1' ;
%  'tor'     , 'Elp'    ,'aa_in'   , 'NH4'   , 'Eas', ;
%  'Egy'     , 'Efe'    , 'Eat'    , 'Esp'   , 'Esd'  ; 
%  'JNH4'    , 'precursor' , 'RNA' ,  'carbo',  'R'}' ; 

figcabbi.figure_output = ...
{'glucose', 'ethanol', 'cells'  ,'Egn'    , 'atp'  ; 
'protein' , 'lipid'  , 'REZ'    ,'Emt'    , 'aa_in' ;
'snf1'    ,  'tor'   , 'R'      ,'Elp'   , 'Elo'    ;
'Eas'     , 'Egy'    , 'Efe'    , 'NH4'   , 'Esd'  ; 
'JNH4'    , 'precursor' , 'RNA' ,  'carbo',  'Eat'}' ; 

% clim_lerp
% figcabbi.figure_output = ...
%  {'NH4'    , 'glucose' , 'cells'     ,'protein' , 'snf1'  ; 
%  'aa_in'   , 'ethanol' ,  'precursor', 'lipid', 'tor'  ;
%  'R'       , 'Egy'     ,  'Efe'      , 'REZ'   , 'RNA',   ;
%  'Emt'     , 'Eas'     , 'Egn'       , 'REZ_nlim', 'Eat'   ; 
%  'Elp'     , 'Elo'     , 'Z'        , 'carbo' ,  'atp'}'; 

figcabbi.y_label.glu  = {'glu. uptake (\muM/h)'};
figcabbi.y_label.eth  = {'eth. prod. (\muM/h)'};
figcabbi.y_label.jnh4 = {'J_{NH4} (\muM/h)'};
figcabbi.y_label.cell = {'cells (g/L)'}; % 'number of cells'
figcabbi.y_label.pc   = {'norm. precursor'};
figcabbi.y_label.aa   = {'norm. amino acid_{in}'};
figcabbi.y_label.ae   = {'norm. ATP'};
figcabbi.y_label.sc   = {'storage carb. (/muM)'}; % 'storage carb. (/muM)'
figcabbi.y_label.snf1 = {'act. SNF1 fr.'};
figcabbi.y_label.tor  = {'act. TORC1 fr.'};
figcabbi.y_label.rez  = {'log_{2}(norm. REZ)'};
figcabbi.y_label.gy   = {'log_{2}(norm. E_{gy})'};
figcabbi.y_label.fe   = {'log_{2}(norm. E_{fe})'};
figcabbi.y_label.gn   = {'log_{2}(norm. E_{gn})'};
figcabbi.y_label.mt   = {'log_{2}(norm. E_{mt})'};
figcabbi.y_label.as   = {'log_{2}(norm. E_{as})'};
figcabbi.y_label.sp   = {'log_{2}(norm. E_{sp})'};
figcabbi.y_label.sd   = {'log_{2}(norm. E_{sd})'};
figcabbi.y_label.at   = {'log_{2}(norm. E_{at})'};
figcabbi.y_label.lipid = {'lipid (\muM)'};
figcabbi.y_label.nh4   = {'NH4 (\muM)'};
figcabbi.y_label.lp    = {'log_{2}(norm. E_{lp})'};
figcabbi.y_label.lo    = {'log_{2}(norm. E_{lo})'};
figcabbi.y_label.r     = {'log_{2}(norm. R)'};
figcabbi.y_label.z     = {'log_{2}(norm. Z)'};


fig_batch.figure_output = ...
{'glucose', 'precursor', 'ethanol', 'aa_in' , 'atp'  ,  'cells'   , 'NH4'   ; 
'protein',  'lipid'   , 'snf1'   , 'tor'   , 'R'    ,  'Z'       , 'Egy'   ;
'Efe'   , 'Eas'      ,  'Emt'   , 'Egn'   , 'Elp'  , 'Elo'      ,  ''     ;
'Eat'    , 'JNH4'     ,  'REZ'   , 'RNA'   , 'carbo',  ''       ,  'Esd'   ; }'; 


% abbvie 
fig_batch.figure_output = ...
{'glucose', 'NH4', 'ethanol', 'cells' , 'protein'  ,  'aa_in'   , 'precursor'   ; 
'snf1',  'R'   , 'Emt'   , 'Egn'   , 'lipid'    ,  'tor'       , 'Egy'   ;
'Efe'   , 'Eas'      ,  ''   , ''   , 'Elp'  , 'Elo'      ,  'atp'     ;
'Eat'    , 'JNH4'     ,  'REZ'   , 'RNA'   , 'carbo',  'Z'       ,  'Esd'   ; }'; 

%  {'glucose', 'precursor', 'ethanol', 'aa_in' , 'atp'  ; 
%  'cells'   , 'NH4'      , 'protein',  'lipid', 'snf1'  ;
%  'tor'     , 'Elp'      ,  'Egn'   , 'Emt'   , 'Eas',   ;
%  'Egy'     , 'Efe'      , 'Eat'    , 'Esp'   ,  'Esd'   ; 
%  'JNH4'    , 'REZ'      , 'RNA'    ,  'carbo',  'R'}'; 

% fig_batch.figure_output = ...
%  {'glucose', 'ethanol', 'cells', 'Egn' , 'atp'  ; 
%  'protein' , 'lipid'  , 'REZ',  'Emt', 'snf1'  ;
%  'tor'     , 'Esp'      ,  'aa_in'   , 'NH4'   , 'Eas',   ;
%  'Egy'     , 'Efe'      , 'Eat'    , 'Elp'   ,  'Elo'   ; 
%  'JNH4'    , 'precursor'      , 'RNA'    ,  'carbo',  'R'}'; 


fig_batch.y_label.glu  = {'glucose (\muM)'};
fig_batch.y_label.eth  = {'ethanol (\muM)'};
fig_batch.y_label.jnh4 = {'J_{NH4} (\muM/h)'};
fig_batch.y_label.cell = {'cell (g/L)'}; % 'number of cells'
fig_batch.y_label.pc   = {'norm. precursor'};
fig_batch.y_label.aa   = {'norm. amino acid_{in}'};
fig_batch.y_label.ae   = {'norm. ATP'};
fig_batch.y_label.sc   = {'storage carb. (/muM)'}; % 'storage carb. (/muM)'
fig_batch.y_label.snf1 = {'act. SNF1 fr.'};
fig_batch.y_label.tor  = {'act. TORC1 fr.'};
fig_batch.y_label.rez  = {'protein fractions'};
fig_batch.y_label.gy   = {'log_{2}(norm. E_{gy})'};
fig_batch.y_label.fe   = {'log_{2}(norm. E_{fe})'};
fig_batch.y_label.gn   = {'log_{2}(norm. E_{gn})'};
fig_batch.y_label.mt   = {'log_{2}(norm. E_{mt})'};
fig_batch.y_label.as   = {'log_{2}(norm. E_{as})'};
fig_batch.y_label.sp   = {'log_{2}(norm. E_{sp})'};
fig_batch.y_label.sd   = {'log_{2}(norm. E_{sd})'};
fig_batch.y_label.at   = {'log_{2}(norm. E_{at})'};
fig_batch.y_label.lipid = {'lipid (\muM)'};
fig_batch.y_label.nh4   = {'NH4 (\muM)'};
fig_batch.y_label.lp    = {'log_{2}(norm. E_{lp})'};
fig_batch.y_label.lo    = {'log_{2}(norm. E_{lo})'};
fig_batch.y_label.r     = {'log_{2}(norm. R)'};
fig_batch.y_label.z     = {'log_{2}(norm. Z)'};

