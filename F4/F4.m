cd('\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3\F4');

saveFigFol = '\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3\F4';

%% V vs VO - Ai32/PV-Cre mice - VC

lims = [-1 2];
indRec = ~avoidRec_VO_FM_perCell & VO.Ai32Mouse;

curInd = visRespCell_VC & FRch_base_to_V > 0 & indRec;

FR_noOpto = log10(meanFR_noOpto);
FR_opto = log10(meanFR_opto);

FR_noOpto_lim = FR_noOpto;
FR_noOpto_lim(FR_noOpto_lim<lims(1)) = lims(1);
FR_noOpto_lim(FR_noOpto_lim>lims(2)) = lims(2);
FR_opto_lim = FR_opto;
FR_opto_lim(FR_opto_lim<lims(1)) = lims(1);
FR_opto_lim(FR_opto_lim>lims(2)) = lims(2);

figure

axes('Position',[.1 .1 .65 .65])
hold all
scatter(FR_noOpto_lim(curInd & VO.BSunits),FR_opto_lim(curInd & VO.BSunits),...
    100,'filled',"^",'MarkerFaceColor','m',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',1);
scatter(FR_noOpto_lim(curInd & ~VO.BSunits),FR_opto_lim(curInd & ~VO.BSunits),...
    100,'filled','MarkerFaceColor','c',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',1);
set(gca,'XLim',lims,'YLim',lims,'tickdir','out','Color','none');
plot(lims,lims,'k:')

bins = linspace(lims(1),lims(2),15);
binCenters = ( bins(1:end-1) + bins(2:end) ) / 2;

h_lOff_BS = histcounts(FR_noOpto_lim(curInd & VO.BSunits),bins);
h_lOn_BS = histcounts(FR_opto_lim(curInd & VO.BSunits),bins);
h_lOff_NS = histcounts(FR_noOpto_lim(curInd & ~VO.BSunits),bins);
h_lOn_NS = histcounts(FR_opto_lim(curInd & ~VO.BSunits),bins);

h_lOff_BS_norm = h_lOff_BS / max(h_lOff_BS);
h_lOn_BS_norm = h_lOn_BS / max(h_lOn_BS);
h_lOff_NS_norm = h_lOff_NS / max(h_lOff_NS);
h_lOn_NS_norm = h_lOn_NS / max(h_lOn_NS);

axes('Position',[.77 .1 .15 .65])
hold all
plot(h_lOn_BS_norm,binCenters,'m','linewidth',2)
plot(h_lOn_NS_norm,binCenters,'c','linewidth',2)
set(gca,'XLim',[0 1],'YLim',lims,'Color','none','box','off','ytick',[],...
    'tickdir','out')

axes('Position',[.1 .77 .65 .15])
hold all
plot(binCenters,h_lOff_BS_norm,'m','linewidth',2)
plot(binCenters,h_lOff_NS_norm,'c','linewidth',2)
set(gca,'XLim',lims,'YLim',[0 1],'Color','none','box','off','xtick',[],...
    'tickdir','out')

% saveFig(fullfile(saveFigFol,'VC_FR_V_vs_VO'),'-dpdf')

numCellsTemp = sum(curInd);
numRecTemp = length(unique(VO.recNumber(curInd)));
numMouseTemp = length(unique(VO.mouseID(curInd)));

%% Effects of light on spontaneous firing in AC

lims = [-1 2];
indRec = ~avoidRec_VO_FM_perCell & VO.Ai32Mouse;

curInd = fil & VO.probeNum == 1 & indRec;

FR_noOpto = log10(meanFR_noOpto_baseline);
FR_opto = log10(meanFR_opto_baseline);

FR_noOpto_lim = FR_noOpto;
FR_noOpto_lim(FR_noOpto_lim<lims(1)) = lims(1);
FR_noOpto_lim(FR_noOpto_lim>lims(2)) = lims(2);
FR_opto_lim = FR_opto;
FR_opto_lim(FR_opto_lim<lims(1)) = lims(1);
FR_opto_lim(FR_opto_lim>lims(2)) = lims(2);

figure
hold all
scatter(FR_noOpto_lim(curInd & VO.BSunits),FR_opto_lim(curInd & VO.BSunits),...
    100,'filled',"^",'MarkerFaceColor','m',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',.5);
scatter(FR_noOpto_lim(curInd & ~VO.BSunits),FR_opto_lim(curInd & ~VO.BSunits),...
    100,'filled','MarkerFaceColor','c',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',.5);
set(gca,'XLim',lims,'YLim',lims,'tickdir','out','Color','none');
plot(lims,lims,'k:')

signrank(FR_noOpto(curInd),FR_opto(curInd))

% saveFig(fullfile(saveFigFol,'sponFR_AC_lightVC'),'-dpdf')

numCellsTemp = sum(curInd);
numRecTemp = length(unique(VO.recNumber(curInd)));
numMouseTemp = length(unique(VO.mouseID(curInd)));

%% SU example - VC

timeWin = [-2 2.5];
cellInd = 673;
scatterSize = 2;
yLim_PSTH = [0 1.5];

figure('Position',[1993 824.3333 560 420])

axes('Position',[0.15 0.7 0.8 0.2]);
indSpikes = VO.raster_times_V{cellInd} > timeWin(1) & VO.raster_times_V{cellInd} < timeWin(2);
scatter(VO.raster_times_V{cellInd}(indSpikes),VO.raster_trInd_V{cellInd}(indSpikes),scatterSize,'k','filled','MarkerFaceAlpha',0.5);
set(gca,'Color','none','XColor','none','fontsize',16,'ytick',[50 100],...
    'linewidth',1.5,'tickdir','out','XLim',timeWin);
ylabel('# trials')
title(sprintf('Unit #: %s',num2str(cellInd)));

axes('Position',[0.15 0.47 0.8 0.2]);
indSpikes = VO.raster_times_VO{cellInd} > timeWin(1) & VO.raster_times_VO{cellInd} < timeWin(2);
scatter(VO.raster_times_VO{cellInd}(indSpikes),VO.raster_trInd_VO{cellInd}(indSpikes),scatterSize,'b','filled','MarkerFaceAlpha',0.5);
set(gca,'Color','none','XColor','none','fontsize',16,'ytick',[50 100],...
    'linewidth',1.5,'tickdir','out','XLim',timeWin);
ylabel('# trials');

axes('Position',[0.15 0.17 0.8 0.28]);
indTime = VO.binCenters_VO > timeWin(1) & VO.binCenters_VO < timeWin(2);
plot(VO.binCenters_VO(indTime),VO.PSTHout_V(cellInd,indTime),'k','linewidth',2);
hold on
plot(VO.binCenters_VO(indTime),VO.PSTHout_VO(cellInd,indTime),'b','linewidth',2);
set(gca,'Color','none','box','off','tickdir','out','fontsize',16,...
    'linewidth',1.5,'XLim',timeWin,'YLim',yLim_PSTH);
ylabel('Spikes/trial')
xlabel('Time relative to CMN start (sec)')

% saveFig(fullfile(saveFigFol,'SUex_rasterPSTH_VC'),'-dpdf')

%% SU example - AC

timeWin = [-2 2.5];
cellInd = 177;
scatterSize = 10;
yLim_PSTH = [0 1.5];

figure('Position',[1993 824.3333 560 420])

axes('Position',[0.15 0.7 0.8 0.2]);
indSpikes = VO.raster_times_V{cellInd} > timeWin(1) & VO.raster_times_V{cellInd} < timeWin(2);
scatter(VO.raster_times_V{cellInd}(indSpikes),VO.raster_trInd_V{cellInd}(indSpikes),scatterSize,'k','filled','MarkerFaceAlpha',0.5);
set(gca,'Color','none','XColor','none','fontsize',16,'ytick',[50 100],...
    'linewidth',1.5,'tickdir','out','XLim',timeWin);
ylabel('# trials')
title(sprintf('Unit #: %s',num2str(cellInd)));

axes('Position',[0.15 0.47 0.8 0.2]);
indSpikes = VO.raster_times_VO{cellInd} > timeWin(1) & VO.raster_times_VO{cellInd} < timeWin(2);
scatter(VO.raster_times_VO{cellInd}(indSpikes),VO.raster_trInd_VO{cellInd}(indSpikes),scatterSize,'b','filled','MarkerFaceAlpha',0.5);
set(gca,'Color','none','XColor','none','fontsize',16,'ytick',[50 100],...
    'linewidth',1.5,'tickdir','out','XLim',timeWin);
ylabel('# trials');

axes('Position',[0.15 0.17 0.8 0.28]);
indTime = VO.binCenters_VO > timeWin(1) & VO.binCenters_VO < timeWin(2);
plot(VO.binCenters_VO(indTime),VO.PSTHout_V(cellInd,indTime),'k','linewidth',2);
hold on
plot(VO.binCenters_VO(indTime),VO.PSTHout_VO(cellInd,indTime),'b','linewidth',2);
set(gca,'Color','none','box','off','tickdir','out','fontsize',16,...
    'linewidth',1.5,'XLim',timeWin,'YLim',yLim_PSTH);
ylabel('Spikes/trial')
xlabel('Time relative to CMN start (sec)')

% saveFig(fullfile(saveFigFol,'SUex_rasterPSTH_AC'),'-dpdf')

%% V vs VO - Ai32/PV-Cre mice - AC

lims = [-1 4];

indRec = ~avoidRec_VO_FM_perCell & VO.Ai32Mouse;

curInd = visRespCell_AC & FRch_base_to_V > 0 & indRec;

FRch_base_to_V_lims = FRch_base_to_V;
FRch_base_to_V_lims(FRch_base_to_V_lims < lims(1)) = lims(1);
FRch_base_to_V_lims(FRch_base_to_V_lims > lims(2)) = lims(2);
FRch_base_to_V_opto_lims = FRch_base_to_V_opto;
FRch_base_to_V_opto_lims(FRch_base_to_V_opto_lims < lims(1)) = lims(1);
FRch_base_to_V_opto_lims(FRch_base_to_V_opto_lims > lims(2)) = lims(2);

figure
hold all
scatter(FRch_base_to_V_lims(curInd & VO.BSunits),FRch_base_to_V_opto_lims(curInd & VO.BSunits),...
    50,'filled',"^",'MarkerFaceColor','k',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
scatter(FRch_base_to_V_lims(curInd & ~VO.BSunits),FRch_base_to_V_opto_lims(curInd & ~VO.BSunits),...
    50,'filled','MarkerFaceColor','k',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
set(gca,'XLim',lims,'YLim',lims,'tickdir','out');
plot(lims,lims,'k:')

% saveFig(fullfile(saveFigFol,'scatter_V_vs_VO'),'-dpdf')

numCellsTemp = sum(curInd)
numRecTemp = length(unique(VO.recNumber(curInd)))
numMouseTemp = length(unique(VO.mouseID(curInd)))

% Group-level stats
FRchange_noOpto = getFRchangeVals_log2(VO.FRs_baseline,VO.FRs_Vonset,indexToNonOptoCells,[-inf inf]);
FRchange_opto = getFRchangeVals_log2(VO.FRs_baseline,VO.FRs_Vonset,indexToOptoCells,[-inf inf]);
pVal = signrank(FRchange_noOpto(curInd),FRchange_opto(curInd));
meanSD = [ mean(FRchange_noOpto(curInd)) std(FRchange_noOpto(curInd)) ; ...
    mean(FRchange_opto(curInd)) std(FRchange_opto(curInd)) ]
medDiff_optoEffects_FRchange = median(FRchange_opto(curInd) - FRchange_noOpto(curInd));

FRchange_noOpto = getFRdiff(VO.FRs_baseline, VO.FRs_Vonset, indexToNonOptoCells);
FRchange_opto = getFRdiff(VO.FRs_baseline, VO.FRs_Vonset, indexToOptoCells);
medDiff_optoEffects_FRdiff = median(FRchange_opto(curInd) - FRchange_noOpto(curInd));

