close all

saveFigFol = '\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3\F2_S5_S7\figSave';

cd('\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3\F2_S5_S7');

%% Input variables

timeBounds = [ -2 2.5 ];
exUnit = 28;
exUnit_FMmod = 181;
exUnit_FMnotMod = 493;

miceToUse = true(length(mouseIDforEachRec_VO),1);
% miceToUse( contains(mouseIDforEachRec_VO,'TRO49') ...
%     | contains(mouseIDforEachRec_VO,'TRO51')) = false;
miceToUse_perCell = repelem(miceToUse,numCellPerRec_VO,1);

% Look at primary AC recordings only, or not
primarySecondaryInd = true(length(miceToUse_perCell),1);
primarySecondaryInd = tunedRec_AC_perCell;

indRec = ~avoidRec_VO_FM_perCell & miceToUse_perCell & primarySecondaryInd;

%% Starter variables

tAx = VO.binCenters_VO; % time axis for imagesc
timeInd = tAx > timeBounds(1) & tAx < timeBounds(2);

% Visually responsive cells with increased FRs
indToLook = find(visRespCell_AC & FRch_base_to_V > 0 & ~avoidRec_VO_FM_perCell);

%% Plot - raster and PSTH from example unit

YLIM = [0 2];

figure

subplot(2,1,1)
scatter(VO.raster_times_V{indToLook(exUnit)},VO.raster_trInd_V{indToLook(exUnit)},...
    4,'k','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0)
set(gca,'XLim',timeBounds,'fontsize',14,'linewidth',1,'tickdir','out',...
    'Color','none','box','off')
ylabel('Trial #')

subplot(2,1,2)
plot(tAx,VO.PSTHout_V(indToLook(exUnit),:),'k','linewidth',1.5)
set(gca,'XLim',timeBounds,'fontsize',14,'linewidth',1,'tickdir','out',...
    'Color','none','box','off'); %,'YLim',YLIM)
yFill = [ YLIM(1) YLIM(1) YLIM(2) YLIM(2) YLIM(1) ];
xFill = [ -1 -.5 -.5 -1 -1 ];
hold on
fill(xFill,yFill,'k','FaceAlpha',0.15,'EdgeAlpha',0);
xFill = [ 0 .5 .5 0 0 ];
fill(xFill,yFill,'r','FaceAlpha',0.15,'EdgeAlpha',0);
xlabel('Time relative to CMN start (sec)')
ylabel('Spikes/trial')

drawnow

% saveFig(fullfile(saveFigFol,'SUex_rasterAndPSTH'),'-dpdf');

%% Plot - waveform from example unit

figure

plot(VO.waveforms_tAx*1000,VO.waveforms_mean(indToLook(exUnit),:),'k','linewidth',1.5)
xFill = [ VO.waveforms_tAx fliplr(VO.waveforms_tAx) VO.waveforms_tAx(1) ] * 1000;
yFill = [ VO.waveforms_mean(indToLook(exUnit),:) + VO.waveforms_SD(indToLook(exUnit),:) ...
    fliplr(VO.waveforms_mean(indToLook(exUnit),:) - VO.waveforms_SD(indToLook(exUnit),:)) ...
    VO.waveforms_mean(indToLook(exUnit),1) + VO.waveforms_SD(indToLook(exUnit),1) ];
hold on
fill(xFill,yFill,'k','FaceAlpha',0.15,'EdgeAlpha',0)
set(gca,'fontsize',16,'linewidth',1','tickdir','out','box','off',...
    'Color','none','XLim',[0.8 2.5])

% saveFig(fullfile(saveFigFol,'SUex_waveform'),'-dpdf');

%% Plot - FR change scatter

tickLabels = [ {sprintf('\x2264%s',num2str(2^FRchangeLims(1)))} cellfun(@num2str,num2cell(2.^(FRchangeLims(1)+1:1:FRchangeLims(2)-1)),'un',0) {sprintf('\x2265%s',num2str(2^FRchangeLims(2)))} ];
tickLabels = cellfun(@(x) sprintf('%sx',x),tickLabels,'UniformOutput',false);

figure('Position',[2392 701 493 512])
hold all

% Get scatter point sizes
zScoreScatterSizeRange = [ 30 140 ] ;

lowerSizeBound = min(abs(FRzScore_Vresp( (nonVisRespCell_AC | visRespCell_AC) & indRec) ) );
upperSizeBound = max(abs(FRzScore_Vresp( (nonVisRespCell_AC | visRespCell_AC) & indRec) ) );

scatSizes = getZScoreScatterSizes(abs(FRzScore_Vresp),zScoreScatterSizeRange,[lowerSizeBound upperSizeBound]);
scatSizes_lims = getZScoreScatterSizes([lowerSizeBound upperSizeBound],zScoreScatterSizeRange,[lowerSizeBound upperSizeBound]);
% Plot scatter points - not visually responsive
curInd = nonVisRespCell_AC & VO.BSunits & indRec;
scatter(FRch_base_to_V(curInd),VO.fracDepths(curInd),scatSizes(curInd),'filled',"^",'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',0.35);
curInd = nonVisRespCell_AC & ~VO.BSunits & indRec;
scatter(FRch_base_to_V(curInd),VO.fracDepths(curInd),scatSizes(curInd),'filled','MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',0.35);
% Plot scatter points - visually responsive
curInd = visRespCell_AC & VO.BSunits & indRec;
scatter(FRch_base_to_V(curInd),VO.fracDepths(curInd),scatSizes(curInd),'filled',"^",'MarkerFaceColor','r','MarkerEdgeColor','none','MarkerFaceAlpha',0.6)
curInd = visRespCell_AC & ~VO.BSunits &indRec;
scatter(FRch_base_to_V(curInd),VO.fracDepths(curInd),scatSizes(curInd),'filled','MarkerFaceColor','r','MarkerEdgeColor','none','MarkerFaceAlpha',0.6)

% Plot additional things
hold on
plot([0 0],fracDepthLims,'k:','linewidth',2)
hold on
plot(FRchangeLims,[depthBins(2) depthBins(2)],'k:','linewidth',1)
hold on
plot(FRchangeLims,[depthBins(3) depthBins(3)],'k:','linewidth',1)
hold on
scatter(2.5,0.05,scatSizes_lims(1),'b','filled')
hold on
scatter(3,0.05,scatSizes_lims(2),'b','filled')
% other
set(gca,'xtick',-1:1:3,'xticklabels',tickLabels,'YDir','reverse',...
    'fontsize',16,'linewidth',1,'tickdir','out','YLim',fracDepthLims,...
    'XLim',FRchangeLims,'Color','none','clipping','off')
xlabel({'FR change during V onset';'(V trials)'})
ylabel('Fractional depth');

% Get the number to cells/recs/mice and show in title
indCells = (nonVisRespCell_AC | visRespCell_AC) & indRec;
numCellsTemp = sum(indCells);
numRecTemp = length(unique(VO.recNumber(indCells)));
numMouseTemp = length(unique(VO.mouseID(indCells)));
title(sprintf('# cells": %s, # recs: %s, # mouse: %s',num2str(numCellsTemp),num2str(numRecTemp),num2str(numMouseTemp) ))

drawnow

% saveFig(fullfile(saveFigFol,'scatter_FRchange_V'),'-dpdf');

%% Plot - example FM (recs) and spiking (SU units) imagescs - movement dependent

% % Unit modulated by FM
% curRec = 4;
% curUnitChoice = 'FMmod';
% clims_FM = [-15 15];
% clims_spike = [-200 200];
% figName = 'FMandSpike_FMmod';

% Unit not modulated by FM
curRec = 9;
curUnitChoice = 'FMnotMod';
clims_FM = [-10 10];
clims_spike = [-50 50];
figName = 'FMandSpike_FMnotMod';

if strcmp(curUnitChoice,'FMmod')
    curUnit = exUnit_FMmod;
elseif strcmp(curUnitChoice,'FMnotMod')
    curUnit = exUnit_FMnotMod;
end

% Plot timing
tAx_faceMovements = VO.binCenters_faceMvmts(plotWinInd_VO);
% Spike timing
binSize = 0.1;
bins = plotWin_VO(1):binSize:plotWin_VO(2);
binCenters = (bins(1:end-1) + bins(2:end)) / 2;
binBaseInd = binCenters > baseWin_VO(1) & binCenters < baseWin_VO(2);
% Other
xTickPoints = -1:1:1;

% Index to non-opto trials for current recording
indToRecInCellData = find(VO.recNumber == VCrecIDs_VO(curRec),1,'first');
indNotOpto = find(~VO.FRs_trialType_optoTrials{indToRecInCellData});
% Get mvmt data and z-score
curMvmtData_zScore = movementsStore_VO_CMN_noABS{curRec}(indNotOpto,:);
orderVal = mean(curMvmtData_zScore(:,CMNactWin_VO,1),2);
% Sort movement data by max z-score
[~,orderSort] = sortrows(orderVal);
curDataSort = curMvmtData_zScore(orderSort,plotWinInd_VO);
% Determine which trials have FM above threshold
overThreshold_cur = trialsOverThresh_CMN_VO_noABS{curRec}(indNotOpto,:,1);

numTrials = length(overThreshold_cur);

% Visually responsive units from the current recording
indexToCurVunits = find(VO.recNumber == VCrecIDs_VO(curRec) & visRespCell_AC & FRch_base_to_V > 0);

fig1 = figure('Position',[1490 914 1012 338]);

% Plot imagesc - FM
sub1 = subplot(1,2,1);
hold all
imagesc(curDataSort);
numXsquares = size(curDataSort,2);
colorbar
colormap(sub1,colMap_MG_exp)
clim(clims_FM);

% Get x-axis tick points and limits
xTickPoints_ind = nan(1,length(xTickPoints));
for j = 1:length(xTickPoints)
    xTickPoints_ind(j) = find( abs(tAx_faceMovements-xTickPoints(j)) == min(abs(tAx_faceMovements-xTickPoints(j))),1,'first' );
end

set(gca,'xtick',xTickPoints_ind,'xticklabel',xTickPoints,'fontsize',16,...
    'box','off','Color','none','linewidth',1.5,'tickdir','out',...
    'YLim',[0 numTrials],'XLim',[0 numXsquares],'Clipping','off');

% Plot scatter points to indicate high FM trials
pointsToPlot = find(overThreshold_cur(orderSort));
for j = 1:sum(overThreshold_cur)
    scatter(numXsquares+2 , pointsToPlot(j) ,7,'r','filled');
end
% Plot time zero line
plot([xTickPoints_ind(2) xTickPoints_ind(2)],[1 numTrials],'k:','linewidth',1) ;
% other
xlabel('Time relative to CMN start (sec)')
ylabel('Trial #')

% Get spike times and associated trial IDs
curRastTimes = VO.raster_times_V{curUnit}; % spike times
curRastTrInds = VO.raster_trInd_V{curUnit}; % spike time trials indices\
% Get the binned spike counts for the current unit
binnedSpkCounts = nan(numTrials,length(binCenters));
for j = 1:numTrials
    binnedSpkCounts(j,:) = histcounts(curRastTimes(curRastTrInds==j),bins);
end
binnedSpkCounts_bs = binnedSpkCounts - mean(binnedSpkCounts(:,binBaseInd),2);
binnedSpkCounts_FRs = binnedSpkCounts_bs / binSize;
binnedSpkCounts_FRs_order = binnedSpkCounts_FRs(orderSort,:);

% Plot Imagesc - spike counts
sub2 = subplot(1,2,2);
hold all
imagesc(binnedSpkCounts_FRs_order)
numXsquares = size(binnedSpkCounts_FRs_order,2);
% colorbar
colormap(sub2,colMap_RB_exp)
clim(clims_spike);

xTickPoints_ind = nan(1,length(xTickPoints));
for j = 1:length(xTickPoints)
    xTickPoints_ind(j) = find( abs(binCenters-xTickPoints(j)) == min(abs(binCenters-xTickPoints(j))),1,'first' );
end

set(gca,'xtick',xTickPoints_ind,'xticklabel',xTickPoints,'fontsize',16,...
    'box','off','Color','none','linewidth',1.5,'tickdir','out',...
    'YLim',[0 numTrials],'XLim',[0 numXsquares],'Clipping','off');

% Plot scatter points to indicate high FM trials
for j = 1:sum(overThreshold_cur)
    scatter(numXsquares+1 , pointsToPlot(j) ,7,'r','filled');
end
% Plot time zero line
plot([xTickPoints_ind(2) xTickPoints_ind(2)],[1 numTrials],'k:','linewidth',1) ;
% other
xlabel('Time relative to CMN start (sec)')
ylabel('Trial #')

% saveFig(fullfile(saveFigFol,figName),'-dpdf');

%% Plot waveform for this unit

figure

plot(VO.waveforms_tAx*1000,VO.waveforms_mean(curUnit ,:),'k','linewidth',1.5)
xFill = [ VO.waveforms_tAx fliplr(VO.waveforms_tAx) VO.waveforms_tAx(1) ] * 1000;
yFill = [ VO.waveforms_mean(curUnit,:) + VO.waveforms_SD(curUnit,:) ...
    fliplr(VO.waveforms_mean(curUnit,:) - VO.waveforms_SD(curUnit,:)) ...
    VO.waveforms_mean(curUnit,1) + VO.waveforms_SD(curUnit,1) ];
hold on
fill(xFill,yFill,'k','FaceAlpha',0.15,'EdgeAlpha',0)
set(gca,'fontsize',16,'linewidth',1','tickdir','out','box','off',...
    'Color','none','XLim',[0.8 2.5])

% saveFig(fullfile(saveFigFol,[figName '_waveform']),'-dpdf');

%% Plot - FM scatter - low vs high

curInd = visRespCell_AC & FRch_base_to_V > 0 & indRec;

% Data
curPvals = pVals_lowVsHighMvmt;
indexNoMvmt = index_noOpto_noMvmt;
indexMvmt = index_noOpto_mvmt;
FRchange_noMvmt = FRch_base_to_V_lowMvmt;
FRchange_mvmt = FRch_base_to_V_highMvmt;

FRchange_noMvmt_lims = FRchange_noMvmt;
FRchange_noMvmt_lims(FRchange_noMvmt_lims<FRchangeLims(1)) = FRchangeLims(1);
FRchange_noMvmt_lims(FRchange_noMvmt_lims>FRchangeLims(2)) = FRchangeLims(2);
FRchange_mvmt_lims = FRchange_mvmt;
FRchange_mvmt_lims(FRchange_mvmt_lims<FRchangeLims(1)) = FRchangeLims(1);
FRchange_mvmt_lims(FRchange_mvmt_lims>FRchangeLims(2)) = FRchangeLims(2);

tickLabels = [ {sprintf('\x2264%s',num2str(2^FRchangeLims(1)))} cellfun(@num2str,num2cell(2.^(FRchangeLims(1)+1:1:FRchangeLims(2)-1)),'un',0) {sprintf('\x2265%s',num2str(2^FRchangeLims(2)))} ];
tickLabels = cellfun(@(x) sprintf('%sx',x),tickLabels,'UniformOutput',false);

% Only looking at cells with increased firing in response to the visual stimulus (all V trials)

% Cells with significant differences between low and high movement trials
pVals_lowVsHighMvmt_FDR = curPvals;
pVals_lowVsHighMvmt_FDR(curInd) = mafdr(pVals_lowVsHighMvmt_FDR(curInd),'BHFDR',true);
sigInd = pVals_lowVsHighMvmt_FDR < 0.05;

figure('Position',[1893 647.6667 654 493.3333])
hold all

% Get scatter sizes
scatterSizeRange = [ 65 125 ];
perScatSizes = linspace(10,50,3);
numTrials_lowMvmt = cellfun(@sum,indexNoMvmt);
numTrials_highMvmt = cellfun(@sum,indexMvmt);
ratioOfAll = (numTrials_highMvmt ./ (numTrials_lowMvmt+numTrials_highMvmt))*100;
lowerSizeBound = min( [ min(ratioOfAll(curInd)) perScatSizes(1) ] ) ;
upperSizeBound = max( [ max(ratioOfAll(curInd)) perScatSizes(end) ] ) ;
scatSizes_nonsig_BS = getZScoreScatterSizes(ratioOfAll(curInd & ~sigInd & VO.BSunits),scatterSizeRange,[lowerSizeBound upperSizeBound]);
scatSizes_sig_BS = getZScoreScatterSizes(ratioOfAll(curInd & sigInd & VO.BSunits),scatterSizeRange,[lowerSizeBound upperSizeBound]);
scatSizes_nonsig_NS = getZScoreScatterSizes(ratioOfAll(curInd & ~sigInd & ~VO.BSunits),scatterSizeRange,[lowerSizeBound upperSizeBound]);
scatSizes_sig_NS = getZScoreScatterSizes(ratioOfAll(curInd & sigInd & ~VO.BSunits),scatterSizeRange,[lowerSizeBound upperSizeBound]);
scatSizes_lims = getZScoreScatterSizes( perScatSizes , scatterSizeRange , [ perScatSizes(1) perScatSizes(end) ] );

% Non-significant
scatter(FRchange_noMvmt_lims(curInd & ~sigInd & VO.BSunits),FRchange_mvmt_lims(curInd & ~sigInd & VO.BSunits),...
    scatSizes_nonsig_BS,'filled',"^",'MarkerFaceColor','k',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
scatter(FRchange_noMvmt_lims(curInd & ~sigInd & ~VO.BSunits),FRchange_mvmt_lims(curInd & ~sigInd & ~VO.BSunits),...
    scatSizes_nonsig_NS,'filled','MarkerFaceColor','k',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
% Significant
scatter(FRchange_noMvmt_lims(curInd & sigInd & VO.BSunits),FRchange_mvmt_lims(curInd & sigInd & VO.BSunits),...
    scatSizes_sig_BS,'filled',"^",'MarkerFaceColor','r',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
scatter(FRchange_noMvmt_lims(curInd & sigInd & ~VO.BSunits),FRchange_mvmt_lims(curInd & sigInd & ~VO.BSunits),...
    scatSizes_sig_NS,'filled','MarkerFaceColor','r',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
% Plot the seconary units
scatter(FRchange_noMvmt_lims(curInd & ~tunedRec_AC_perCell),FRchange_mvmt_lims(curInd & ~tunedRec_AC_perCell),...
    20,'filled',"^",'MarkerFaceColor','w',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',1);

% Plot the circle around the movement independent cell
hold on
scatter(FRchange_noMvmt_lims(exUnit_FMnotMod),FRchange_mvmt_lims(exUnit_FMnotMod),...
    scatterSizeRange(2)+scatterSizeRange(2)*0.2,'b');
% Plot the circle around the movement dependent cell
hold on
scatter(FRchange_noMvmt_lims(exUnit_FMmod),FRchange_mvmt_lims(exUnit_FMmod),...
    scatterSizeRange(2)+scatterSizeRange(2)*0.2,'r');
hold on
set(gca,'xtick',FRchangeLims(1):1:FRchangeLims(2),'xticklabels',tickLabels,...
    'fontsize',16,'linewidth',1,'tickdir','out','YLim',FRchangeLims,...
    'XLim',FRchangeLims,'ytick',FRchangeLims(1):1:FRchangeLims(2),'yticklabels',tickLabels,...
    'Color','none','clipping','off');
hold on
plot(FRchangeLims,FRchangeLims,'k:','linewidth',1);
hold on
plot([0 0],FRchangeLims,'k:','linewidth',0.5)
hold on
plot(FRchangeLims,[0 0],'k:','linewidth',0.5)
hold on
scatter(2.4,-0.05,scatSizes_lims(1),'b','filled')
hold on
scatter(2.7,-0.05,scatSizes_lims(2),'b','filled')
hold on
scatter(3,-0.05,scatSizes_lims(3),'b','filled')
xlabel({'FR change during VIS';'(no twitch)'})
ylabel({'FR change during VIS';'(twitch)'});

numCellsTemp = sum(curInd);
numRecTemp = length(unique(VO.recNumber(curInd)));
numMouseTemp = length(unique(VO.mouseID(curInd)));

title(sprintf('# cells": %s, # recs: %s, # mouse: %s',num2str(numCellsTemp),num2str(numRecTemp),num2str(numMouseTemp) ))

daspect([1 1 1])

% Group-level p-value
FRch_lowMvmt_temp = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset,index_noOpto_noMvmt, [-inf inf]);
FRch_highMvmt_temp = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset,index_noOpto_mvmt, [-inf inf]);
signrank(FRch_lowMvmt_temp(curInd),FRch_highMvmt_temp(curInd))

% saveFig(fullfile(saveFigFol,'scatter_lowVsHighMvmt'),'-dpdf');

%% Plot - FM scatter - low vs high - split by recording

% FRchange_noMvmt = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_noMvmt,[-inf inf]);
% FRchange_mvmt = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_mvmt,[-inf inf]);
FRchange_noMvmt = getFRchangeVals(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_noMvmt);
FRchange_mvmt = getFRchangeVals(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_mvmt);

FRchangeLims_split = [-5 5];
yticks = [FRchangeLims_split(1) 0 FRchangeLims_split(2)];
yTicks = { sprintf('%s%s','\leq',num2str(FRchangeLims_split(1))) '0' sprintf('%s%s','\geq',num2str(FRchangeLims_split(2))) };

% Cells to plot
curInd = visRespCell_AC & FRch_base_to_V > 0 & indRec;
% Details for these cells
FRchangeDiff = FRchange_mvmt(curInd) - FRchange_noMvmt(curInd) ;
sigIndCur = sigInd(curInd);
BScur = VO.BSunits(curInd);
curRecs = VO.recNumber(curInd);
curMice = VO.mouseID(curInd);
recInd = sum((curRecs == unique(curRecs)') .* repmat(1:length(unique(curRecs)),sum(curInd),1),2); % x-axis points
xInd = recInd + (rand(sum(curInd),1)/5); % add a little jitter to the x-axis points
ratioOfAllCur = ratioOfAll(curInd); % for getting the scatter point sizes

% Data within limits - FR change
FRchangeDiff_lim = FRchangeDiff;
FRchangeDiff_lim(FRchangeDiff_lim < FRchangeLims_split(1)) = FRchangeLims_split(1);
FRchangeDiff_lim(FRchangeDiff_lim > FRchangeLims_split(2)) = FRchangeLims_split(2);

% Determine the recordings for each mouse
uniMice = unique(curMice);
recsMouse = cell(length(uniMice),1);
for i = 1:length(uniMice)
    recsMouse{i} = unique(recInd(strcmp(curMice,uniMice{i})));
end

figure
hold all

% BS - non-sig
curPlotInd = ~sigIndCur & BScur;
curScatSizes = getZScoreScatterSizes(ratioOfAllCur(curPlotInd),scatterSizeRange,[lowerSizeBound upperSizeBound]);
scatter( xInd(curPlotInd), FRchangeDiff_lim(curPlotInd),...
    curScatSizes,'filled',"^",'MarkerFaceColor','k',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
% NS - non-sig
curPlotInd = ~sigIndCur & ~BScur;
curScatSizes = getZScoreScatterSizes(ratioOfAllCur(curPlotInd),scatterSizeRange,[lowerSizeBound upperSizeBound]);
scatter( xInd(curPlotInd), FRchangeDiff_lim(curPlotInd),...
    curScatSizes,'filled','MarkerFaceColor','k',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
% BS - sig
curPlotInd = sigIndCur & BScur;
curScatSizes = getZScoreScatterSizes(ratioOfAllCur(curPlotInd),scatterSizeRange,[lowerSizeBound upperSizeBound]);
scatter( xInd(curPlotInd), FRchangeDiff_lim(curPlotInd),...
    curScatSizes,'filled',"^",'MarkerFaceColor','r',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
% NS - sig
curPlotInd = sigIndCur & ~BScur;
curScatSizes = getZScoreScatterSizes(ratioOfAllCur(curPlotInd),scatterSizeRange,[lowerSizeBound upperSizeBound]);
scatter( xInd(curPlotInd), FRchangeDiff_lim(curPlotInd),...
    curScatSizes,'filled','MarkerFaceColor','r',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
% Extras
plot([0 length(unique(curRecs))+1],[0 0],'k:','linewidth',1)

scatter(0.5,2,scatSizes_lims(1),'b','filled')
hold on
scatter(1.5,2,scatSizes_lims(2),'b','filled')

set(gca,'XLim',[0 length(unique(curRecs))+1],'YLim',FRchangeLims_split,...
    'xtick',1:length(unique(curRecs)),'tickdir','out','linewidth',1,...
    'Color','none','ytick',yticks,'yticklabel',yTicks);

% Some stats
% p = kruskalwallis(FRchangeDiff,recInd);

numCellsPerRec = sum(recInd == unique(recInd)');
mouseIDforEachRec = repelem(1:length(recsMouse),cellfun(@length,recsMouse));
mouseIDforEachCell = repelem(mouseIDforEachRec,numCellsPerRec);

% p = kruskalwallis(FRchangeDiff,mouseIDforEachCell);

% saveFig(fullfile(saveFigFol,'scatter_lowVsHighMvmt_byRecording'),'-dpdf');

%% Plot - FM scatter - low vs high - suppressed cells (supplementary figure 7)

FRchange_sup = [-1 0];

curInd = visRespCell_AC & FRch_base_to_V < 0 & indRec; % visually suppressed cells

% Data
curPvals = pVals_lowVsHighMvmt;
indexNoMvmt = index_noOpto_noMvmt;
indexMvmt = index_noOpto_mvmt;
FRchange_noMvmt = FRch_base_to_V_lowMvmt;
FRchange_mvmt = FRch_base_to_V_highMvmt;

FRchange_noMvmt_lim = FRchange_noMvmt;
FRchange_noMvmt_lim(FRchange_noMvmt_lim<FRchange_sup(1)) = FRchange_sup(1);
FRchange_noMvmt_lim(FRchange_noMvmt_lim>FRchange_sup(2)) = FRchange_sup(2);
FRchange_mvmt_lim = FRchange_mvmt;
FRchange_mvmt_lim(FRchange_mvmt_lim<FRchange_sup(1)) = FRchange_sup(1);
FRchange_mvmt_lim(FRchange_mvmt_lim>FRchange_sup(2)) = FRchange_sup(2);

tickLabels = [ {sprintf('\x2264%s',num2str(2^FRchange_sup(1)))} cellfun(@num2str,num2cell(2.^(FRchange_sup(1)+1:1:FRchange_sup(2)-1)),'un',0) {sprintf('\x2265%s',num2str(2^FRchange_sup(2)))} ];
tickLabels = cellfun(@(x) sprintf('%sx',x),tickLabels,'UniformOutput',false);

% Cells with significant differences between low and high movement trials
pVals_lowVsHighMvmt_FDR = curPvals;
pVals_lowVsHighMvmt_FDR(curInd) =  mafdr(pVals_lowVsHighMvmt_FDR(curInd),'BHFDR',true);
sigInd = pVals_lowVsHighMvmt_FDR < 0.05;

figure('Position',[1893 647.6667 654 493.3333])
hold all

% Get scatter sizes
scatterSizeRange = [ 65 125 ];
perScatSizes = linspace(10,50,3);
numTrials_lowMvmt = cellfun(@sum,indexNoMvmt);
numTrials_highMvmt = cellfun(@sum,indexMvmt);
ratioOfAll = (numTrials_highMvmt ./ (numTrials_lowMvmt+numTrials_highMvmt))*100;
lowerSizeBound = min( [ min(ratioOfAll(curInd)) perScatSizes(1) ] ) ;
upperSizeBound = max( [ max(ratioOfAll(curInd)) perScatSizes(end) ] ) ;
scatSizes_nonsig_BS = getZScoreScatterSizes(ratioOfAll(curInd & ~sigInd & VO.BSunits),scatterSizeRange,[lowerSizeBound upperSizeBound]);
scatSizes_sig_BS = getZScoreScatterSizes(ratioOfAll(curInd & sigInd & VO.BSunits),scatterSizeRange,[lowerSizeBound upperSizeBound]);
scatSizes_nonsig_NS = getZScoreScatterSizes(ratioOfAll(curInd & ~sigInd & ~VO.BSunits),scatterSizeRange,[lowerSizeBound upperSizeBound]);
scatSizes_sig_NS = getZScoreScatterSizes(ratioOfAll(curInd & sigInd & ~VO.BSunits),scatterSizeRange,[lowerSizeBound upperSizeBound]);
scatSizes_lims = getZScoreScatterSizes( perScatSizes , scatterSizeRange , [ perScatSizes(1) perScatSizes(end) ] );

% Non-significant
scatter(FRchange_noMvmt_lim(curInd & ~sigInd & VO.BSunits),FRchange_mvmt_lim(curInd & ~sigInd & VO.BSunits),...
    scatSizes_nonsig_BS,'filled',"^",'MarkerFaceColor','k',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
scatter(FRchange_noMvmt_lim(curInd & ~sigInd & ~VO.BSunits),FRchange_mvmt_lim(curInd & ~sigInd & ~VO.BSunits),...
    scatSizes_nonsig_NS,'filled','MarkerFaceColor','k',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
% Significant
scatter(FRchange_noMvmt_lim(curInd & sigInd & VO.BSunits),FRchange_mvmt_lim(curInd & sigInd & VO.BSunits),...
    scatSizes_sig_BS,'filled',"^",'MarkerFaceColor','r',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
scatter(FRchange_noMvmt_lim(curInd & sigInd & ~VO.BSunits),FRchange_mvmt_lim(curInd & sigInd & ~VO.BSunits),...
    scatSizes_sig_NS,'filled','MarkerFaceColor','r',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);

set(gca,'xtick',FRchange_sup(1):1:FRchange_sup(2),'xticklabels',tickLabels,...
    'fontsize',16,'linewidth',1,'tickdir','out','YLim',FRchange_sup,...
    'XLim',FRchange_sup,'ytick',FRchange_sup(1):1:FRchange_sup(2),'yticklabels',tickLabels,...
    'Color','none','clipping','off');

plot(FRchange_sup,FRchange_sup,'k:','linewidth',1);
scatter(-.1,-0.05,scatSizes_lims(1),'b','filled')
hold on
scatter(0,-0.05,scatSizes_lims(2),'b','filled')
hold on
scatter(1,-0.05,scatSizes_lims(3),'b','filled')
xlabel({'FR change during VIS';'(no FM)'})
ylabel({'FR change during VIS';'(FM)'});

numCellsTemp = sum(curInd);
numRecTemp = length(unique(VO.recNumber(curInd)));
numMouseTemp = length(unique(VO.mouseID(curInd)));

title(sprintf('# cells": %s, # recs: %s, # mouse: %s',num2str(numCellsTemp),num2str(numRecTemp),num2str(numMouseTemp) ))

daspect([1 1 1])

% Group-level stat
FRch_lowMvmt_temp = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_noMvmt,[-inf inf]);
FRch_highMvmt_temp = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_mvmt,[-inf inf]);
signrank(FRch_lowMvmt_temp(curInd),FRch_highMvmt_temp(curInd))

% saveFig(fullfile(saveFigFol,'scatter_lowVsHighMvmt_sup'),'-dpdf');

