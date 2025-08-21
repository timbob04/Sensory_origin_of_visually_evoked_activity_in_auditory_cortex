close all; clc

figSaveFol = '\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3\F3\figSave';
cd('\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3\F3')

%% Input variables

miceToUse = true(length(mouseIDforEachRec_VO),1);
% miceToUse( contains(mouseIDforEachRec_VO,'TRO49') ...
%     | contains(mouseIDforEachRec_VO,'TRO51')) = false;
miceToUse_perCell = repelem(miceToUse,numCellPerRec_VO,1);

% To look at primary AC units only, or not
primarySecondaryInd = true(length(miceToUse_perCell),1);
% primarySecondaryInd = tunedRec_AC_perCell;

binSizeSpk = mean(diff(bins_spkCounts));

%% Index to noise burst responsive cells in VC and vis responsive cells in AC

curInd_off = noiseRespCell_VC_CMNdur & ~avoidRec_VO_FM_perCell & miceToUse_perCell & FRzScore_noiseResp_CMNlen > 0;
curInd_V = visRespCell_AC & ~avoidRec_VO_FM_perCell & miceToUse_perCell & FRzScore_Vresp > 0;

%% Sound responsive cell in VC example

% Inputs
indCell = 613;
clims_FM = [-10 10];
clims_spike = [-50 50];
plotLims = [-.5 1.5];
xTickPoints = -.5:.5:1.5;
actWindow = [0 .5];
meanLims_FM = [-5 50];
meanLims_sp = [-1 6];

% Plot timing
plotInd = VO.binCenters_faceMvmts > plotLims(1) & VO.binCenters_faceMvmts < plotLims(2);
actWinInd = VO.binCenters_faceMvmts > actWindow(1) & VO.binCenters_faceMvmts < actWindow(2);
tAx_faceMovements = VO.binCenters_faceMvmts(plotInd);
binCenterInd = binCenters_spkCounts > plotLims(1) & binCenters_spkCounts < plotLims(2);
binCenterActInd = binCenters_spkCounts > actWindow(1) & binCenters_spkCounts < actWindow(2);
spkBinTax = binCenters_spkCounts(binCenterInd);

% Trials above threshold
overThresh = trialsOverThresh_off_VO_perCell{indCell}(:,:,1);

figure('Position',[447 818 1146 420])

% Plot the probe plot
axes('Position',[.04 .2 .08 .53]);
curRec = VO.recNumber(indCell);
indCurRec = VCrecIDs_VO == curRec;
indToCurRec = find(VO.recNumber == curRec & VO.probeNum == 2);
xJitter = rand(length(indToCurRec),1);
scatter(xJitter,VO.fracDepths(indToCurRec),65,'k','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0)
set(gca,'XLim',[-1 2],'YDir','reverse','Color','none','XColor','none','tickdir','out',...
    'YLim',[0 1],'ytick',[0 1],'linewidth',1.5,'fontsize',15)
hold on
indToCurUnit = indToCurRec == indCell;
curDepths = VO.fracDepths(indToCurRec);
scatter(xJitter(indToCurUnit),curDepths(indToCurUnit),150,'r','LineWidth',1.5)
title({'Probe';'plot'})

% Plot the face movements
curMvmtData_zScore = movementsStore_VO_offset_noABS{indCurRec};
curzScores = mean(curMvmtData_zScore(:,actWinInd,1),2);
[~,orderSort] = sortrows(curzScores);
curDataSort = curMvmtData_zScore(orderSort,plotInd);
[numTrials,numXsquares] = size(curDataSort);
axFM = axes('Position',[.25 .2 .25 .53]);
hold all
imagesc(curDataSort);
colorbar;
colormap(axFM,colMap_MG_exp)
clim(clims_FM);
% Get x-axis tick points and limits
xTickPoints_ind = nan(1,length(xTickPoints));
for j = 1:length(xTickPoints)
    xTickPoints_ind(j) = find( abs(tAx_faceMovements-xTickPoints(j)) == min(abs(tAx_faceMovements-xTickPoints(j))),1,'first' );
end
plot([xTickPoints_ind(2) xTickPoints_ind(2)],[0 numTrials],'k:','linewidth',1.5)
set(gca,'xtick',xTickPoints_ind,'xticklabel',xTickPoints,'fontsize',16,...
    'box','off','Color','none','linewidth',1.5,'tickdir','out',...
    'YLim',[0 numTrials],'XLim',[0 numXsquares],'Clipping','off');
xlabel('Time relative to noise train start (sec)')
ylabel('Trial #')
% Trials above threshold line
overThresh_sort = find(overThresh(orderSort));
for i = 1:length(overThresh_sort)
    plot([numXsquares numXsquares],[overThresh_sort(i)-.5 overThresh_sort(i)+.5],'g','LineWidth',4)
end
% Means
pause(0.1)
width = axFM.Position(3);
meanAx_FM = axes('Position',[.25 .75 width .07]);
hold all
plot(tAx_faceMovements,mean(curMvmtData_zScore(~overThresh,plotInd)),'k','linewidth',1.5);
plot(tAx_faceMovements,mean(curMvmtData_zScore(overThresh,plotInd)),'g','linewidth',1.5);
plot([tAx_faceMovements(1) tAx_faceMovements(end)],[0 0],'k:')
plot([tAx_faceMovements(1) tAx_faceMovements(1)],[0 50],'k')
set(meanAx_FM,'Color','none','XColor','none','YColor','none','YLim',meanLims_FM)
% Diff
diffAx_FM = axes('Position',[.25 .83 width .07]);
hold all
diffFM = abs(mean(curMvmtData_zScore(overThresh,plotInd))-mean(curMvmtData_zScore(~overThresh,plotInd)));
plot(tAx_faceMovements,diffFM,'k','linewidth',1.5);
plot([tAx_faceMovements(1) tAx_faceMovements(end)],[0 0],'k:')
set(diffAx_FM,'Color','none','XColor','none','YColor','none','YLim',meanLims_FM)

% Spike count plot
axSpike = axes('Position',[.6 .2 .25 .53]);
binnedSpkCounts_bs_order = binnedSpkCounts_bs_store_off{indCell}(orderSort,binCenterInd,1) / binSizeSpk; % in FR
hold all
imagesc(binnedSpkCounts_bs_order)
numXsquares = size(binnedSpkCounts_bs_order,2);
colorbar
colormap(axSpike,colMap_RB_exp)
clim(clims_spike);
xTickPoints_ind = nan(1,length(xTickPoints));
for j = 1:length(xTickPoints)
    xTickPoints_ind(j) = find( abs(spkBinTax-xTickPoints(j)) == min(abs(spkBinTax-xTickPoints(j))),1,'first' );
end
% Indicate over-threshold trials
for i = 1:length(overThresh_sort)
    plot([numXsquares numXsquares],[overThresh_sort(i)-.5 overThresh_sort(i)+.5],'g','LineWidth',4)
end
plot([xTickPoints_ind(2) xTickPoints_ind(2)],[0 numTrials],'k:','linewidth',1.5)
set(gca,'xtick',xTickPoints_ind,'xticklabel',xTickPoints,'fontsize',16,...
    'box','off','Color','none','linewidth',1.5,'tickdir','out',...
    'YLim',[0 numTrials],'XLim',[0 numXsquares],'Clipping','off');
xlabel('Time relative to noise train start (sec)')
ylabel('Trial #')
% Means
pause(0.1)
width = axSpike.Position(3);
meanAx_sp = axes('Position',[.6 .74 width .07]);
plot(binCenters_spkCounts(binCenterInd),binnedSpkCounts_bs_noMvmt_off(indCell,binCenterInd),'k','linewidth',1.5);
hold all
plot(binCenters_spkCounts(binCenterInd),binnedSpkCounts_bs_mvmt_off(indCell,binCenterInd),'g','linewidth',1.5);
plot(plotLims,[0 0],'k:')
plot([plotLims(1) plotLims(1)],[0 5],'k')
set(meanAx_sp,'Color','none','XColor','none','YColor','none','YLim',meanLims_sp)
% Diff
diffAx_spk = axes('Position',[.6 .83 width .07]);
hold all
diffSpk = abs(binnedSpkCounts_bs_mvmt_off(indCell,binCenterInd)-binnedSpkCounts_bs_noMvmt_off(indCell,binCenterInd));
plot(binCenters_spkCounts(binCenterInd),diffSpk,'k','linewidth',1.5);
plot(plotLims,[0 0],'k:')
set(diffAx_spk,'Color','none','XColor','none','YColor','none','YLim',meanLims_sp)

% Plot waveform for this unit
axes('Position',[.9 .2 .08 .2]);
plot(VO.waveforms_tAx*1000,VO.waveforms_mean(indCell ,:),'k','linewidth',1.5)
xFill = [ VO.waveforms_tAx fliplr(VO.waveforms_tAx) VO.waveforms_tAx(1) ] * 1000;
yFill = [ VO.waveforms_mean(indCell,:) + VO.waveforms_SD(indCell,:) ...
    fliplr(VO.waveforms_mean(indCell,:) - VO.waveforms_SD(indCell,:)) ...
    VO.waveforms_mean(indCell,1) + VO.waveforms_SD(indCell,1) ];
hold all
fill(xFill,yFill,'k','FaceAlpha',0.15,'EdgeAlpha',0)
set(gca,'fontsize',16,'linewidth',1','tickdir','out','box','off',...
    'Color','none','XLim',[0.8 2.5],'YColor','none','XColor','none')
minPoint = min(yFill) + min(yFill)*0.1;
minPointText = min(yFill) + ((min(yFill)*0.1)*2);
plot([1 1.5],[minPoint minPoint],'k','LineWidth',1)
text(1.25,minPointText,'0.5 ms','HorizontalAlignment','center','VerticalAlignment','top','FontSize',14)

% % Save
% savename = fullfile(figSaveFol,'example_sound_VC');
% saveFig(savename,'-dpdf')

%% Visually responsive cell in AC example

% Inputs
indCell = 214;
clims_FM = [-5 5];
clims_spike = [-70 70];
plotLims = [-.5 1.5];
xTickPoints = -.5:.5:1.5;
actWindow = [0 .5];
meanLims_FM = [-1 25];
meanLims_sp = [-1 5];

% Plot timing
plotInd = VO.binCenters_faceMvmts > plotLims(1) & VO.binCenters_faceMvmts < plotLims(2);
actWinInd = VO.binCenters_faceMvmts > actWindow(1) & VO.binCenters_faceMvmts < actWindow(2);
tAx_faceMovements = VO.binCenters_faceMvmts(plotInd);
binCenterInd = binCenters_spkCounts > plotLims(1) & binCenters_spkCounts < plotLims(2);
spkBinTax = binCenters_spkCounts(binCenterInd);

noOpto = ~VO.FRs_trialType_optoTrials{indCell};

% Trials above threshold
overThresh = trialsOverThresh_CMN_VO_perCell{indCell}(noOpto,:,1);

figure('Position',[447 818 1146 420])

% Plot the probe plot
axPp = axes('Position',[.04 .2 .08 .53]);
curRec = VO.recNumber(indCell);
indCurRec = find(VCrecIDs_VO == curRec);
indToCurRec = find(VO.recNumber == curRec & VO.probeNum == 1);
xJitter = rand(length(indToCurRec),1);
scatter(xJitter,VO.fracDepths(indToCurRec),65,'k','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0)
set(gca,'XLim',[-1 2],'YDir','reverse','Color','none','XColor','none','tickdir','out',...
    'YLim',[0 1],'ytick',[0 1],'linewidth',1.5,'fontsize',15)
hold on
indToCurUnit = indToCurRec == indCell;
curDepths = VO.fracDepths(indToCurRec);
scatter(xJitter(indToCurUnit),curDepths(indToCurUnit),150,'r','LineWidth',1.5)
title({'Probe';'plot'})

% Plot the face movements
curMvmtData_zScore = movementsStore_VO_CMN_noABS{indCurRec}(noOpto,:,1);
curzScores = mean(curMvmtData_zScore(:,actWinInd,1),2);
[~,orderSort] = sortrows(curzScores);
curDataSort = curMvmtData_zScore(orderSort,plotInd);
[numTrials,numXsquares] = size(curDataSort);
axFM = axes('Position',[.25 .2 .25 .53]);
hold all
imagesc(curDataSort);
colorbar;
colormap(axFM,colMap_MG_exp)
clim(clims_FM);
% Get x-axis tick points and limits
xTickPoints_ind = nan(1,length(xTickPoints));
for j = 1:length(xTickPoints)
    xTickPoints_ind(j) = find( abs(tAx_faceMovements-xTickPoints(j)) == min(abs(tAx_faceMovements-xTickPoints(j))),1,'first' );
end
plot([xTickPoints_ind(2) xTickPoints_ind(2)],[0 numTrials],'k:','linewidth',1.5)
set(gca,'xtick',xTickPoints_ind,'xticklabel',xTickPoints,'fontsize',16,...
    'box','off','Color','none','linewidth',1.5,'tickdir','out',...
    'YLim',[0 numTrials],'XLim',[0 numXsquares],'Clipping','off');
xlabel('Time relative to CMN start (sec)')
ylabel('Trial #')
% Trials above threshold line
overThresh_sort = find(overThresh(orderSort));
for i = 1:length(overThresh_sort)
    plot([numXsquares numXsquares],[overThresh_sort(i)-.5 overThresh_sort(i)+.5],'g','LineWidth',4)
end
% Means
pause(0.1)
width = axFM.Position(3);
meanAx_FM = axes('Position',[.25 .75 width .07]);
hold all
plot(tAx_faceMovements,mean(curMvmtData_zScore(~overThresh,plotInd)),'k','linewidth',1.5);
plot(tAx_faceMovements,mean(curMvmtData_zScore(overThresh,plotInd)),'g','linewidth',1.5);
plot([tAx_faceMovements(1) tAx_faceMovements(end)],[0 0],'k:')
plot([tAx_faceMovements(1) tAx_faceMovements(1)],[0 40],'k')
set(meanAx_FM,'Color','none','XColor','none','YColor','none','YLim',meanLims_FM)
% Diff
diffAx_FM = axes('Position',[.25 .83 width .07]);
hold all
diffFM = abs(mean(curMvmtData_zScore(overThresh,plotInd))-mean(curMvmtData_zScore(~overThresh,plotInd)));
plot(tAx_faceMovements,diffFM,'k','linewidth',1.5);
plot([tAx_faceMovements(1) tAx_faceMovements(end)],[0 0],'k:')
set(diffAx_FM,'Color','none','XColor','none','YColor','none','YLim',meanLims_FM)

% Spike count plot
axSpike = axes('Position',[.6 .2 .25 .53]);
binnedSpkCounts_bs_order = binnedSpkCounts_bs_store_V{indCell}(orderSort,binCenterInd,1) / binSizeSpk; % in FR
hold all
imagesc(binnedSpkCounts_bs_order)
numXsquares = size(binnedSpkCounts_bs_order,2);
colorbar
colormap(axSpike,colMap_RB_exp)
clim(clims_spike);
xTickPoints_ind = nan(1,length(xTickPoints));
for j = 1:length(xTickPoints)
    xTickPoints_ind(j) = find( abs(spkBinTax-xTickPoints(j)) == min(abs(spkBinTax-xTickPoints(j))),1,'first' );
end
% Indicate over-threshold trials
for i = 1:length(overThresh_sort)
    plot([numXsquares numXsquares],[overThresh_sort(i)-.5 overThresh_sort(i)+.5],'g','LineWidth',4)
end
plot([xTickPoints_ind(2) xTickPoints_ind(2)],[0 numTrials],'k:','linewidth',1.5)
set(gca,'xtick',xTickPoints_ind,'xticklabel',xTickPoints,'fontsize',16,...
    'box','off','Color','none','linewidth',1.5,'tickdir','out',...
    'YLim',[0 numTrials],'XLim',[0 numXsquares],'Clipping','off');
xlabel('Time relative to CMN start (sec)')
ylabel('Trial #')
% Means
pause(0.1)
width = axSpike.Position(3);
meanAx_sp = axes('Position',[.6 .74 width .07]);
plot(binCenters_spkCounts(binCenterInd),binnedSpkCounts_bs_noMvmt_V(indCell,binCenterInd),'k','linewidth',1.5);
hold all
plot(binCenters_spkCounts(binCenterInd),binnedSpkCounts_bs_mvmt_V(indCell,binCenterInd),'g','linewidth',1.5);
plot(plotLims,[0 0],'k:')
plot([plotLims(1) plotLims(1)],[0 5],'k')
set(meanAx_sp,'Color','none','XColor','none','YColor','none','YLim',meanLims_sp)
% Diff
diffAx_spk = axes('Position',[.6 .83 width .07]);
hold all
diffSpk = abs(binnedSpkCounts_bs_mvmt_V(indCell,binCenterInd)-binnedSpkCounts_bs_noMvmt_V(indCell,binCenterInd));
plot(binCenters_spkCounts(binCenterInd),diffSpk,'k','linewidth',1.5);
plot(plotLims,[0 0],'k:')
set(diffAx_spk,'Color','none','XColor','none','YColor','none','YLim',meanLims_sp)

% Plot waveform for this unit
axes('Position',[.9 .2 .08 .2]);
plot(VO.waveforms_tAx*1000,VO.waveforms_mean(indCell ,:),'k','linewidth',1.5)
xFill = [ VO.waveforms_tAx fliplr(VO.waveforms_tAx) VO.waveforms_tAx(1) ] * 1000;
yFill = [ VO.waveforms_mean(indCell,:) + VO.waveforms_SD(indCell,:) ...
    fliplr(VO.waveforms_mean(indCell,:) - VO.waveforms_SD(indCell,:)) ...
    VO.waveforms_mean(indCell,1) + VO.waveforms_SD(indCell,1) ];
hold all
fill(xFill,yFill,'k','FaceAlpha',0.15,'EdgeAlpha',0)
set(gca,'fontsize',16,'linewidth',1','tickdir','out','box','off',...
    'Color','none','XLim',[0.8 2.5],'YColor','none','XColor','none')
minPoint = min(yFill) + min(yFill)*0.1;
minPointText = min(yFill) + ((min(yFill)*0.1)*2);
plot([1 1.5],[minPoint minPoint],'k','LineWidth',1)
text(1.25,minPointText,'0.5 ms','HorizontalAlignment','center','VerticalAlignment','top','FontSize',14)

% % Save
% savename = fullfile(figSaveFol,'example_visual_AC');
% saveFig(savename,'-dpdf')

%% Difference plots

clims_spike = [-5 5];
meanLims = [-.2 2];

load('PuOr.mat')

diffs_V = binnedSpkCounts_bs_mvmt_V - binnedSpkCounts_bs_noMvmt_V;
diffs_off = binnedSpkCounts_bs_mvmt_off  -binnedSpkCounts_bs_noMvmt_off;

% Get differences and sort them
[~,sortOrder_V] = sort(delta_V_noABS(curInd_V),'ascend');
dataInd = diffs_V(curInd_V,:);
diffs_V_ind_order = dataInd(sortOrder_V,:);
[~,sortOrder_off] = sort(delta_off_noABS(curInd_off),'ascend');
dataInd = diffs_off(curInd_off,:);
diffs_off_ind_order = dataInd(sortOrder_off,:);

figure('Position',[447 818 1146 420])

% Visual
numRec = size(diffs_V_ind_order,1);
axSpike_V = axes('Position',[.25 .2 .25 .53]);
hold all
datToPlot = diffs_V_ind_order(:,binCenterInd);
imagesc(datToPlot)
numXsquares = size(datToPlot,2);
colorbar
colormap(axSpike_V,PuOr/255)
clim(clims_spike);
xTickPoints_ind = nan(1,length(xTickPoints));
for j = 1:length(xTickPoints)
    xTickPoints_ind(j) = find( abs(spkBinTax-xTickPoints(j)) == min(abs(spkBinTax-xTickPoints(j))),1,'first' );
end
plot([xTickPoints_ind(2) xTickPoints_ind(2)],[0 numRec],'k:','linewidth',1.5)
set(gca,'xtick',xTickPoints_ind,'xticklabel',xTickPoints,'fontsize',16,...
    'box','off','Color','none','linewidth',1.5,'tickdir','out',...
    'YLim',[0 numRec],'XLim',[0 numXsquares],'Clipping','off');
xlabel('Time relative to CMN start (sec)')
ylabel('Trial #')
% Means
pause(0.1)
width = axSpike_V.Position(3);
meanAx_sp = axes('Position',[.25 .76 width .07]);
hold all
plot(binCenters_spkCounts(binCenterInd),mean(datToPlot,'omitnan'),'k','linewidth',1.5);
plot(plotLims,[0 0],'k:')
plot([plotLims(1) plotLims(1)],meanLims,'k')
set(meanAx_sp,'Color','none','XColor','none','YColor','none','YLim',meanLims)

% Noise burst
numRec = size(diffs_off_ind_order,1);
axSpike_off = axes('Position',[.6 .2 .25 .53]);
hold all
datToPlot = diffs_off_ind_order(:,binCenterInd);
imagesc(datToPlot)
numXsquares = size(datToPlot,2);
colorbar
colormap(axSpike_off,PuOr/255)
clim(clims_spike);
xTickPoints_ind = nan(1,length(xTickPoints));
for j = 1:length(xTickPoints)
    xTickPoints_ind(j) = find( abs(spkBinTax-xTickPoints(j)) == min(abs(spkBinTax-xTickPoints(j))),1,'first' );
end
plot([xTickPoints_ind(2) xTickPoints_ind(2)],[0 numRec],'k:','linewidth',1.5)
set(gca,'xtick',xTickPoints_ind,'xticklabel',xTickPoints,'fontsize',16,...
    'box','off','Color','none','linewidth',1.5,'tickdir','out',...
    'YLim',[0 numRec],'XLim',[0 numXsquares],'Clipping','off');
xlabel('Time relative to CMN start (sec)')
ylabel('Trial #')
% Means
pause(0.1)
width = axSpike_off.Position(3);
meanAx_sp = axes('Position',[.6 .76 width .07]);
hold all
plot(binCenters_spkCounts(binCenterInd),mean(datToPlot,'omitnan'),'k','linewidth',1.5);
plot(plotLims,[0 0],'k:')
plot([plotLims(1) plotLims(1)],meanLims,'k')
set(meanAx_sp,'Color','none','XColor','none','YColor','none','YLim',meanLims)

% % Save
% savename = fullfile(figSaveFol,'differencePlots');
% saveFig(savename,'-dpdf')

%% Plot - FM scatter - low vs high - AC and VC

plotLims = [-1 3];

FRch_off_lowMv = FRch_base_to_off_lowMvmt_CMNdur(curInd_off);
FRch_off_highMv = FRch_base_to_off_highMvmt_CMNdur(curInd_off);
FRch_V_lowMv = FRch_base_to_V_lowMvmt(curInd_V);
FRch_V_highMv = FRch_base_to_V_highMvmt(curInd_V);
FRch_lowMv = [ FRch_V_lowMv ; FRch_off_lowMv ];
FRch_highMv = [ FRch_V_highMv ; FRch_off_highMv ];

tickLabels = [ {sprintf('\x2264%s',num2str(2^plotLims(1)))} cellfun(@num2str,num2cell(2.^(plotLims(1)+1:1:FRchangeLims(2)-1)),'un',0) {sprintf('\x2265%s',num2str(2^plotLims(2)))} ];
tickLabels = cellfun(@(x) sprintf('%sx',x),tickLabels,'UniformOutput',false);

% Percentage of no-FM trials
perNoFM_off = cellfun(@(x) (sum(~x(:,:,1))/length(x(:,:,1)))*100, trialsOverThresh_off_VO_perCell);
perNoFM_V = cellfun(@(x,y) (sum(~x(y,:,1))/length(x(y,:,1)))*100, trialsOverThresh_CMN_VO_perCell, VO.FRs_trialType_optoTrials);

% Get scatter sizes for plot
scatSizeLims = [ 50 130 ]; % scatter sizes for the number limits below
numLims = [ 10 90 ]; % number limits for the scatter sizes above
scatSizes_V = getZScoreScatterSizes(100-perNoFM_V,scatSizeLims,numLims);
scatSizes_off = getZScoreScatterSizes(100-perNoFM_off,scatSizeLims,numLims);
scatSizes = [ scatSizes_V(curInd_V) ; scatSizes_off(curInd_off) ];

cols = [ repmat([1 0 0],sum(curInd_V),1) ; repmat([0 0 1],sum(curInd_off),1) ];
BSind = [ VO.BSunits(curInd_V) ; VO.BSunits(curInd_off) ];

% Scatter
figure
hold all

FRch_lowMv_lims = FRch_lowMv;
FRch_lowMv_lims(FRch_lowMv_lims<plotLims(1)) = plotLims(1);
FRch_lowMv_lims(FRch_lowMv_lims>plotLims(2)) = plotLims(2);
FRch_highMv_lims = FRch_highMv;
FRch_highMv_lims(FRch_highMv_lims<plotLims(1)) = plotLims(1);
FRch_highMv_lims(FRch_highMv_lims>plotLims(2)) = plotLims(2);

scatter(FRch_lowMv_lims(BSind),FRch_highMv_lims(BSind),scatSizes(BSind),cols(BSind,:),'filled','^','MarkerFaceAlpha',0.5)
scatter(FRch_lowMv_lims(~BSind),FRch_highMv_lims(~BSind),scatSizes(~BSind),cols(~BSind,:),'filled','MarkerFaceAlpha',0.5)
plot(plotLims,plotLims,'k:','linewidth',1);
plot([0 0],plotLims,'k:','linewidth',0.5)
plot(plotLims,[0 0],'k:','linewidth',0.5)
plot(plotLims,plotLims,'k:','linewidth',1);
plot([0 0],plotLims,'k:','linewidth',0.5)
plot(plotLims,[0 0],'k:','linewidth',0.5)
scatter(2.5,0.5,scatSizeLims(1),'k','filled')
scatter(3,0.5,scatSizeLims(2),'k','filled')

set(gca,'xtick',FRchangeLims(1):1:FRchangeLims(2),'xticklabels',tickLabels,...
    'fontsize',16,'linewidth',1,'tickdir','out','YLim',plotLims,...
    'XLim',plotLims,'ytick',plotLims(1):1:plotLims(2),'yticklabels',tickLabels,...
    'Color','none','clipping','off');
daspect([1 1 1])

% savename = fullfile(figSaveFol,'low_vs_high_V_and_off');
% saveFig(savename,'-dpdf')

%% Plot - x-y differences (between FR changes)

% FR change for low and high movement trials
FRch_temp_lowMvmt = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_noMvmt, [-inf inf]);
FRch_temp_highMvmt = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_mvmt, [-inf inf]);
FRch_temp_lowMvmt_off = getFRchangeVals_log2(VO.FRs_off_baseline_CMNlen, VO.FRs_off_act_CMNlen, index_noMvmt_off, [-inf inf]);
FRch_temp_highMvmt_off = getFRchangeVals_log2(VO.FRs_off_baseline_CMNlen,VO.FRs_off_act_CMNlen,index_mvmt_off, [-inf inf]);

yLims = [-3 3];
yTick = -3:3:3;

diff_FRchange_off = FRch_temp_highMvmt_off(curInd_off) - FRch_temp_lowMvmt_off(curInd_off);
diff_FRchange_V = FRch_temp_highMvmt(curInd_V) - FRch_temp_lowMvmt(curInd_V);
diff_all = [diff_FRchange_V ; diff_FRchange_off];
diff_all_lims = diff_all;
diff_all_lims(diff_all_lims<yLims(1)) = yLims(1);
diff_all_lims(diff_all_lims>yLims(2)) = yLims(2);

pval = ranksum(diff_FRchange_V, diff_FRchange_off);

xPoints = [ ones(length(diff_FRchange_V),1) ; ones(length(diff_FRchange_off),1)*2 ];
numPoints = length(xPoints);
divider = 5;
jitterX = rand(1,numPoints) / divider - (1/(divider*2));
xPoints = xPoints + jitterX';

% Get scatter sizes for plot
scatSizeLims = [ 50 140 ]; % scatter sizes for the number limits below
numLims = [ 10 90 ]; % number limits for the scatter sizes above
scatSizes_V = getZScoreScatterSizes(100-perNoFM_V,scatSizeLims,numLims);
scatSizes_off = getZScoreScatterSizes(100-perNoFM_off,scatSizeLims,numLims);
scatSizes = [ scatSizes_V(curInd_V) ; scatSizes_off(curInd_off) ];

cols = [ repmat([1 0 0],sum(curInd_V),1) ; repmat([0 0 1],sum(curInd_off),1) ];
BSind = [ VO.BSunits(curInd_V) ; VO.BSunits(curInd_off) ];

figure('Position',[2806 766 162 386])
hold all

scatter(xPoints(BSind),diff_all_lims(BSind),scatSizes(BSind),cols(BSind,:),'filled','^','MarkerFaceAlpha',0.5)
scatter(xPoints(~BSind),diff_all_lims(~BSind),scatSizes(~BSind),cols(~BSind,:),'filled','MarkerFaceAlpha',0.5)
med_V = median(diff_FRchange_V);
plot([.7 1.3],[med_V med_V],'k','linewidth',3)
med_off = median(diff_FRchange_off);
plot([1.7 2.3],[med_off med_off],'k','linewidth',3)
plot([.5 2.5],[0 0],'k:','linewidth',1)

scatter(1.5,-2.5,scatSizeLims(1),'k','filled')
scatter(2,-2.5,scatSizeLims(2),'k','filled')

set(gca,'xtick',[1 2],'XLim',[.5 2.5],...
    'fontsize',16,'linewidth',1,'tickdir','out','YLim',yLims,...
    'ytick',yTick,...
    'Color','none','clipping','off');

mean_V = mean(diff_FRchange_V,'omitnan');
SD_V = std(diff_FRchange_V);
fprintf('\nV mean +- SD: %s +- %s',num2str(mean_V),num2str(SD_V))

mean_A = mean(diff_FRchange_off);
SD_A = std(diff_FRchange_off);
fprintf('\nA mean +- SD: %s +- %s',num2str(mean_A),num2str(SD_A))

% savename = fullfile(figSaveFol,'low_vs_high_V_and_off_diff');
% saveFig(savename,'-dpdf')

