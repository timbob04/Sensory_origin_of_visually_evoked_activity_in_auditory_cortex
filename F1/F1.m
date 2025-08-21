close all

saveFigFol = '\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3\F1\figSave';

cd('\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3\F1');

%% Input variables

recNum = 18; % example recording number
exNum_CMN = 2; % 37; % example trial number from recording number above
exNum_off = 7; % 7 example trial number from recording number above
cLimVal = 40; % color limits for the imagesc
imagescLims_sinRec = [-5 50];

miceToUse = true(length(mouseIDforEachRec_VO),1);
% miceToUse( contains(mouseIDforEachRec_VO,'TRO49') ...
%     | contains(mouseIDforEachRec_VO,'TRO51')) = false;

optoTrials_curRec = VO.FRs_trialType_optoTrials{recFirstPositions_VO(recNum)}; % find opto trials

recInd = ~avoidRec_VO_AvsV & miceToUse;

%% # of recordings and mice

numMice = length(unique(mouseIDforEachRec_VO(recInd)));
fprintf('\nNumber of mice for this analysis: %s',num2str(numMice));
numRec = sum(recInd);
fprintf('\nNumber of rec for this analysis: %s',num2str(numRec));

%% Plot - FM for example trial - CMN

YLIM = [-2 15];
YTICK = [0 15];

figure
hold all

% Get the CMN traces
curDat_zScore_CMN = movementsStore_VO_CMN_noABS{recNum}(exNum_CMN,plotWinInd_VO,1);
plot(tToPlot_VO,curDat_zScore_CMN,'k','linewidth',2)
plot([tToPlot_VO(1) tToPlot_VO(end)],[zScoreThresh zScoreThresh],'c:','linewidth',1.5)
plot([tToPlot_VO(1) tToPlot_VO(end)],[0 0],'k:','linewidth',1.5)
meanResp = mean(movementsStore_VO_CMN_noABS{recNum}(exNum_CMN,CMNactWin_VO,1));
plot(CMNwin_VO,[meanResp meanResp],'Color',colMean,'LineWidth',1.5)

set(gca,'XLim',plotWin_VO,'Color','none',...
    'box','off','linewidth',1,'tickdir','out','fontsize',16,...
    'xtick',xTicks_VO,'YLim',YLIM,'ytick',YTICK)

% saveFig(fullfile(saveFigFol,'FM_singleTrial_CMN'),'-dpdf');

%% Plot - FM for example trial - noise burst

YLIM = [-10 100];
YTICK = [0 100];

figure
hold all

curDat_zScore_CMN = movementsStore_VO_offset_noABS{recNum}(exNum_off,plotWinInd_VO,1);

plot(tToPlot_VO,curDat_zScore_CMN,'k','linewidth',2)
plot([tToPlot_VO(1) tToPlot_VO(end)],[zScoreThresh zScoreThresh],'c:','linewidth',1.5)
plot([tToPlot_VO(1) tToPlot_VO(end)],[0 0],'k:','linewidth',1.5)
meanResp = mean(movementsStore_VO_offset_noABS{recNum}(exNum_off,CMNactWin_VO,1));
plot(CMNwin_VO,[meanResp meanResp],'Color',colMean,'LineWidth',1.5)
set(gca,'XLim',plotWin_VO,'Color','none',...
    'box','off','linewidth',1,'tickdir','out','fontsize',16,...
    'xtick',xTicks_VO,'YLim',YLIM,'ytick',YTICK)

% saveFig(fullfile(saveFigFol,'FM_singleTrial_noiseBurst'),'-dpdf');

%% Plot - imagesc and mean FM PSTH for example recording - CMN

figure('Position',[1719 627.6667 560 695.3333])

% Get PSTH data - mean and SD
curDat_zScore_CMN = movementsStore_VO_CMN_noABS{recNum}(~optoTrials_curRec,:,:);
curDat_zScore_CMN_mean = mean(curDat_zScore_CMN);
curDat_zScore_CMN_SD = std(curDat_zScore_CMN) / sqrt(size(curDat_zScore_CMN,1));
curDat_zScore_CMN_upBound = curDat_zScore_CMN_mean + curDat_zScore_CMN_SD;
curDat_zScore_CMN_lowBound = curDat_zScore_CMN_mean - curDat_zScore_CMN_SD;

orderVal = mean(curDat_zScore_CMN(:,CMNactWin_VO),2);

overThreshold_cur = sum(trialsOverThresh_CMN_VO_noABS{recNum}(~optoTrials_curRec,:,1:numPCsToUse),3) > 0;

% Imagesc
subplot(3,1,2:3)
[~,curOrder] = sort(orderVal); % order of imagesc
imagesc(flipud(curDat_zScore_CMN(curOrder,plotWinInd_VO))) % make imagesc
clim([-cLimVal cLimVal]);
colormap(colMap_MG);
% colorbar
xTickPoints_VO = nan(1,length(xTicks_VO));
for i = 1:length(xTicks_VO)
    xTickPoints_VO(i) = find( abs(tToPlot_VO - xTicks_VO(i)) == min(abs(tToPlot_VO - xTicks_VO(i))) ,1,'first');
end
set(gca,'Color','none','box','off','linewidth',1.5,'fontsize',16,...
    'tickdir','out','xtick',xTickPoints_VO,...
    'xticklabel',xTicks_VO,'clipping','off')
% Indicate trials over threshold
overThresh = flipud(overThreshold_cur(curOrder));
sizeX = sum(plotWinInd_VO);
hold all
for i = 1:size(curDat_zScore_CMN,1)
    if overThresh(i)
        scatter(sizeX+3,i,20,'k','filled');
    end
end
% Indicate example trial
newExPos = find(find(~optoTrials_curRec) == exNum_CMN);
scatter(sizeX+4,find(flipud(curOrder) == newExPos),20,'c','filled');
% other
xlabel('Time relative to VIS start (sec)')
% ylabel('Trial')

% Mean PSTH
subplot(3,1,1)
plot(tToPlot_VO,curDat_zScore_CMN_mean(plotWinInd_VO),'k','linewidth',2);
hold on
xFill = [ tToPlot_VO fliplr(tToPlot_VO) tToPlot_VO(1) ];
yFill = [ curDat_zScore_CMN_upBound(plotWinInd_VO) ...
    fliplr(curDat_zScore_CMN_lowBound(plotWinInd_VO)) ...
    curDat_zScore_CMN_upBound(find(plotWinInd_VO,1,'first')) ];
fill(xFill,yFill,'k','FaceAlpha',0.1,'EdgeAlpha',0)
set(gca,'Color','none','box','off','fontsize',15,'tickdir','out','linewidth',1,...
    'XLim',plotWin_VO,'YLim',imagescLims_sinRec,'ytick',[0 imagescLims_sinRec(2)])
% ylabel({'Movement'; '(z-score)'})

% saveFig(fullfile(saveFigFol,'FM_singleRecording_CMN'),'-dpdf');

%% Plot - imagesc and mean FM PSTH for example recording - noise burst

figure('Position',[1719 627.6667 560 695.3333])

% Get PSTH data - mean and SD

curDat_zScore_noiseBurst = movementsStore_VO_offset_noABS{recNum};
curDat_zScore_noiseBurst_mean = mean(curDat_zScore_noiseBurst);
curDat_zScore_noiseBurst_SD = std(curDat_zScore_noiseBurst) / sqrt(size(curDat_zScore_noiseBurst,1));
curDat_zScore_noiseBurst_upBound = curDat_zScore_noiseBurst_mean + curDat_zScore_noiseBurst_SD;
curDat_zScore_noiseBurst_lowBound = curDat_zScore_noiseBurst_mean - curDat_zScore_noiseBurst_SD;

orderVal = mean(curDat_zScore_noiseBurst(:,CMNactWin_VO),2);

overThreshold_cur = sum(trialsOverThresh_offset_VO_noABS{recNum}(:,:,1:numPCsToUse),3) > 0;

% Imagesc
subplot(3,1,2:3)
[~,curOrder] = sort(orderVal); % order of imagesc
imagesc(flipud(curDat_zScore_noiseBurst(curOrder,plotWinInd_VO))) % make imagesc
clim([-cLimVal cLimVal]);
colormap(colMap_MG);
colorbar
xTickPoints_VO = nan(1,length(xTicks_VO));
for i = 1:length(xTicks_VO)
    xTickPoints_VO(i) = find( abs(tToPlot_VO - xTicks_VO(i)) == min(abs(tToPlot_VO - xTicks_VO(i))) ,1,'first');
end
set(gca,'Color','none','box','off','linewidth',1.5,'fontsize',16,...
    'tickdir','out','xtick',xTickPoints_VO,...
    'xticklabel',xTicks_VO,'clipping','off')
% Indicate trials over threshold
overThresh = flipud(overThreshold_cur(curOrder));
sizeX = sum(plotWinInd_VO);
hold all
for i = 1:size(curDat_zScore_noiseBurst,1)
    if overThresh(i)
        scatter(sizeX+3,i,20,'k','filled');
    end
end
% Indicate example trial
scatter(sizeX+4,find(flipud(curOrder) == exNum_off),20,'c','filled');
% other
xlabel('Time relative to AUD start (sec)')
% ylabel('Trial')

% Mean PSTH
subplot(3,1,1)
plot(tToPlot_VO,curDat_zScore_noiseBurst_mean(plotWinInd_VO),'k','linewidth',2);
hold on
xFill = [ tToPlot_VO fliplr(tToPlot_VO) tToPlot_VO(1) ];
yFill = [ curDat_zScore_noiseBurst_upBound(plotWinInd_VO) ...
    fliplr(curDat_zScore_noiseBurst_lowBound(plotWinInd_VO)) ...
    curDat_zScore_noiseBurst_upBound(find(plotWinInd_VO,1,'first')) ];
fill(xFill,yFill,'k','FaceAlpha',0.1,'EdgeAlpha',0)
set(gca,'Color','none','box','off','fontsize',15,'tickdir','out','linewidth',1,...
    'XLim',plotWin_VO,'YLim',imagescLims_sinRec,'ytick',[0 imagescLims_sinRec(2)])
% ylabel({'Movement'; '(z-score)'})

% saveFig(fullfile(saveFigFol,'FM_singleRecording_noiseBurst'),'-dpdf');

%% Plot - mean movement for all recordings

yLim_linePlot = [-2 30];
ColAxisLims = [-15 15]; % colorbar axis limits
xTicks = -1:1:2; % put xticks at these locations

recIDs = 1:size(movementsStore_VO_CMN_noABS,1);

recIDs_taken = recIDs(recInd);

% Get mean trace data

meanMovements = cellfun(@(x,y) mean(x(~y,:,:)), movementsStore_VO_CMN_noABS,...
    VO.FRs_trialType_optoTrials(recFirstPositions_VO),'UniformOutput',false);
meanMovements_cat_CMN = cat(1,meanMovements{:});
meanMovements = cellfun(@mean,movementsStore_VO_offset_noABS,'UniformOutput',false);
meanMovements_cat_off = cat(1,meanMovements{:});

% Vector for ordering by maximum sound response
orderVal = mean(meanMovements_cat_off(:,CMNactWin_VO,1),2);

% Time point for indicating example recording
tInd = find(tAx_VO > plotWin_VO(1)+0.1,1,'first');

figure('Position',[1928 796 917 451])

subplot(3,2,[3 5])
[~,orderOut] = plotImagesc(meanMovements_cat_CMN(:,:,1),recInd,orderVal,ColAxisLims,tAx_VO,plotWin_VO,xTicks_VO,colMap_MG_exp);
recNumPos = find(orderOut==recNum);
hold on
scatter(tInd,recNumPos,50,'k','filled')
colorbar
pause(1)
xPos_imWithColorbar = get(gca,'Position');

subplot(3,2,[4 6])
plotImagesc(meanMovements_cat_off(:,:,1),recInd,orderVal,ColAxisLims,tAx_VO,plotWin_VO,xTicks_VO,colMap_MG_exp);
colorbar

xFill = [ tAx_VO(plotWinInd_VO) fliplr(tAx_VO(plotWinInd_VO)) tAx_VO(find(plotWinInd_VO,1,'first')) ];

subplot(3,2,1)
resp_mean = mean(meanMovements_cat_CMN(recInd,:));
resp_SD = std(meanMovements_cat_CMN(recInd,:)) / sqrt(sum(recInd));
resp_upper = resp_mean + resp_SD;
resp_lower = resp_mean - resp_SD;
plot(tAx_VO(plotWinInd_VO),resp_mean(plotWinInd_VO),'k','linewidth',2);
hold on
yFill = [ resp_upper(plotWinInd_VO) fliplr(resp_lower(plotWinInd_VO)) resp_upper(find(plotWinInd_VO,1,'first')) ];
fill(xFill,yFill,'k','FaceAlpha',0.15,'EdgeAlpha',0);
set(gca,'fontsize',16,'Xlim',plotWin_VO,'YLim',yLim_linePlot,'Color','none',...
    'box','off','tickdir','out','linewidth',1,'Ytick',[0 yLim_linePlot(2)])
xPos_plot = get(gca,'Position');
set(gca,'Position',[xPos_plot(1) xPos_plot(2) xPos_imWithColorbar(3) xPos_plot(4)])

subplot(3,2,2)
resp_mean = mean(meanMovements_cat_off(recInd,:));
resp_SD = std(meanMovements_cat_off(recInd,:)) / sqrt(sum(recInd));
resp_upper = resp_mean + resp_SD;
resp_lower = resp_mean - resp_SD;
plot(tAx_VO(plotWinInd_VO),resp_mean(plotWinInd_VO),'k','linewidth',2);
hold on
yFill = [ resp_upper(plotWinInd_VO) fliplr(resp_lower(plotWinInd_VO)) resp_upper(find(plotWinInd_VO,1,'first')) ];
fill(xFill,yFill,'k','FaceAlpha',0.15,'EdgeAlpha',0);

set(gca,'fontsize',16,'Xlim',plotWin_VO,'YLim',yLim_linePlot,'Color','none',...
    'box','off','tickdir','out','linewidth',1,'Ytick',[0 yLim_linePlot(2)])
xPos_plot_2 = get(gca,'Position');
set(gca,'Position',[xPos_plot_2(1) xPos_plot(2) xPos_imWithColorbar(3) xPos_plot(4)])

title(sprintf('# mice: %s, # recs: %s',num2str(numMice),num2str(numRec)))

% saveFig(fullfile(saveFigFol,'FM_allRecordings'),'-dpdf');

%% Plot - stats - % of trials above threshold

lims = [ 0 100 ]; % min and max limits for plot

figure
hold all
xTickLab = [ {sprintf('%s%s','\leq',num2str(lims(1)))} {sprintf('%s%s','\geq',num2str(lims(2)))} ];

% Get data
perWin_CMN_VO_lim = perOver_CMN_noABS;
perWin_off_VO_lim = perOver_offset_noABS;

% Statistical test
pVal_per = signrank(perWin_CMN_VO_lim(recInd),perWin_off_VO_lim(recInd))

% Bind data within min and max limits
perWin_CMN_VO_lim(perWin_CMN_VO_lim < lims(1)) = lims(1);
perWin_CMN_VO_lim(perWin_CMN_VO_lim > lims(2)) = lims(2);
perWin_off_VO_lim(perWin_off_VO_lim < lims(1)) = lims(1);
perWin_off_VO_lim(perWin_off_VO_lim > lims(2)) = lims(2);

% Index to example recording
recNumInd = false(length(perWin_off_VO_lim),1);
recNumInd(recNum) = true;

% Plot all data except example recording
scatter(perWin_off_VO_lim(~recNumInd & recInd),perWin_CMN_VO_lim(~recNumInd & recInd),120,'k','filled',...
    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0);
% Plot example recording
scatter(perWin_off_VO_lim(recNumInd & recInd),perWin_CMN_VO_lim(recNumInd & recInd),120,'r','filled',...
    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0);
% other
plot(lims,lims,'k:','linewidth',1)
set(gca,'YLim',lims,'XLim',lims,'xtick',lims,...
    'ytick',lims,'fontsize',16,'linewidth',1.5,...
    'Color','none','box','off','tickdir','out',...
    'xticklabel',xTickLab,'yticklabel',xTickLab)
xlabel('Auditory');
ylabel('Visual window');

title(sprintf('p val: %s',num2str(pVal_per)))

mean_V = mean(perOver_CMN_noABS(recInd));
SD_V = std(perOver_CMN_noABS(recInd));
fprintf('\nV mean +- SD: %s +- %s',num2str(mean_V),num2str(SD_V))

mean_A = mean(perOver_offset_noABS(recInd));
SD_A = std(perOver_offset_noABS(recInd));
fprintf('\nA mean +- SD: %s +- %s',num2str(mean_A),num2str(SD_A))

% saveFig(fullfile(saveFigFol,'scatter_per'),'-dpdf');

%%  Plot - stats - FM median

lims = [ 10 50 ]; % min and max limits for plot

curTrialsOverThresh_CMN = trialsOverThresh_CMN_VO_noABS; 
curTrialsOverThresh_offset = trialsOverThresh_offset_VO_noABS;
curMaxWin_CMN = maxInWin_VO_noABS;
curMaxWin_offset = maxInWin_VO_offset_noABS;

xTickLab = [ {sprintf('%s%s','\leq',num2str(lims(1)))} {sprintf('%s%s','\geq',num2str(lims(2)))} ];

optoTrialsForEachRec = VO.FRs_trialType_optoTrials(recFirstPositions_VO);

curNumPCs = numPCsToUse;
% Get CMN medians
dataConcat = cellfun(@(x) reshape(x(:,:,1:curNumPCs),[],1),curMaxWin_CMN,'UniformOutput',false); % for each recording, concatenate the 'curMaxWin_CMN' for all the PCs required (curNumPCs)
indicesConcat = cellfun(@(x) reshape(x(:,:,1:curNumPCs),[],1),curTrialsOverThresh_CMN,'UniformOutput',false); % for each recording, concatenate the 'curTrialsOverThresh_CMN' for all the PCs required (curNumPCs)
curFunc_CMN = @(x,y,z) median(x(repmat(~y,curNumPCs,1) & z),'all');
medians_CMN = cellfun(curFunc_CMN , dataConcat , optoTrialsForEachRec , indicesConcat);
% Get offset medians
dataConcat = cellfun(@(x) reshape(x(:,:,1:curNumPCs),[],1),curMaxWin_offset,'UniformOutput',false); 
indicesConcat = cellfun(@(x) reshape(x(:,:,1:curNumPCs),[],1),curTrialsOverThresh_offset,'UniformOutput',false);
curFunc_offset = @(x,z) median(x(z),'all');
medians_offset = cellfun(curFunc_offset , dataConcat , indicesConcat);

curPval = signrank(medians_offset(recInd),medians_CMN(recInd));

medians_CMN_lim = medians_CMN; 
medians_CMN_lim(medians_CMN_lim<lims(1)) = lims(1);
medians_CMN_lim(medians_CMN_lim>lims(2)) = lims(2);
medians_offset_lim = medians_offset; 
medians_offset_lim(medians_offset_lim<lims(1)) = lims(1);
medians_offset_lim(medians_offset_lim>lims(2)) = lims(2);

% Index to example recording
recNumInd = false(length(medians_CMN_lim),1);
recNumInd(recNum) = true;

figure
hold all

% Plot all data except example recording
scatter(medians_offset_lim(~recNumInd & recInd),medians_CMN_lim(~recNumInd & recInd),120,'k','filled',...
    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0);
% Plot example recording
scatter(medians_offset_lim(recNumInd & recInd),medians_CMN_lim(recNumInd & recInd),120,'r','filled',...
    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0);
% other
plot(lims,lims,'k:','linewidth',1)
set(gca,'YLim',lims,'XLim',lims,'xtick',lims,...
    'ytick',lims,'fontsize',16,'linewidth',1.5,...
    'Color','none','box','off','tickdir','out',...
    'xticklabel',xTickLab,'yticklabel',xTickLab)
xlabel('Auditory');
ylabel('Visual window');
title(sprintf('P val =  %s',num2str(curPval)))

mean_V = mean(medians_CMN(recInd));
SD_V = std(medians_CMN(recInd));
fprintf('\nV mean +- SD: %s +- %s',num2str(mean_V),num2str(SD_V))

mean_A = mean(medians_offset(recInd));
SD_A = std(medians_offset(recInd));
fprintf('\nA mean +- SD: %s +- %s',num2str(mean_A),num2str(SD_A))

% saveFig(fullfile(saveFigFol,'scatter_amplitude_median'),'-dpdf');
