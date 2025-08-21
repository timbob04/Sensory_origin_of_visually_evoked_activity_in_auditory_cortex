close all

saveFigFol = '\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3\S3\figSave';

cd('\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3\S3')

%% Input variables

recNum = 3; % example recording number 3
exNum = 14; % example trial number from recording number above 14

miceToUse = true(numRecIn_AVO,1);
% miceToUse( contains(mouseIDforEachRec_AVO,'TRO34') ...
%     | contains(mouseIDforEachRec_AVO,'TS036') ...
%     | contains(mouseIDforEachRec_AVO,'NaN') ) = false;

recInd = miceToUse & ~avoidRec_AVO;

AVtrials = AVO.AVtrials{recNum};

%% Plot - FM for example trial

YLIM = [-5 100];
YTICK = [0 100];

figure
hold all

curDat_zScore = movementsStore_AVO_noABS{recNum}(AVtrials,:,1);
plot(tAx_AVO(plotWinInd_AVO),curDat_zScore(exNum,plotWinInd_AVO),'k','linewidth',2)
plot([plotWin_AVO(1) plotWin_AVO(2)],[zScoreThresh zScoreThresh],'c:','linewidth',1.5)
plot([plotWin_AVO(1) plotWin_AVO(end)],[0 0],'k:','linewidth',1.5)
set(gca,'XLim',plotWin_AVO,'Color','none',...
    'box','off','linewidth',1,'tickdir','out','ytick',YTICK,'fontsize',16,...
    'xtick',xTicks_AVO,'YLim',YLIM)

% saveFig(fullfile(saveFigFol,'FM_singleTrial'),'-dpdf');

%% Plot - imagesc - example recording

cLimVal = 35; % color limits for the imagesc
YLIM = [-2 30];
ytick = [0 30];

figure('Position',[1719 627.6667 560 695.3333])

% Get PSTH data - mean and SD

curDat_zScore = movementsStore_AVO_noABS{recNum}(AVtrials,:,1);
curDat_AVO_zScore_mean_ex = mean(curDat_zScore);
curDat_AVO_zScore_SD_ex = std(curDat_zScore) / sqrt(numTrials_AVO);
curDat_AVO_zScore_upBound_ex = curDat_AVO_zScore_mean_ex + curDat_AVO_zScore_SD_ex;
curDat_AVO_zScore_lowBound_ex = curDat_AVO_zScore_mean_ex - curDat_AVO_zScore_SD_ex;

mxWin_CMN = meanInWin_AVO_noABS{recNum,1}(AVtrials);
mxWin_RDS = meanInWin_AVO_noABS{recNum,2}(AVtrials);

maxWin_AVO = max([mxWin_CMN mxWin_RDS],[],2);

% Imagesc
subplot(3,1,2:3)
[~,curOrder] = sort(maxWin_AVO); % order of imagesc
imagesc(flipud(curDat_zScore(curOrder,plotWinInd_AVO))) % make imagesc
clim([-cLimVal cLimVal]);
% colorbar
colormap(colMap_MG);
xTickPoints_AVO = nan(1,length(xTicks_AVO));
for i = 1:length(xTicks_AVO)
    xTickPoints_AVO(i) = find( abs(plotWinTime_AVO - xTicks_AVO(i)) == min(abs(plotWinTime_AVO - xTicks_AVO(i))) ,1,'first');
end
set(gca,'Color','none','box','off','linewidth',1.5,'fontsize',16,...
    'tickdir','out','xtick',xTickPoints_AVO,...
    'xticklabel',xTicks_AVO,'clipping','off')
% Indicate trials with FM over threshold
overThreshold_cur_CMN = trialsOverThresh_AVO_noABS{recNum,1}(AVtrials,:,1);
overThreshold_cur_RDS = trialsOverThresh_AVO_noABS{recNum,2}(AVtrials,:,1);
overThresh_CMN_store_order = flipud(overThreshold_cur_CMN(curOrder));
overThresh_RDS_store_order = flipud(overThreshold_cur_RDS(curOrder));
sizeX = sum(plotWinInd_AVO);
hold all
for i = 1:sum(AVtrials)
    if overThresh_CMN_store_order(i)
        scatter(sizeX+3,i,20,'k','filled');
    end
    if overThresh_RDS_store_order(i)
        scatter(sizeX+7,i,20,'k','filled');
    end
end
% Indicate example trial
scatter(sizeX+4,find(flipud(curOrder) == exNum),20,'c','filled');
% other
xlabel('Time relative to sound start (sec)')
% ylabel('Trial')

% Mean PSTH
subplot(3,1,1)
plot(plotWinTime_AVO,curDat_AVO_zScore_mean_ex(plotWinInd_AVO),'k','linewidth',2);
hold on
xFill = [ plotWinTime_AVO fliplr(plotWinTime_AVO) plotWinTime_AVO(1) ];
yFill = [ curDat_AVO_zScore_upBound_ex(plotWinInd_AVO) ...
    fliplr(curDat_AVO_zScore_lowBound_ex(plotWinInd_AVO)) ...
    curDat_AVO_zScore_upBound_ex(find(plotWinInd_AVO,1,'first')) ];
fill(xFill,yFill,'k','FaceAlpha',0.1,'EdgeAlpha',0)
set(gca,'Color','none','box','off','fontsize',15,'tickdir','out',...
    'linewidth',1,'XLim',plotWin_AVO,'YLim',YLIM,'ytick',ytick)
% ylabel({'Movement'; '(z-score)'})

% saveFig(fullfile(saveFigFol,'FM_singleRecording'),'-dpdf');

%% Plot - imagesc - all recordings

cLimVal = 10; % color limits for the imagesc
YLIM = [-2 10];
ytick = [0 10];

% Get data for plotting
data_bs = cellfun(@(x,y) x(y,:,:) - mean(x(y,baseInd_AVO,:),2) ,...
    movementsStore_AVO_noABS(recInd) , AVO.AVtrials(recInd) , 'UniformOutput',false);
data_zscore = cellfun(@(x) x ./ std(x(:,baseInd_AVO,:),[],2) , ...
    data_bs , 'UniformOutput',false);
meanResp_rec_temp = cellfun(@mean,data_zscore,'UniformOutput',false);
meanResp_rec = cat(1,meanResp_rec_temp{:});
allRec_mean = mean(meanResp_rec);
allRec_SD = std(meanResp_rec) / sqrt(size(meanResp_rec,1)) ;
allRec_upperBound = allRec_mean + allRec_SD;
allRec_lowerBound = allRec_mean - allRec_SD;

actTimeInd = plotWinTime_AVO > -0.5 & plotWinTime_AVO < 0.5; % for ordering the rows by the max response

figure('Position',[1719 627.6667 560 695.3333])

% Imagesc
subplot(3,1,2:3)
[~,curOrder] = sort(mean(meanResp_rec(:,actTimeInd),2));
imagesc(flipud(meanResp_rec(curOrder,plotWinInd_AVO,1))) % make imagesc
hold all
clim([-cLimVal cLimVal]);
colormap(colMap_MG);
% colorbar
xTickPoints_AVO = nan(1,length(xTicks_AVO));
for i = 1:length(xTicks_AVO)
    xTickPoints_AVO(i) = find( abs(plotWinTime_AVO - xTicks_AVO(i)) == min(abs(plotWinTime_AVO - xTicks_AVO(i))) ,1,'first');
end
set(gca,'Color','none','box','off','linewidth',1.5,'fontsize',16,...
    'tickdir','out','xtick',xTickPoints_AVO,...
    'xticklabel',xTicks_AVO,'clipping','off')
% Indicate example recording
recsAllID = 1:numRecIn_AVO; % all recordings
recsAllID = recsAllID(recInd); % only recordings in plot
recsAllID = recsAllID(curOrder); % reorder by imagesc order
exRecInd = find(fliplr(recsAllID) == recNum); % find position of example recording in all recording imagesc
scatter(sizeX+4,exRecInd,20,'c','filled');
% other
xlabel('Time relative to sound start (sec)')
% ylabel('Trial')

% Mean PSTH
subplot(3,1,1)
plot(plotWinTime_AVO,allRec_mean(plotWinInd_AVO),'k','linewidth',2);
hold on
xFill = [ plotWinTime_AVO fliplr(plotWinTime_AVO) plotWinTime_AVO(1) ];
yFill = [ allRec_upperBound(plotWinInd_AVO) ...
    fliplr(allRec_lowerBound(plotWinInd_AVO)) ...
    allRec_upperBound(find(plotWinInd_AVO,1,'first')) ];
fill(xFill,yFill,'k','FaceAlpha',0.1,'EdgeAlpha',0)
set(gca,'Color','none','box','off','fontsize',15,'tickdir','out',...
    'linewidth',1,'XLim',plotWin_AVO,'YLim',YLIM,'ytick',ytick)
% ylabel({'Movement'; '(z-score)'})

% saveFig(fullfile(saveFigFol,'FM_allRecordings'),'-dpdf');

%% Plot - stats - % of trials above threshold

figure
hold all

lims = [ 0 50 ];
xTickLab = [ {sprintf('%s%s','\leq',num2str(lims(1)))} {sprintf('%s%s','\geq',num2str(lims(2)))} ];

% Bind data within min and max limits
perWin_CMN_AVO_lim = perOver_AVO_noABS(:,1);
perWin_RDS_AVO_lim = perOver_AVO_noABS(:,2);

pVal_per = signrank(perWin_CMN_AVO_lim(recInd),perWin_RDS_AVO_lim(recInd));

perWin_RDS_AVO_lim(perWin_RDS_AVO_lim < lims(1)) = lims(1);
perWin_RDS_AVO_lim(perWin_RDS_AVO_lim > lims(2)) = lims(2);
perWin_CMN_AVO_lim(perWin_CMN_AVO_lim < lims(1)) = lims(1);
perWin_CMN_AVO_lim(perWin_CMN_AVO_lim > lims(2)) = lims(2);

% Index to example recording
recNumInd = false(numRecIn_AVO,1);
recNumInd(recNum) = true;

% Plot all data except example recording
scatter(perWin_RDS_AVO_lim(~recNumInd & recInd),perWin_CMN_AVO_lim(~recNumInd & recInd),120,'k','filled',...
    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0);
% Plot example recording
scatter(perWin_RDS_AVO_lim(recNumInd & recInd),perWin_CMN_AVO_lim(recNumInd & recInd),120,'r','filled',...
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

mean_V = mean(perOver_AVO_noABS(recInd,1));
SD_V = std(perOver_AVO_noABS(recInd,1));
fprintf('\nV mean +- SD: %s +- %s',num2str(mean_V),num2str(SD_V))

mean_A = mean(perOver_AVO_noABS(recInd,2));
SD_A = std(perOver_AVO_noABS(recInd,2));
fprintf('\nA mean +- SD: %s +- %s',num2str(mean_A),num2str(SD_A))

% saveFig(fullfile(saveFigFol,'scatter_per'),'-dpdf');

%%  Plot - stats - FM median

lims = [ 0 40 ]; % min and max limits for plot

curTrialsOverThresh_CMN = trialsOverThresh_AVO_noABS(:,1); 
curTrialsOverThresh_RDS = trialsOverThresh_AVO_noABS(:,2);
curMaxWin_CMN = maxInWin_AVO_noABS(:,1);
curMaxWin_RDS = maxInWin_AVO_noABS(:,2);

xTickLab = [ {sprintf('%s%s','\leq',num2str(lims(1)))} {sprintf('%s%s','\geq',num2str(lims(2)))} ];

curNumPCs = numPCsToUse;
% Function to get the data for the curNumPCs, in one column
dataConcat_func = @(x) reshape(x(:,:,1:curNumPCs),[],1);
% Function to get the trials to take, for the curNumPCs, in one column
indicesConcat_func = @(x,y) reshape(x(:,:,1:curNumPCs),[],1) & repmat(y,curNumPCs,1);
% Get the data and indices for the CMN
dataConcat_CMN = cellfun(dataConcat_func,curMaxWin_CMN(recInd),'UniformOutput',false); 
indicesConcat_CMN = cellfun(indicesConcat_func,curTrialsOverThresh_CMN(recInd),AVO.AVtrials(recInd),'UniformOutput',false);
% Get the data and indices for the RDS
dataConcat_RDS = cellfun(dataConcat_func,curMaxWin_RDS(recInd),'UniformOutput',false); 
indicesConcat_RDS = cellfun(indicesConcat_func,curTrialsOverThresh_RDS(recInd),AVO.AVtrials(recInd),'UniformOutput',false);
% Function to get the median for each recording
median_func = @(x,y) median(x(y));
% Median of maximal amplitudes for all AV trials with FM above threshold 
medians_CMN = nan(numRecIn_AVO,1);
medians_CMN(recInd) = cellfun(median_func,dataConcat_CMN,indicesConcat_CMN);
medians_RDS = nan(numRecIn_AVO,1);
medians_RDS(recInd) = cellfun(median_func,dataConcat_RDS,indicesConcat_RDS);

curPval = signrank(medians_RDS(recInd),medians_CMN(recInd));

medians_CMN_lim = medians_CMN; 
medians_CMN_lim(medians_CMN_lim<lims(1)) = lims(1);
medians_CMN_lim(medians_CMN_lim>lims(2)) = lims(2);
medians_RDS_lim = medians_RDS; 
medians_RDS_lim(medians_RDS_lim<lims(1)) = lims(1);
medians_RDS_lim(medians_RDS_lim>lims(2)) = lims(2);

% Index to example recording
recNumInd = false(numRecIn_AVO,1);
recNumInd(recNum) = true;

figure
hold all

% Plot all data except example recording
scatter(medians_RDS_lim(~recNumInd & recInd),medians_CMN_lim(~recNumInd & recInd),120,'k','filled',...
    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0);
% Plot example recording
scatter(medians_RDS_lim(recNumInd & recInd),medians_CMN_lim(recNumInd & recInd),120,'r','filled',...
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

mean_A = mean(medians_RDS(recInd));
SD_A = std(medians_RDS(recInd));
fprintf('\nA mean +- SD: %s +- %s',num2str(mean_A),num2str(SD_A))

% saveFig(fullfile(saveFigFol,'scatter_med'),'-dpdf');

%% # of recordings and mice

numMice = length(unique(mouseIDforEachRec_AVO));
fprintf('\nNumber of mice for this analysis: %s',num2str(numMice));
numRec = sum(recInd);
fprintf('\nNumber of rec for this analysis: %s',num2str(numRec));

