
saveFigFol = '\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3\S1\figSave';

cd('\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3\S1');

%% Input variables

recNum = 18; % example recording number
exNum_CMN = 135; % example trial number from recording number above
exNum_off = 7; % 7 example trial number from recording number above

optoTrials_curRec = VO.FRs_trialType_optoTrials{recFirstPositions_VO(recNum)}; % find opto trials  

miceToUse = true(length(mouseIDforEachRec_VO),1);
% miceToUse( contains(mouseIDforEachRec_VO,'TRO49') ...
%     | contains(mouseIDforEachRec_VO,'TRO51')) = false;

recInd = ~avoidRec_VO_AvsV & miceToUse;

% X-tick
tickLabelPoints = -2:2:6;
binCenters = ( bins_thSel_mean(1:end-1) + bins_thSel_mean(2:end) ) / 2;
tickLabels = [ {sprintf('%s%s','\leq',num2str(tickLabelPoints(1)))} cellfun(@num2str,num2cell(tickLabelPoints(2:end-1)),'UniformOutput',false)  {sprintf('%s%s','\geq',num2str(tickLabelPoints(end)))} ];

%% # of recordings and mice

numMice = length(unique(mouseIDforEachRec_VO(recInd)));
fprintf('\nNumber of mice for this analysis: %s',num2str(numMice));
numRec = sum(recInd);
fprintf('\nNumber of rec for this analysis: %s',num2str(numRec));

%% Plot - FM for example trial - CMN

YLIM = [-1 10];
YTICK = [0 10];

figure
hold all

% Get the CMN traces

curDat_zScore_CMN = movementsStore_VO_CMN_noABS{recNum}(exNum_CMN,plotWinInd_VO,1);

plot(tToPlot_VO,curDat_zScore_CMN,'k','linewidth',2)
plot([tToPlot_VO(1) tToPlot_VO(end)],[zScoreThresh zScoreThresh],'c:','linewidth',1.5)
plot([tToPlot_VO(1) tToPlot_VO(end)],[0 0],'k:','linewidth',1.5)
meanResp = mean(movementsStore_VO_CMN_noABS{recNum}(exNum_CMN,CMNactWin_VO,1));
plot([0 .5],[meanResp meanResp],'c','LineWidth',1)
set(gca,'XLim',plotWin_VO,'Color','none',...
    'box','off','linewidth',1,'tickdir','out','fontsize',16,...
    'xtick',xTicks_VO,'YLim',YLIM,'ytick',YTICK);
title(sprintf('exNum_CMN: %s',num2str(exNum_CMN)))

% saveFig(fullfile(saveFigFol,'FM_singleTrial_CMN'),'-dpdf');

%% Plot - FM for example trial - noise burst

YLIM = [-5 80];
YTICK = [0 80];

figure
hold all

curDat_zScore_CMN = movementsStore_VO_offset_noABS{recNum}(exNum_off,plotWinInd_VO,1);

plot(tToPlot_VO,curDat_zScore_CMN,'k','linewidth',2)
plot([tToPlot_VO(1) tToPlot_VO(end)],[zScoreThresh zScoreThresh],'c:','linewidth',1.5)
plot([tToPlot_VO(1) tToPlot_VO(end)],[0 0],'k:','linewidth',1.5)
meanResp = mean(movementsStore_VO_offset_noABS{recNum}(exNum_off,CMNactWin_VO,1));
plot([0 .5],[meanResp meanResp],'c','LineWidth',1)
set(gca,'XLim',plotWin_VO,'Color','none',...
    'box','off','linewidth',1,'tickdir','out','fontsize',16,...
    'xtick',xTicks_VO,'YLim',YLIM,'ytick',YTICK);
title(sprintf('exNum_off: %s',num2str(exNum_off)))

% saveFig(fullfile(saveFigFol,'FM_singleTrial_off'),'-dpdf');

%% Get normalized counts - percentage of total trials

h_CMN_norm = (h_CMN_mean ./ sum(h_CMN_mean,2))*100; % Percentage of total count
h_CMN_norm_mean = mean( h_CMN_norm(recInd,:,:) , 1 );
h_off_norm = (h_off_mean ./ sum(h_off_mean,2))*100; % Percentage of total count
h_off_norm_mean = mean( h_off_norm(recInd,:,:) , 1 );

%% Plot - normalized counts - example recording

figure

b = bar( binCenters, [ h_CMN_norm(recNum,:,1) ; h_off_norm(recNum,:,1) ] );

b(1).FaceColor = 'r';
b(2).FaceColor = 'b';
b(1).EdgeColor = 'none';
b(2).EdgeColor = 'none';

set(gca,'Box','off','linewidth',1,'fontsize',15,'tickdir','out',...
    'xtick',tickLabelPoints,'XTickLabel',tickLabels,'Color','none',...
    'XTickLabelRotation',90,'Xlim',[tickLabelPoints(1) tickLabelPoints(end)]);
xlabel('Max FM (bs; z-score) ')
ylabel('Normalized trial count')

% saveFig(fullfile(saveFigFol,'thresholdCounts_exRec_mean'),'-dpdf');

%% Plot - normalized counts - all recordings - mean

figure
hold all

% Plot CMN individual recordings
for i = 1:numRec_VO
    if recInd(i)
        plot(binCenters,h_CMN_norm(i,:,1),'r','linewidth',0.5);
    end
end
% Plot CMN recording mean
plot(binCenters,h_CMN_norm_mean(:,:,1),'r-o','linewidth',2)
% Plot noise burst individual recordings
for i = 1:numRec_VO
    if recInd(i)
        plot(binCenters,h_off_norm(i,:,1),'b','linewidth',0.5);
    end
end
% Plot noise burst recording mean
plot(binCenters,h_off_norm_mean(:,:,1),'b-o','linewidth',2)
% Setting and save
set(gca,'Box','off','linewidth',1,'fontsize',15,'tickdir','out',...
    'xtick',tickLabelPoints,'XTickLabel',tickLabels,'Color','none',...
    'XTickLabelRotation',90,'XLim',[tickLabelPoints(1) tickLabelPoints(end)]);

% saveFig(fullfile(saveFigFol,'thresholdCounts_allRec_mean'),'-dpdf');

