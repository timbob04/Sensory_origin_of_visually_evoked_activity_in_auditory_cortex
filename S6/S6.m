close all

saveFigFol = '\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3\S6\figSave';

cd('\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3\S6');

%% Input variables

miceToUse = true(length(mouseIDforEachRec_VO),1);
miceToUse_perCell = repelem(miceToUse,numCellPerRec_VO,1);

primarySecondaryInd = true(length(miceToUse_perCell),1);

indRec = ~avoidRec_VO_FM_perCell & miceToUse_perCell & primarySecondaryInd;

recInd = ~avoidRec_VO_AvsV;

minPerVFMtrials = 5;

%% Percentage of FM trials for A vs V

lims = [ 0 100 ]; % min and max limits for plot

figure
hold all
xTickLab = [ {sprintf('%s%s','\leq',num2str(lims(1)))} {sprintf('%s%s','\geq',num2str(lims(2)))} ];

% Get data
perWin_CMN_VO_lim = perOver_CMN_longer;
perWin_off_VO_lim = perOver_off_longer;

% Statistical test
pVal_per = signrank(perWin_CMN_VO_lim(recInd),perWin_off_VO_lim(recInd));

% Bind data within min and max limits
perWin_CMN_VO_lim(perWin_CMN_VO_lim < lims(1)) = lims(1);
perWin_CMN_VO_lim(perWin_CMN_VO_lim > lims(2)) = lims(2);
perWin_off_VO_lim(perWin_off_VO_lim < lims(1)) = lims(1);
perWin_off_VO_lim(perWin_off_VO_lim > lims(2)) = lims(2);
% Plot all data except example recording
scatter(perWin_off_VO_lim(recInd),perWin_CMN_VO_lim(recInd),120,'k','filled',...
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

% saveFig(fullfile(saveFigFol,'scatter_per'),'-dpdf');

%% Plot - FM scatter - low vs high

% Data
curPvals = pVals_lowVsHighMvmt_longerB;
indexNoMvmt = index_noOpto_noMvmt_longerB;
indexMvmt = index_noOpto_mvmt_longerB;
FRchange_noMvmt = FRch_base_to_V_lowMvmt_longerB;
FRchange_mvmt = FRch_base_to_V_highMvmt_longerB;

tickLabels = [ {sprintf('\x2264%s',num2str(2^FRchangeLims(1)))} cellfun(@num2str,num2cell(2.^(FRchangeLims(1)+1:1:FRchangeLims(2)-1)),'un',0) {sprintf('\x2265%s',num2str(2^FRchangeLims(2)))} ];
tickLabels = cellfun(@(x) sprintf('%sx',x),tickLabels,'UniformOutput',false);

% Only looking at cells with increased firing in response to the visual stimulus (all V trials)
perOver_CMN_longer_perCell = repelem(perOver_CMN_longer,numCellPerRec_VO,1);
curInd = visRespCell_AC & FRch_base_to_V > 0 & indRec & perOver_CMN_longer_perCell >= minPerVFMtrials;

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
scatter(FRchange_noMvmt(curInd & ~sigInd & VO.BSunits),FRchange_mvmt(curInd & ~sigInd & VO.BSunits),...
    scatSizes_nonsig_BS,'filled',"^",'MarkerFaceColor','k',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
scatter(FRchange_noMvmt(curInd & ~sigInd & ~VO.BSunits),FRchange_mvmt(curInd & ~sigInd & ~VO.BSunits),...
    scatSizes_nonsig_NS,'filled','MarkerFaceColor','k',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
% Significant
scatter(FRchange_noMvmt(curInd & sigInd & VO.BSunits),FRchange_mvmt(curInd & sigInd & VO.BSunits),...
    scatSizes_sig_BS,'filled',"^",'MarkerFaceColor','r',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
scatter(FRchange_noMvmt(curInd & sigInd & ~VO.BSunits),FRchange_mvmt(curInd & sigInd & ~VO.BSunits),...
    scatSizes_sig_NS,'filled','MarkerFaceColor','r',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);

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

% Group-level stat
FRch_lowMvmt = getFRchangeVals_log2(VO.FRs_baseline,VO.FRs_Vonset,index_noOpto_noMvmt_longerB,[-inf inf]);
FRch_highMvmt = getFRchangeVals_log2(VO.FRs_baseline,VO.FRs_Vonset,index_noOpto_mvmt_longerB,[-inf inf]);
signrank(FRch_lowMvmt(curInd),FRch_highMvmt(curInd))

% saveFig(fullfile(saveFigFol,'scatter_lowVsHighMvmt'),'-dpdf');
