
saveFigFol = '\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3\S4\figSave';

cd('\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3\S4')

miceToUse = true(length(mouseIDforEachRec_VO),1);
% miceToUse( contains(mouseIDforEachRec_VO,'TRO49') ...
%     | contains(mouseIDforEachRec_VO,'TRO51')) = false;
miceToUse_perCell = repelem(miceToUse,numCellPerRec_VO,1);

primarySecondaryInd = true(length(miceToUse_perCell),1);
% primarySecondaryInd = tunedRec_AC_perCell;

indRec = ~avoidRec_VO_FM_perCell & miceToUse_perCell & primarySecondaryInd;

%% Plot - FM scatter - low vs high - 3 PCs

% Data
curPvals = pVals_lowVsHighMvmt_morePCs;
indexNoMvmt = index_noOpto_noMvmt_morePCs;
indexMvmt = index_noOpto_mvmt_morePCs;
FRchange_noMvmt = FRch_base_to_V_lowMvmt_morePCs;
FRchange_mvmt = FRch_base_to_V_highMvmt_morePCs;

tickLabels = [ {sprintf('\x2264%s',num2str(2^FRchangeLims(1)))} cellfun(@num2str,num2cell(2.^(FRchangeLims(1)+1:1:FRchangeLims(2)-1)),'un',0) {sprintf('\x2265%s',num2str(2^FRchangeLims(2)))} ];
tickLabels = cellfun(@(x) sprintf('%sx',x),tickLabels,'UniformOutput',false);
 
% Only looking at cells with increased firing in response to the visual stimulus (all V trials)
curInd = visRespCell_AC & FRch_base_to_V > 0 & indRec;

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
upperSizeBound = max( [ max(ratioOfAll(curInd)) perScatSizes(2) ] ) ;
scatSizes_nonsig_BS = getZScoreScatterSizes(ratioOfAll(curInd & ~sigInd & VO.BSunits),scatterSizeRange,[lowerSizeBound upperSizeBound]);
scatSizes_sig_BS = getZScoreScatterSizes(ratioOfAll(curInd & sigInd & VO.BSunits),scatterSizeRange,[lowerSizeBound upperSizeBound]);
scatSizes_nonsig_NS = getZScoreScatterSizes(ratioOfAll(curInd & ~sigInd & ~VO.BSunits),scatterSizeRange,[lowerSizeBound upperSizeBound]);
scatSizes_sig_NS = getZScoreScatterSizes(ratioOfAll(curInd & sigInd & ~VO.BSunits),scatterSizeRange,[lowerSizeBound upperSizeBound]);
scatSizes_lims = getZScoreScatterSizes(perScatSizes,scatterSizeRange,[perScatSizes(1) perScatSizes(end)]);

fprintf('\nMin ratio: %s, max ratio',num2str(min(ratioOfAll(curInd))),num2str(max(ratioOfAll(curInd))))

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

plot(FRchangeLims,FRchangeLims,'k:','linewidth',1);
plot([0 0],FRchangeLims,'k:','linewidth',0.5)
plot(FRchangeLims,[0 0],'k:','linewidth',0.5)
scatter(2.4,-0.05,scatSizes_lims(1),'b','filled',"^")
scatter(2.7,-0.05,scatSizes_lims(2),'b','filled',"^")
scatter(3,-0.05,scatSizes_lims(3),'b','filled',"^")
scatter(2.4,-0.4,scatSizes_lims(1),'b','filled')
scatter(2.7,-0.4,scatSizes_lims(2),'b','filled')
scatter(3,-0.4,scatSizes_lims(3),'b','filled')
xlabel({'FR change during VIS';'(no twitch)'})
ylabel({'FR change during VIS';'(twitch)'});

numCellsTemp = sum(curInd);
numRecTemp = length(unique(VO.recNumber(curInd)));
numMouseTemp = length(unique(VO.mouseID(curInd)));

title(sprintf('# cells": %s, # recs: %s, # mouse: %s',num2str(numCellsTemp),num2str(numRecTemp),num2str(numMouseTemp) ))

daspect([1 1 1])

% Group stats
FRch_temp_lowMvmt = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_noMvmt_morePCs, [-inf inf]);
FRch_temp_highMvmt = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_mvmt_morePCs, [-inf inf]);
ranksum(FRch_temp_lowMvmt(curInd), FRch_temp_highMvmt(curInd))

% saveFig(fullfile(saveFigFol,'scatter_lowVsHighMvmt_morePCs'),'-dpdf');

