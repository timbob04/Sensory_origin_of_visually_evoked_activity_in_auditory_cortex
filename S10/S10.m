close all

saveFigFol = '\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3\S10\figFol';

%% WT and TRO44 data

% WT data
indRec = ~avoidRec_VO_FM_perCell & VO.wildTypeMouse;
curInd = visRespCell_AC & FRch_base_to_V > 0 & indRec;
FRchange_noOpto = getFRchangeVals(VO.FRs_baseline,VO.FRs_Vonset,indexToNonOptoCells);
FRchange_opto = getFRchangeVals(VO.FRs_baseline,VO.FRs_Vonset,indexToOptoCells);
FRch_WT_V = FRchange_noOpto(curInd);
FRch_WT_opto = FRchange_opto(curInd);
BS_WT = VO.BSunits(curInd);

% Include TRO44 data
indTRO44 = TRO44.FRchanges(:,1) > 0;
FRch_V = [FRch_WT_V ; 2.^TRO44.FRchanges(indTRO44,1)];
FRch_VO = [FRch_WT_opto ; 2.^TRO44.FRchanges(indTRO44,2)];
BS_cur = [BS_WT ; TRO44.BSunit(indTRO44)];

[ median(FRch_V) std(FRch_V) ; ...
    median(FRch_VO) std(FRch_VO) ] 

pVal_WT = signrank(FRch_V, FRch_VO);

%% Scatter - V vs VO

lims = [-1 3];

FRch_V_lims = log2(FRch_V);
FRch_V_lims(FRch_V_lims<lims(1)) = lims(1);
FRch_V_lims(FRch_V_lims>lims(2)) = lims(2);
FRch_VO_lims = log2(FRch_VO);
FRch_VO_lims(FRch_VO_lims<lims(1)) = lims(1);
FRch_VO_lims(FRch_VO_lims>lims(2)) = lims(2);

figure
hold all
scatter(FRch_V_lims(BS_cur),FRch_VO_lims(BS_cur),...
    50,'filled',"^",'MarkerFaceColor','k',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
scatter(FRch_V_lims(~BS_cur),FRch_VO_lims(~BS_cur),...
    50,'filled','MarkerFaceColor','k',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
set(gca,'XLim',lims,'YLim',lims,'tickdir','out');
plot(lims,lims,'k:')

% saveFig(fullfile(saveFigFol,'scatter_V_vs_VO'),'-dpdf')

numCellsTemp = sum(curInd) + sum(indTRO44)
numRecTemp = length(unique(VO.recNumber(curInd))) + length(unique(TRO44.recNum(indTRO44)))
numMouseTemp = length(unique(VO.mouseID(curInd))) + 1

%% Box plot

diffLims = [-2 2];
xLim = [.6 1.4]; 

diffs = FRch_VO - FRch_V;
medDiffs = median(diffs(~isinf(diffs)),'omitnan');
diffs_lim = diffs;
diffs_lim(diffs_lim<diffLims(1)) = diffLims(1);
diffs_lim(diffs_lim>diffLims(2)) = diffLims(2);
numPoints = length(diffs);

yticks = [diffLims(1) 0 diffLims(2)];
yTicksLabels = { sprintf('%s%s','\leq',num2str(diffLims(1))) '0' sprintf('%s%s','\geq',num2str(diffLims(2))) };

figure('Position',[1440 918 170 420])
hold all
divider = 5;
jitter = (rand(numPoints,1) / divider) - (1/divider/2);
xData = ones(numPoints,1) + jitter;
scatter(xData(BS_cur),diffs_lim(BS_cur),100,'filled',"^",'MarkerFaceColor','k',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
scatter(xData(~BS_cur),diffs_lim(~BS_cur),100,'filled','MarkerFaceColor','k',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
plot(xLim,[0 0],'k:')
plot(xLim,[medDiffs medDiffs],'r')
set(gca,'XLim',xLim,'YLim',diffLims,'tickdir','out',...
    'LineWidth',1,'fontsize',15,'xtick','','Color','none',...
    'ytick',yticks,'yticklabel',yTicksLabels);
% 
% saveFig(fullfile(saveFigFol,'boxPlot'),'-dpdf')






