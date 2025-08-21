
saveFigFol = '\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3\S8\figSave';

cd('\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3\S8')

miceToUse = true(length(mouseIDforEachRec_VO),1);
% miceToUse( contains(mouseIDforEachRec_VO,'TRO49') ...
%     | contains(mouseIDforEachRec_VO,'TRO51')) = false;
miceToUse_perCell = repelem(miceToUse,numCellPerRec_VO,1);

primarySecondaryInd = true(length(miceToUse_perCell),1);

curInd = GLMfil & visRespCells_GLM & miceToUse_perCell & primarySecondaryInd;
curInd_notVisResp = GLMfil & ~visRespCells_GLM & miceToUse_perCell & primarySecondaryInd;

%% Stats for figures

% Significance values for each predictor, for visually responsive cells 
pVals_allPredictors_col = [ pVal_V_FDR(curInd) ; pVal_FM_FDR(curInd) ; pVal_VxFM_FDR(curInd) ];
pVals_allPredictors_cell = {pVal_V_FDR(curInd) ...
    pVal_FM_FDR(curInd) pVal_VxFM_FDR(curInd) };
% Significance below alpha
pSig_col = pVals_allPredictors_col < alphaLevel;
pSig_cell = cellfun(@(x) x < alphaLevel , pVals_allPredictors_cell,'UniformOutput',false);
% Number of significant cells for each predictor
numSig_V = sum(pSig_cell{1});
numSig_FM = sum(pSig_cell{2});
numSig_VxFM = sum(pSig_cell{3});
% Percentage of significant cells relative to total number of visually responsive cells, for each predictor
numCells_GLM = sum(curInd);
perSig_V = ( numSig_V / numCells_GLM ) * 100;
perSig_FM = ( numSig_FM / numCells_GLM ) * 100;
perSig_VxFM = ( numSig_VxFM / numCells_GLM ) * 100;

% Numbers of cells/rec/mice - visually responsive cells
numCellsTemp = sum(curInd);
numRecTemp = length(unique(VO.recNumber(curInd)));
numMouseTemp = length(unique(VO.mouseID(curInd)));
fprintf('\n# visually responsive cells/rec/mice: %s, %s, %s',num2str(numCellsTemp),num2str(numRecTemp),num2str(numMouseTemp));
numCells_GLM_V = numCellsTemp; numRecs_GLM_V = numRecTemp; numMice_GLM_V = numMouseTemp;

% Numbers of cells/rec/mice - all possible GLM cells
numCellsTemp = sum(curInd | curInd_notVisResp);
numRecTemp = length(unique(VO.recNumber(curInd | curInd_notVisResp)));
numMouseTemp = length(unique(VO.mouseID(curInd | curInd_notVisResp)));
fprintf('\n# possible GLM cells/rec/mice: %s, %s, %s',num2str(numCellsTemp),num2str(numRecTemp),num2str(numMouseTemp))

%% Plot - imagesc of visually responsive cells - according to the GLM

figure

orderVec = abs(FRzScore_Vresp); % order results by these numbers (low to high)
ColAxisLims = [-1 1]; % colorbar axis limits
tAx = VO.binCenters_VO; % time axis for imagesc
plotBetween = [ -2.5 2.5 ]; % plot between these times
xTicks = -2:1:2; % put xticks at these locations
[sortedImsc,outVecOrder] = plotImagesc(PSTH_Vtrial_norm,curInd,orderVec,ColAxisLims,tAx,plotBetween,xTicks,colMap_RB_exp);
title({'V resp. units (increased FR)';'ordered by Z-score'})
xlabel('Time relative to CMN start (sec)')
ylabel('Unit #');
colorbar
drawnow

% saveFig(fullfile(saveFigFol,'GLM_visRespImsc'),'-dpdf');

%% Plot - beta weights and significance

% y-data for plotting - beta values
yData_raw_GLM = [ beta_V_full_lim(curInd) ; beta_FM_full_lim(curInd) ; beta_VxFM_full_lim(curInd) ];
yData_raw_GLM_paired = [ beta_V_full_lim(curInd) , beta_FM_full_lim(curInd) , beta_VxFM_full_lim(curInd) ];
yData_raw_GLM_cell = { beta_V_full_lim(curInd) beta_FM_full_lim(curInd) beta_VxFM_full_lim(curInd) };

BS_GLM = repmat(VO.BSunits(curInd),3,1);

% x-data for plotting - predictor group
divider = 5;
jitterX = rand(1,sum(curInd)*3) / divider - (1/(divider*2)); % add some x-axis jitter
xDataInd_GLM = repelem(1:3,repmat(sum(curInd),3,1));
xData_GLM = xDataInd_GLM + jitterX;

figure('Position',[2214 657 560 525])
hold all

% non-significant
scatter(xData_GLM(~pSig_col & BS_GLM)',yData_raw_GLM(~pSig_col & BS_GLM),100,'k','filled','^','MarkerFaceAlpha',0.3);
scatter(xData_GLM(~pSig_col & ~BS_GLM)',yData_raw_GLM(~pSig_col & ~BS_GLM),100,'k','filled','MarkerFaceAlpha',0.3);
% significant
scatter(xData_GLM(pSig_col & BS_GLM)',yData_raw_GLM(pSig_col & BS_GLM),100,'r','filled','^','MarkerFaceAlpha',0.3);
scatter(xData_GLM(pSig_col & ~BS_GLM)',yData_raw_GLM(pSig_col & ~BS_GLM),100,'r','filled','MarkerFaceAlpha',0.3);
% other
for i = 1:3
    plot([i-.25 i+.25],[median(yData_raw_GLM_cell{i}) median(yData_raw_GLM_cell{i})],'b','Linewidth',3)
end
plot([0.5 3.5],[0 0],'k:','linewidth',1)
set(gca,'Color','none','tickdir','out','linewidth',1,'fontsize',16,...
    'YLim',[-1 1],'XLim',[.5 3.5],'xtick',1:3)

% Some details for the title
cellCountTx = sprintf('# cells: %s, # recs: %s, # mouse: %s',num2str(numCells_GLM_V),num2str(numRecs_GLM_V),num2str(numMice_GLM_V) );
sigTx = sprintf('%% significant: V: %s, FM: %s, VxFM: %s',num2str(round(perSig_V,1)),num2str(round(perSig_FM,1)),num2str(round(perSig_VxFM,1)));
% ANOVA p
pANOVA = friedman(yData_raw_GLM_paired,1,'off');
p_V_FM_bon = signrank(yData_raw_GLM_paired(:,1),yData_raw_GLM_paired(:,2)) * 2;
p_V_VFM_bon = signrank(yData_raw_GLM_paired(:,1),yData_raw_GLM_paired(:,3)) * 2;

% title
title( { sprintf('ANOVA P: %s',num2str(pANOVA)) ; cellCountTx ; sigTx } );

% saveFig(fullfile(saveFigFol,'GLM_betaWeightPlot'),'-dpdf');

%% Plot - V beta vs FM beta scatter

tickLab = {'0' '\geq1'};

curPval = signrank(abs(beta_V_full(curInd)),abs(beta_FM_full(curInd)));

figure
hold all

scatter(abs(beta_V_full_lim(curInd & VO.BSunits)),abs(beta_FM_full_lim(curInd & VO.BSunits)),...
    140,'k','filled',"^",'MarkerFaceAlpha',0.5)
scatter(abs(beta_V_full_lim(curInd & ~VO.BSunits)),abs(beta_FM_full_lim(curInd & ~VO.BSunits)),...
    140,'k','filled','MarkerFaceAlpha',0.5)
plot([0 1],[0 1],'k:','LineWidth',1)

set(gca,'XLim',[0 1],'YLim',[0 1],'fontsize',16,'tickdir','out',...
    'xtick',[0 1],'XTickLabel',tickLab,'Ytick',[0 1],'YTickLabel',tickLab,...
    'Color','none')

xlabel('\beta1 (V)');
ylabel('\beta2 (FM)');
title(sprintf('P = %s',num2str(curPval)))

% saveFig(fullfile(saveFigFol,'scatter_beta_V_FM'),'-dpdf');
