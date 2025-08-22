close all

cd('\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3\S11');

saveFigFol = '\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3\S11\figFol';

%% Inputs

% The median difference between the V and VO conditions (observed condition)
med_VvsVO = -0.3826; % from the FR change values
med_VvsVO_diff = -2.16; % from the differences in FRs

iterToUse = 4;
exUnit = 206;

PSTH_timeWin = [-.5 1];
binSize = 0.05;

bins = PSTH_timeWin(1):binSize:PSTH_timeWin(2);
binCenters = (bins(1:end-1) + bins(2:end)) / 2;

%% Example PSTHs from same cell for ran and other trials

indRanHalf_forRaster = indsRanHalf{exUnit,iterToUse}(indexToNonOptoCells{exUnit});
indOtherHalf_forRaster = ~indRanHalf_forRaster;

% Ran half PSTH
spkInd = ismember(VO.raster_trInd_V{exUnit},find(indRanHalf_forRaster));
spkTimes = VO.raster_times_V{exUnit}(spkInd);
h_ran = histcounts(spkTimes,bins);
% Other half PSTH
spkInd = ismember(VO.raster_trInd_V{exUnit},find(indOtherHalf_forRaster));
spkTimes = VO.raster_times_V{exUnit}(spkInd);
h_other = histcounts(spkTimes,bins);

figure
hold all
plot(binCenters,h_ran,'k','linewidth',1.5)
plot(binCenters,h_other,'r','linewidth',1.5)
set(gca,'XLim',PSTH_timeWin,'tickdir','out','fontsize',15,...
    'Color','none')

saveFig(fullfile(saveFigFol,'PSTHs_exUnit'),'-dpdf')

%% Scatter - FR change for random and other half - example iteration

lims = [-1 3];
xLim = [.6 1.4];
diffLims = [-1.5 1.5];
        
incFR = FRch_base_to_V_ranHalf(:,iterToUse) > 0;
curInd = visRespCell_AC_ranHalf(:,iterToUse) & incFR;
group1 = FRch_base_to_V_ranHalf(curInd,iterToUse);
group2 = FRch_base_to_V_otherHalf(curInd,iterToUse);
diffs = group2 - group1;
diffs_lim = diffs;
diffs_lim(diffs_lim<diffLims(1)) = diffLims(1);
diffs_lim(diffs_lim>diffLims(2)) = diffLims(2);

group1_lims = group1;
group1_lims(group1_lims<lims(1)) = lims(1);
group1_lims(group1_lims>lims(2)) = lims(2);
group2_lims = group2;
group2_lims(group2_lims<lims(1)) = lims(1);
group2_lims(group2_lims>lims(2)) = lims(2);

curBS = VO.BSunits(curInd);

indExUnit = find(curInd) == exUnit;

figure('Position',[1820 872 835 420])

subplot(1,3,1:2)
hold all
scatter(group1_lims(curBS),group2_lims(curBS),...
    50,'filled',"^",'MarkerFaceColor','k',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
scatter(group1_lims(~curBS),group2_lims(~curBS),...
    50,'filled','MarkerFaceColor','k',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
scatter(group1_lims(indExUnit),group2_lims(indExUnit),150,'b')
set(gca,'XLim',lims,'YLim',lims,'tickdir','out','Color','none');
plot(lims,lims,'k:')

% Difference boxplot
subplot(1,3,3)
group1 = FRch_base_to_V_ranHalf_notLog(curInd,iterToUse);
group2 = FRch_base_to_V_otherHalf_notLog(curInd,iterToUse);
diffs = group2 - group1;
diffs_lim = diffs;
diffs_lim(diffs_lim<diffLims(1)) = diffLims(1);
diffs_lim(diffs_lim>diffLims(2)) = diffLims(2);
hold all
divider = 5;
jitter = (rand(length(group1_lims),1) / divider) - (1/divider/2);
xData = ones(length(group1_lims),1) + jitter;
scatter(xData(curBS),diffs_lim(curBS),100,'filled',"^",'MarkerFaceColor','k',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
scatter(xData(~curBS),diffs_lim(~curBS),100,'filled','MarkerFaceColor','k',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
plot(xLim,[0 0],'k:')
plot(xLim,[diffFR_ranOther_inc(iterToUse) diffFR_ranOther_inc(iterToUse)],'r')
set(gca,'XLim',xLim,'YLim',diffLims,'tickdir','out',...
    'LineWidth',1,'fontsize',15,'xtick','','Color','none');

saveFig(fullfile(saveFigFol,'scatter_ranOther'),'-dpdf')

numCellsTemp = sum(curInd);
numRecTemp = length(unique(VO.recNumber(curInd)))
numMouseTemp = length(unique(VO.mouseID(curInd)))

%% Histogram - FR change

diff_vals = diffFR_ranOther_inc_notLog;

figure
hold all
xLim = 0.5;
xtick = [-xLim 0 xLim];
xtickLabs = cellfun(@(x) num2str(x,'%.1f'),num2cell(xtick),'UniformOutput',false);
limPoint = max(abs(diff_vals));
if limPoint > xLim; error('Histogram xLim not correct'); end
bins = -limPoint:0.025:limPoint;
binCenters = ( bins(1:end-1) + bins(2:end) ) / 2 ;
h = histcounts(diff_vals,bins);
bar(binCenters,h);
YLIM = get(gca,'YLim');
plot([0 0],YLIM,'k:','linewidth',1)
plot([med_VvsVO med_VvsVO],YLIM,'b:','linewidth',1)
set(gca,'XLim',[-xLim xLim],'Color','none','box','off','linewidth',1,...
    'tickdir','out','fontsize',14,'xtick',xtick,'xticklabel',xtickLabs)
xlabel('Mean (\DeltaFR_o_t_h_e_r - \DeltaFR_r_a_n_d_o_m)')
ylabel('Number of iterations (1000 total)')
title('Cells with increased V onset FRs')

saveFig(fullfile(saveFigFol,'histogram_change'),'-dpdf')

%% Histogram - difference in FR

diff_vals = diffFR_sub_ranOther_inc;
binSize_diff = 0.1;
binSizePlot = 0.5;

figure
hold all
xLim = 2.5;
xtick = -xLim:binSizePlot:xLim;
xtickLabs = cellfun(@(x) num2str(x,'%.1f'),num2cell(xtick),'UniformOutput',false);
limPoint = max(abs(diff_vals));
if limPoint > xLim; error('Histogram xLim not correct'); end
bins = -limPoint:binSize_diff:limPoint;
binCenters = ( bins(1:end-1) + bins(2:end) ) / 2 ;
h = histcounts(diff_vals,bins);
bar(binCenters,h);
YLIM = get(gca,'YLim');
plot([0 0],YLIM,'k:','linewidth',1)
plot([med_VvsVO_diff med_VvsVO_diff],YLIM,'b:','linewidth',2)
set(gca,'XLim',[-xLim xLim],'Color','none','box','off','linewidth',1,...
    'tickdir','out','fontsize',14,'xtick',xtick,'xticklabel',xtickLabs)
xlabel('Mean (\DeltaFR_o_t_h_e_r - \DeltaFR_r_a_n_d_o_m)')
ylabel('Number of iterations (1000 total)')
title('Cells with increased V onset FRs')

saveFig(fullfile(saveFigFol,'histogram_diff'),'-dpdf')




