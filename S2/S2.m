close all

saveFigFol = '\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3\S2\saveFig';

cd('\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3\S2');

miceToUse = true(length(mouseIDforEachRec_VO),1);
% miceToUse( contains(mouseIDforEachRec_VO,'TRO49') ...
%     | contains(mouseIDforEachRec_VO,'TRO51')) = false;

recNum_perRow = 19; % example recording number

recInd = ~avoidRec_VO_AvsV & miceToUse;

%% Plot - stats - % of trials above threshold

lims = [ 0 100 ]; % min and max limits for plot

figure
hold all
xTickLab = [ {sprintf('%s%s','\leq',num2str(lims(1)))} {sprintf('%s%s','\geq',num2str(lims(2)))} ];

% Get data
perWin_CMN_VO_lim = perOver_CMN_noABS_morePCs;
perWin_off_VO_lim = perOver_offset_noABS_morePCs;

% Statistical test
pVal_per = signrank(perWin_CMN_VO_lim(recInd),perWin_off_VO_lim(recInd));

% Bind data within min and max limits
perWin_CMN_VO_lim(perWin_CMN_VO_lim < lims(1)) = lims(1);
perWin_CMN_VO_lim(perWin_CMN_VO_lim > lims(2)) = lims(2);
perWin_off_VO_lim(perWin_off_VO_lim < lims(1)) = lims(1);
perWin_off_VO_lim(perWin_off_VO_lim > lims(2)) = lims(2);

% Index to example recording
recNumInd = false(length(perWin_off_VO_lim),1);
recNumInd(recNum_perRow) = true;

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

mean_V = mean(perOver_CMN_noABS_morePCs(recInd));
SD_V = std(perOver_CMN_noABS_morePCs(recInd));
fprintf('\nV mean +- SD: %s +- %s',num2str(mean_V),num2str(SD_V))

mean_A = mean(perOver_offset_noABS_morePCs(recInd));
SD_A = std(perOver_offset_noABS_morePCs(recInd));
fprintf('\nA mean +- SD: %s +- %s',num2str(mean_A),num2str(SD_A))

% saveFig(fullfile(saveFigFol,'scatter_per'),'-dpdf');
