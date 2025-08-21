close all; clc

cd('\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3\S9')
saveFigFol = '\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3\S9\figSave';

%% Inputs

miceToUse = true(length(mouseIDforEachRec_VO),1);

recInd = ~avoidRec_VO_AvsV & miceToUse;
recIndPerCell = repelem(recInd,numCellPerRec_VO);

miceToUse = true(length(mouseIDforEachRec_VO),1);
% miceToUse( contains(mouseIDforEachRec_VO,'TRO49') ...
%     | contains(mouseIDforEachRec_VO,'TRO51')) = false;
miceToUse_perCell = repelem(miceToUse,numCellPerRec_VO,1);

% Look at primary AC recordings only, or not
primarySecondaryInd = true(length(miceToUse_perCell),1);
% primarySecondaryInd = tunedRec_AC_perCell;

indRec = ~avoidRec_VO_FM_perCell & miceToUse_perCell & primarySecondaryInd;

%% Percentage of FM trials

allPerOver_CMN = [ perOver_CMN_noABS(recInd) perOver_CMN_ROI(recInd,:) perOver_CMN_ball(recInd) ];
allPerOver_off = [ perOver_offset_noABS(recInd) perOver_off_ROI(recInd,:) perOver_off_ball(recInd) ];

xCMN_points = 1:3:16;
xOff_points = xCMN_points + 1;

numR = size(allPerOver_CMN,1);
medSide = .3;

colScat = [ repmat([1 0 0],numR,1) ; repmat([0 0 1],numR,1) ];
    
figure('Position', [1817 1000 838 286])
hold all
for i = 1:size(allPerOver_CMN,2)
    x1 = ones(numR,1) * xCMN_points(i);
    x2 = ones(numR,1) * xOff_points(i);
    scatter([x1 ; x2],[allPerOver_CMN(:,i) ; allPerOver_off(:,i)],50, colScat )
    for j = 1:numR
        plot( [x1(j) x2(j)] , [allPerOver_CMN(j,i) ; allPerOver_off(j,i)] ,'k')
    end
    % Median lines
    medCMN = median(allPerOver_CMN(:,i));
    medoff = median(allPerOver_off(:,i));
    plot([x1(i)-medSide x1(i)+medSide],[medCMN medCMN],'r','linewidth',4);
    plot([x2(i)-medSide x2(i)+medSide],[medoff medoff],'b','linewidth',4);
end

set(gca,'Color','none','YLim',[0 100],'XLim',[0 18],'tickdir','out',...
    'fontsize',15,'xtick',xCMN_points+.5,'xtickLabel','','linewidth',1)

% saveFig(fullfile(saveFigFol,'perFMtrials_ROI'),'-dpdf');

%% Difference in visually evoked firing

minPer = 5;

curInd = visRespCell_AC & FRch_base_to_V > 0 & indRec;

perOver_CMN_noABS_perCell = repelem(perOver_CMN_noABS,numCellPerRec_VO);
perOver_CMN_ROI_perCell = repelem(perOver_CMN_ROI,numCellPerRec_VO,1);
perOver_CMN_ball_perCell = repelem(perOver_CMN_ball,numCellPerRec_VO);

allPerOver_CMN_perCell = [ perOver_CMN_noABS_perCell perOver_CMN_ROI_perCell perOver_CMN_ball_perCell ];
perOverThresh = allPerOver_CMN_perCell >= minPer;

% FR change for low and high movement trials - unbound
FRch_base_to_V_lowMvmt_temp = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_noMvmt, [-inf inf]);
FRch_base_to_V_highMvmt_temp = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_mvmt, [-inf inf]);
FRch_base_to_V_lowMvmt_ROI_temp = nan(size(FRch_base_to_V_lowMvmt_ROI));
FRch_base_to_V_highMvmt_ROI_temp = nan(size(FRch_base_to_V_highMvmt_ROI));
for i = 1:4
    FRch_base_to_V_lowMvmt_ROI_temp(:,i) = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_noMvmt_ROI, [-inf inf]);
    FRch_base_to_V_highMvmt_ROI_temp(:,i) = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_mvmt_ROI, [-inf inf]);
end
FRch_base_to_V_lowMvmt_ball_temp = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_noMvmt_ball, [-inf inf]);
FRch_base_to_V_highMvmt_ball_temp = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_mvmt_ball, [-inf inf]);

% Difference in change in FR for high mvmt trials minus low mvmt trials
allFRCh_lowMvmt = [ FRch_base_to_V_lowMvmt_temp FRch_base_to_V_lowMvmt_ROI_temp FRch_base_to_V_lowMvmt_ball_temp ];
allFRCh_highMvmt = [ FRch_base_to_V_highMvmt_temp FRch_base_to_V_highMvmt_ROI_temp FRch_base_to_V_highMvmt_ball_temp ];

allPVals_lowVsHighMvmt = [ pVals_lowVsHighMvmt pVals_lowVsHighMvmt_ROI pVals_lowVsHighMvmt_ball ];
allPVals_lowVsHighMvmt_FDR = allPVals_lowVsHighMvmt;
for i = 1:size(allPVals_lowVsHighMvmt,2)
    allPVals_lowVsHighMvmt_FDR(curInd,i) = mafdr(allPVals_lowVsHighMvmt_FDR(curInd,i),'BHFDR',true);
end

sigInd = allPVals_lowVsHighMvmt_FDR < 0.05;

numSigInd = sum( sigInd & perOverThresh & repmat(curInd,1,6) );
numToTake = sum( perOverThresh & repmat(curInd,1,6) );
sumSigInd = (numSigInd ./ numToTake) * 100

FRchangeLims_split = [-3 3];
yticks = [FRchangeLims_split(1) 0 FRchangeLims_split(2)];
yTicks = { sprintf('%s%s','\leq','-2.5') '0' sprintf('%s%s','\geq','2.5') };

FRchange_all = allFRCh_highMvmt - allFRCh_lowMvmt;

divi = 5;
jitter = (rand(numCell,6)/divi) - (1/(divi*2));

figure('Position', [1817 1000 838 286])
hold all

for i = 1:6

    curDiffInd = curInd & perOverThresh(:,i);

    curFRdiff = FRchange_all(:,i);
    curSig = sigInd(:,i);

    FRchangeDiff_lim = curFRdiff;
    FRchangeDiff_lim(FRchangeDiff_lim < FRchangeLims_split(1)) = FRchangeLims_split(1);
    FRchangeDiff_lim(FRchangeDiff_lim > FRchangeLims_split(2)) = FRchangeLims_split(2);

    curX = ones(numCell,1) * i + jitter(:,i);

    % BS - non-sig
    curPlotInd = curDiffInd & ~curSig & VO.BSunits;
    scatter( curX(curPlotInd),  FRchangeDiff_lim(curPlotInd),...
        100,'filled',"^",'MarkerFaceColor','k',...
        'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
    % NS - non-sig
    curPlotInd = curDiffInd & ~curSig & ~VO.BSunits;
    scatter( curX(curPlotInd),  FRchangeDiff_lim(curPlotInd),...
        100,'filled','MarkerFaceColor','k',...
        'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
    % BS - sig
    curPlotInd = curDiffInd & curSig & VO.BSunits;
    scatter( curX(curPlotInd),  FRchangeDiff_lim(curPlotInd),...
        100,'filled',"^",'MarkerFaceColor','r',...
        'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
    % NS - sig
    curPlotInd = curDiffInd & curSig & ~VO.BSunits;
    scatter( curX(curPlotInd),  FRchangeDiff_lim(curPlotInd),...
        100,'filled','MarkerFaceColor','r',...
        'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
end

plot([.5 6.5],[0 0],'k:','linewidth',1)

set(gca,'XLim',[.5 6.5],'YLim',FRchangeLims_split,'fontsize',15,...
    'xtick',1:6,'tickdir','out','linewidth',1,...
    'Color','none','ytick',yticks,'yticklabel',yTicks);

% saveFig(fullfile(saveFigFol,'scatter_lowVsHighMvmt_ROI'),'-dpdf');

%% Example traces

exRec = 10;
timePeriod = [480 540];

startTime = VO.store_frameTimes{exRec}(1);

figure('Position',[2218 73 1046 1278])

timeAx = VO.store_frameTimes{exRec};
timeInd = timeAx > startTime+timePeriod(1) & timeAx < startTime+timePeriod(2);
% FM - all
sub{1} = subplot(6,1,1);
curData = VO.store_FM_PC1{exRec};
zScoreData = (curData - mean(curData)) / std(curData);
plot(timeAx(timeInd),zScoreData(timeInd))
% FM - ROI 1
sub{2} = subplot(6,1,2);
curData = VO.store_FM_ROI{exRec}(:,1);
zScoreData = (curData - mean(curData)) / std(curData);
plot(timeAx(timeInd),zScoreData(timeInd))
% FM - ROI 2
sub{3} = subplot(6,1,3);
curData = VO.store_FM_ROI{exRec}(:,2);
zScoreData = (curData - mean(curData)) / std(curData);
plot(timeAx(timeInd),zScoreData(timeInd))
% FM - ROI 3
sub{4} = subplot(6,1,4);
curData = VO.store_FM_ROI{exRec}(:,3);
zScoreData = (curData - mean(curData)) / std(curData);
plot(timeAx(timeInd),zScoreData(timeInd))
% FM - ROI 4
sub{5} = subplot(6,1,5);
curData = VO.store_FM_ROI{exRec}(:,4);
zScoreData = (curData - mean(curData)) / std(curData);
plot(timeAx(timeInd),zScoreData(timeInd))
% Ball movement
sub{6} = subplot(6,1,6);
timeAx = VO.store_ball_time{exRec};
timeInd = timeAx > startTime+timePeriod(1) & timeAx < startTime+timePeriod(2);
curData = VO.store_ball_vel{exRec};
zScoreData = (curData - mean(curData)) / std(curData);
plot(timeAx(timeInd),zScoreData(timeInd))

for i = 1:6
    set(sub{i},'XLim',[startTime+timePeriod(1) startTime+timePeriod(2)],...
        'XTick',[startTime+timePeriod(1) startTime+timePeriod(2)],...
        'xticklabel',[0 diff(timePeriod)],'Color','none','box','off',...
        'tickdir','out')
end

indVOn = find(VO.store_trigOn_VO{exRec} > startTime+timePeriod(1) & VO.store_trigOn_VO{exRec} < startTime+timePeriod(2));
indAOn = find(VO.store_trigOn_off{exRec} > startTime+timePeriod(1) & VO.store_trigOn_off{exRec} < startTime+timePeriod(2));

YLim = get(sub{6},'YLim');
yFill = [ YLim(1) YLim(1) YLim(2) YLim(2) YLim(1) ];
hold all
for i = 1:length(indVOn)
    curTrigOn = VO.store_trigOn_VO{exRec}(indVOn(i));
    xFill = [ curTrigOn curTrigOn+.5 curTrigOn+.5 curTrigOn curTrigOn ];
    fill(xFill,yFill,'r','FaceAlpha',.2,'EdgeAlpha',0)
end
for i = 1:length(indAOn)
    curTrigOn = VO.store_trigOn_off{exRec}(indAOn(i));
    xFill = [ curTrigOn curTrigOn+.5 curTrigOn+.5 curTrigOn curTrigOn ];
    fill(xFill,yFill,'b','FaceAlpha',.2,'EdgeAlpha',0)
end

% saveFig(fullfile(saveFigFol,'exampleTraces_mvmt'),'-dpdf');
