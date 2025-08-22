clear all; close all; clc

tic

cd('\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\Figures\ver3')

usePreWorkspace = false;

%% Load data - either load the data and process it (entire script below), or load a previous workspace with everything already processed

if usePreWorkspace
%     load('\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\workSpaceSave\allVars_2025-02-10_11-51-58.mat'); first submission - (saving because the regression to the mean analysis has a random component, so storing these just to go back to the exact version, if need-be)
%     load('\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\workSpaceSave\allVars_2025-05-29_14-27-31'); % revised submission
    load('\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\workSpaceSave\allVars_2025-08-19_18-15-54.mat'); % final submission
else
    VOdataFileName = 'VOstats_NP_2025-08-20_15-51-13.mat';
    AVOdataFileName = 'AVOstats_probeNum2_2025-05-15_11-51-18.mat';
    VO = load(fullfile('\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\code\VO_CMNonly\summaryData',VOdataFileName));
    AVO = load(fullfile('\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\code\AVO_collectStats',AVOdataFileName)); % the data
    load(fullfile('\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\code\VO_CMNonly\summaryData','FRchange_VandVO_TRO44.mat'),'TRO44')
end

%% Inputs

numPCsToUse = 1; % use this many FM PCs for analysis
numPCsToUse_higher = 3; % use this many FM PCs for more stringent FM analysis
zScoreThresh = 2; % threshold for detecting a FM
zScoreThresh_23 = 4; % a different threshold - for testing
alphaLevel = 0.01;
FRchangeLims = [-1 5]; % FR change limits on all plots
depthBins = [0 .36 .54 1]; % from Ryan's eLife paper
fracDepthLims = [ 0 1 ]; % fractional depth limits - only look at cells within these depth limits

%% Indices to time periods - AVO dataset

% AVO time periods
plotWin_AVO = [-1.7 1.5]; % face movement plotting window limits
xTicks_AVO = -1:1:1; % face movement plotting window tick points
baseWin = [-1.5 -0.5]; % baseline window
CMNwin = [-0.5 0]; % CMN window
RDSwin = [0 0.5]; % RDS window

% Starting data
tAx_AVO = AVO.binCenters_faceMvmtsAll;
% Time periods
baseInd_AVO = tAx_AVO > baseWin(1) & tAx_AVO < baseWin(2);
CMNactWin_AVO = tAx_AVO > CMNwin(1) & tAx_AVO < CMNwin(2);
RDSactWin_AVO = tAx_AVO > RDSwin(1) & tAx_AVO < RDSwin(2);
% Plotting window
plotWinInd_AVO = tAx_AVO > plotWin_AVO(1) & tAx_AVO < plotWin_AVO(2);
plotWinTime_AVO = tAx_AVO(plotWinInd_AVO);

%% Indices to time periods - VO dataset

% AVO time periods
plotWin_VO = [-1.0 1.5]; % face movement plotting window limits
xTicks_VO = -1:1:1; % face movement plotting window tick points
baseWin_VO = [-.5 0]; % baseline window
CMNwin_VO = [0 .5]; % CMN window

% Starting data
tAx_VO = VO.binCenters_faceMvmts;
% Time periods
baseInd_VO = tAx_VO > baseWin_VO(1) & tAx_VO < baseWin_VO(2);
CMNactWin_VO = tAx_VO > CMNwin_VO(1) & tAx_VO < CMNwin_VO(2);
% Plotting window
plotWinInd_VO = tAx_VO > plotWin_VO(1) & tAx_VO < plotWin_VO(2);
tToPlot_VO = tAx_VO(plotWinInd_VO);

%% Indices to time periods - VO - noise bursts

% VO time periods - noise bursts
plotWin_off = [-2.5 2.5]; % face movement plotting window limits
xTicks_off = -2:1:2; % face movement plotting window tick points
baseWin_off = [-1.05 0]; % baseline window
actWin_off = [0 1.05]; % CMN window

% Starting data
tAx_off = VO.binCenters_off;
% Time periods
baseInd_off = tAx_off > baseWin_off(1) & tAx_off < baseWin_off(2);
actInd_off = tAx_off > actWin_off(1) & tAx_off < actWin_off(2);
% Plotting window
plotWinInd_off = tAx_off > plotWin_off(1) & tAx_off < plotWin_off(2);
tToPlot_off = tAx_off(plotWinInd_off);

%% General automatic inputs

% AVO dataset
numRecIn_AVO = size(AVO.trialMvmtsAll,1); % total number of recordings - AVO
numTrials_AVO = sum(AVO.AVtrials{end}); % number of AV trials
[~,recFirstPositions_AVO] = max(AVO.recNum == unique(AVO.recNum)',[],1); % index to the first cell in each recording
mouseIDforEachRec_AVO = cellstr(AVO.mouseID_col(recFirstPositions_AVO,:));

% VO dataset
[~,recFirstPositions_VO] = max(VO.recNumber == unique(VO.recNumber)',[],1); % index to the first cell in each recording
numRec_VO = length(recFirstPositions_VO); % total number of recordings - VO
VCrecIDs_VO = unique(VO.recNumber); % recording IDs
mouseIDforEachRec_VO = VO.mouseID(recFirstPositions_VO);
numCellPerRec_VO = sum(VO.recNumber == VCrecIDs_VO'); % number of cells for each recording
numCell = length(VO.FRs_Vonset); % number of cells, total

%% Primary/secondary designation - VO

latThreshold = 0.014; % latency threshold for discriminating between primary-like and secondary-like onset latencies
percentRequired = 30; % how many significant channels are required to be tuned for a site to be considered primary-like

toneSig = cellfun(@(x) x < alphaLevel, VO.tonePvalsTwoSided_MU, 'UniformOutput', false); % channels with significant responses
varSig = cellfun(@(x) x < alphaLevel, VO.toneVarFrac_MU, 'UniformOutput', false); % channels with significant FTC variances
median_latency = cellfun(@(x,y,z) median(x(y & z),'omitnan'),VO.toneLatencies_MU, toneSig,VO.toneDepthWithin_MU); % median latency for all significant channels, for all channels in the cortex - one latency per recording
tuned_channels = cellfun(@(x,y,z) x & y & z, varSig , toneSig , VO.toneDepthWithin_MU , 'UniformOutput', false ); % channels with significant responses and tuning tests, for all channels in the cortex
tuned_recPer = cellfun(@(x,y,z) (sum(x) / sum(y & z)) * 100 , tuned_channels , toneSig , VO.toneDepthWithin_MU); % percent of significant channels with significant tuning tests

tunedRec = median_latency <= latThreshold & tuned_recPer >= percentRequired; % recordings with primary-like latencies and significantly tuned
tunedOnlyRec = tuned_recPer >= percentRequired; 

tunedRec_AC = tunedRec(1:2:end); % tunedRec, but now for AC recordings only (VC excluded)
tunedRec_AC_perCell = repelem(tunedRec_AC,numCellPerRec_VO); % the same as above, but now one logical value for each cell

medLatency_byCell = repelem(median_latency(1:2:end),numCellPerRec_VO);

% Median latency plotting for tuned and untuned recording sites
cols = zeros(length(median_latency),3);
cols(tunedOnlyRec,:) = repmat([1 0 0],sum(tunedOnlyRec),1);
figure
indTake = 1:2:length(median_latency);
scatter(1:length(median_latency)/2,median_latency(indTake),50,cols(indTake,:),'filled')

%% Recordings to avoid - VO dataset

% VO recordings to avoid for A versus V analysis - figure 1
avoidRec_VO_AvsV = false(numRec_VO,1);
notEnoughTrials_CMN = cellfun(@(x) size(x,1),VO.trialMvmtAll(recFirstPositions_VO)) < 199; % recordings with less than 199 CMN trials
notEnoughTrials_noiseBurst = cellfun(@(x) size(x,1),VO.trialMvmtAll_offset(recFirstPositions_VO)) < 30; % recordings with less than 30 noise burst trials
containsNaNs = cellfun(@(x) sum(isnan(x(1:end,:,:)),'all'),VO.trialMvmtAll(recFirstPositions_VO)) > 0;
avoidRec_VO_AvsV( notEnoughTrials_CMN | notEnoughTrials_noiseBurst | containsNaNs ) = true;
% VO recordings to avoid for w/ and w/o FM - figure 2 onwards
avoidRec_VO_FM = false(numRec_VO,1);
avoidRec_VO_FM( notEnoughTrials_CMN | containsNaNs ) = true;
avoidRec_VO_FM_perCell = repelem(avoidRec_VO_FM,numCellPerRec_VO,1);

%% Recordings to avoid - AVO dataset

avoidRec_AVO = false(numRecIn_AVO,1);
containsNaNs = cellfun(@(x) sum(isnan(x(:))),AVO.trialMvmtsAll) > 0;
isemptyData = cellfun(@isempty ,AVO.trialMvmtsAll) > 0;
avoidRec_AVO(containsNaNs | isemptyData) = true;

%% Color maps

colMap_MG = colMapGen([1 0 1],[0 1 0],30,'midCol',[1 1 1]);
colMap_MG_exp = colMapGenExp([0 1 0],[1 1 1],[1 0 1],15);
colMap_RB = colMapGen([1 0 0],[0 0 1],30,'midCol',[1 1 1]);
colMap_RB_exp = colMapGenExp([0 0 1],[1 1 1],[1 0 0],15);
colMap_RW = [ ones(10,1) linspace(1,0,10)' linspace(1,0,10)' ];
colMap_WM_exp = colMapGenExp([1 1 1],[1 1 1],[1 0 1],16);
colMap_WM_exp = colMap_WM_exp(16:end,:);
numSegments = 15;
colMap_WM = [ones(numSegments,1) linspace(1,0,numSegments)' ones(numSegments,1)];
colMean = [0.68 0.14 1];

%% Normalized PSTHs

% V
baseInd = VO.binCenters_VO > baseWin_VO(1) & VO.binCenters_VO < baseWin_VO(2);
actInd = VO.binCenters_VO > CMNwin_VO(1) & VO.binCenters_VO < CMNwin_VO(2);
PSTH_Vtrial_bs = VO.PSTHout_V  - mean(VO.PSTHout_V(:,baseInd),2);
PSTH_Vtrial_zScore = PSTH_Vtrial_bs ./ std(PSTH_Vtrial_bs(:,baseInd),[],2);
PSTH_Vtrial_norm = PSTH_Vtrial_bs ./ max(abs(PSTH_Vtrial_bs(:,actInd)),[],2);

% Noise bursts
PSTH_offTrial_bs = VO.PSTHout_off - mean(VO.PSTHout_off(:,baseInd_off),2);
PSTH_offTrial_zScore = PSTH_offTrial_bs ./ std(PSTH_offTrial_bs(:,baseInd_off),[],2);
PSTH_offTrial_norm = PSTH_offTrial_bs ./ max(abs(PSTH_offTrial_bs(:,actInd_off)),[],2);

%% AVO - CMM & RDS - % of trials above threshold and FM amplitude

movementsStore_AVO = cell(numRecIn_AVO,1);
movementsStore_AVO_noABS = cell(numRecIn_AVO,1);
trialsOverThresh_AVO_noABS = cell(numRecIn_AVO,2);
perOver_AVO_noABS = nan(numRecIn_AVO,2);
maxInWin_AVO_noABS = cell(numRecIn_AVO,2);
meanInWin_AVO_noABS = cell(numRecIn_AVO,2);

timeWindows = [ CMNactWin_AVO ; RDSactWin_AVO ];

for i = 1:numRecIn_AVO

    AVtials = AVO.AVtrials{i};

    % Get the traces
    if ~isempty(AVO.trialMvmtsAll{i})
        
        curRecMvmts = AVO.trialMvmtsAll{i};
        curRecMvmts_bs = curRecMvmts - mean(curRecMvmts(:,baseInd_AVO,:),2);
        curRecMvmts_zScore = curRecMvmts_bs ./ std(curRecMvmts_bs(:,baseInd_AVO,:),[],2);
        movementsStore_AVO_noABS{i} = curRecMvmts_zScore; % store movement traces
        curRecMvmts_zScore_abs = abs(curRecMvmts_zScore); % z-score and abs
        movementsStore_AVO{i} = curRecMvmts_zScore_abs; % store movement traces 
    
        % CMN and RDS
        for ii = 1:2 % loop through the CMN and RDS windows
            
            curWin = timeWindows(ii,:);
            
            % Trials over threshold
            trialsOverThForEachPC = mean(curRecMvmts_zScore(:,curWin,:),2) > zScoreThresh;
            trialsOverThresh_AVO_noABS{i,ii} = trialsOverThForEachPC;
          
            % Get the percentage of trials over threshold for the curWin
            overThreshold = sum( trialsOverThForEachPC(:,:,1:numPCsToUse) , 3 ) > 0; % now find trials where any of the PCs from 1:numPCsToUse within curWin are above threshold
            perOver_AVO_noABS(i,ii) = ( sum(overThreshold(AVtials)) / length(overThreshold(AVtials)) ) * 100; % percentage of trials over threshold
            
            % Get the maximal response for each trial
            maxInWin_AVO_noABS{i,ii} = max(curRecMvmts_zScore(:,curWin,:),[],2);
        
            % Get the mean response for each trial
            meanInWin_AVO_noABS{i,ii} = mean(curRecMvmts_zScore(:,curWin,:),2);

        end
        
    end

end

%% VO - CMM - % of trials above threshold and FM amplitude

movementsStore_VO_CMN = cell(numRec_VO,1);
movementsStore_VO_CMN_noABS = cell(numRec_VO,1);
perOver_CMN = nan(numRec_VO,1);
perOver_CMN_morePCs = nan(numRec_VO,1);
trialsOverThresh_CMN_VO = cell(numRec_VO,1);
trialsOverThresh_CMN_VO_PC23 = cell(numRec_VO,1);
maxInWin_VO = cell(numRec_VO,1);
meanInWin_VO = cell(numRec_VO,1);
perOver_CMN_noABS = nan(numRec_VO,1);
perOver_CMN_noABS_morePCs = nan(numRec_VO,1);
trialsOverThresh_CMN_VO_noABS = cell(numRec_VO,1);
maxInWin_VO_noABS = cell(numRec_VO,1);
meanInWin_VO_noABS = cell(numRec_VO,1);
trialsOverThresh_CMN_VO_neg = cell(numRec_VO,1);

for i = 1:numRec_VO

    optoTrials = VO.FRs_trialType_optoTrials{recFirstPositions_VO(i)}; % find opto trials

    % Get the traces
    curRecMvmts = VO.trialMvmtAll{recFirstPositions_VO(i)}; % CMN movement data for current recording
    curRecMvmts_bs = curRecMvmts - mean(curRecMvmts(:,baseInd_VO,:),2); % baseline subtract
    curRecMvmts_zScore = curRecMvmts_bs ./ std(curRecMvmts_bs(:,baseInd_VO,:),[],2);
    movementsStore_VO_CMN_noABS{i} = curRecMvmts_zScore; % store movement traces
    curRecMvmts_zScore_abs = abs( curRecMvmts_zScore ); % z-score and abs
    movementsStore_VO_CMN{i} = curRecMvmts_zScore_abs; % store movement traces

    % Trials over threshold
    trialsOverThresh_CMN_VO_noABS{i} = mean(curRecMvmts_zScore(:,CMNactWin_VO,:),2) > zScoreThresh;
    
    % Max in CMN window 
    maxInWin_VO_noABS{i} = max(curRecMvmts_zScore(:,CMNactWin_VO,:),[],2);

    % Mean in CMN window 
    meanInWin_VO_noABS{i} = mean(curRecMvmts_zScore(:,CMNactWin_VO,:),2);

    % Percent over threshold, for numPCsToUse and numPCsToUse_higher
    overThreshold = sum( trialsOverThresh_CMN_VO_noABS{i}(:,:,1:numPCsToUse) , 3 ) > 0;
    perOver_CMN_noABS(i) = ( sum(overThreshold(~optoTrials)) / length(overThreshold(~optoTrials)) ) * 100;
    overThreshold = sum( trialsOverThresh_CMN_VO_noABS{i}(:,:,1:numPCsToUse_higher) , 3 ) > 0;
    perOver_CMN_noABS_morePCs(i) = ( sum(overThreshold(~optoTrials)) / length(overThreshold(~optoTrials)) ) * 100;

    %%%%%%%%%% abs stuff %%%%%%%%%%
    
    % Trials over threshold
    trialsOverThresh_CMN_VO{i} = mean(curRecMvmts_zScore_abs(:,CMNactWin_VO,:),2) > zScoreThresh;
    trialsOverThresh_CMN_VO_PC23{i} = mean(curRecMvmts_zScore_abs(:,CMNactWin_VO,:),2) > zScoreThresh_23;

    % Max in CMN window
    maxInWin_VO{i} = max(curRecMvmts_zScore_abs(:,CMNactWin_VO,:),[],2);

    % Max in CMN window
    meanInWin_VO{i} = mean(curRecMvmts_zScore_abs(:,CMNactWin_VO,:),2);
    
    % Percent over threshold, for numPCsToUse and numPCsToUse_higher
    overThreshold = sum( trialsOverThresh_CMN_VO{i}(:,:,1:numPCsToUse) , 3 ) > 0;
    perOver_CMN(i) = ( sum(overThreshold(~optoTrials)) / length(overThreshold(~optoTrials)) ) * 100; 
    % Get the percentage of trials over threshold - for numPCsToUse_higher
    overThreshold_morePCs = sum( trialsOverThresh_CMN_VO{i}(:,:,1:numPCsToUse_higher) , 3 ) > 0; 
    perOver_CMN_morePCs(i) = ( sum(overThreshold_morePCs(~optoTrials)) / length(overThreshold_morePCs(~optoTrials)) ) * 100; 

    % Neg crossings
    trialsOverThresh_CMN_VO_neg{i} = sum ( -curRecMvmts_zScore(:,CMNactWin_VO,:) > zScoreThresh , 2 ) > 0 ; 

end

%% VO - noise burst - % of trials above threshold and FM amplitude

movementsStore_VO_offset = cell(numRec_VO,1); 
movementsStore_VO_offset_noABS = cell(numRec_VO,1);
perOver_offset = nan(numRec_VO,1);
perOver_offset_morePCs = nan(numRec_VO,1);
trialsOverThresh_offset_VO = cell(numRec_VO,1);
trialsOverThresh_offset_VO_PC23 = cell(numRec_VO,1);
maxInWin_VO_offset = cell(numRec_VO,1);
meanInWin_VO_offset = cell(numRec_VO,1);
perOver_offset_noABS = nan(numRec_VO,1);
perOver_offset_noABS_morePCs = nan(numRec_VO,1);
trialsOverThresh_offset_VO_noABS = cell(numRec_VO,1);
maxInWin_VO_offset_noABS = cell(numRec_VO,1);
meanInWin_VO_offset_noABS = cell(numRec_VO,1);

for i = 1:numRec_VO

    % Get the traces
    curRecMvmts = VO.trialMvmtAll_offset{recFirstPositions_VO(i)}; % noise burst movement data for current recording
    curRecMvmts_bs = curRecMvmts - mean(curRecMvmts(:,baseInd_VO,:),2); % baseline subtract
    curRecMvmts_zScore = curRecMvmts_bs ./ std(curRecMvmts(:,baseInd_VO,:),[],2);
    movementsStore_VO_offset_noABS{i} = curRecMvmts_zScore;
    curRecMvmts_zScore_abs = abs( curRecMvmts_zScore );
    movementsStore_VO_offset{i} = curRecMvmts_zScore_abs;

    % Trials over threshold
    trialsOverThresh_offset_VO_noABS{i} = mean(curRecMvmts_zScore(:,CMNactWin_VO,:),2) > zScoreThresh;
    
    % Max in CMN window 
    maxInWin_VO_offset_noABS{i} = max(curRecMvmts_zScore(:,CMNactWin_VO,:),[],2);

    % Mean in CMN window 
    meanInWin_VO_offset_noABS{i} = mean(curRecMvmts_zScore(:,CMNactWin_VO,:),2);

    % Percent over threshold, for numPCsToUse and numPCsToUse_higher
    overThreshold = sum( trialsOverThresh_offset_VO_noABS{i}(:,:,1:numPCsToUse) , 3 ) > 0; 
    perOver_offset_noABS(i) = ( sum(overThreshold) / length(overThreshold) ) * 100;
    overThreshold = sum( trialsOverThresh_offset_VO_noABS{i}(:,:,1:numPCsToUse_higher) , 3 ) > 0; 
    perOver_offset_noABS_morePCs(i) = ( sum(overThreshold) / length(overThreshold) ) * 100;

    %%%%%%%%%% abs stuff %%%%%%%%%%

    % Trials over threshold
    trialsOverThresh_offset_VO{i} = mean(curRecMvmts_zScore_abs(:,CMNactWin_VO,:),2) > zScoreThresh;
    trialsOverThresh_offset_VO_PC23{i} = mean(curRecMvmts_zScore_abs(:,CMNactWin_VO,:),2) > zScoreThresh_23;

    % Max in window
    maxInWin_VO_offset{i} = max(curRecMvmts_zScore_abs(:,CMNactWin_VO,:),[],2);

    % Mean in window
    meanInWin_VO_offset{i} = mean(curRecMvmts_zScore_abs(:,CMNactWin_VO,:),2);

    % Percent over threshold, for numPCsToUse and numPCsToUse_higher
    overThreshold = sum( trialsOverThresh_offset_VO{i}(:,:,1:numPCsToUse) , 3 ) > 0;
    perOver_offset(i) = ( sum(overThreshold) / length(overThreshold) ) * 100;
    overThreshold_morePCs = sum( trialsOverThresh_offset_VO{i}(:,:,1:numPCsToUse_higher) , 3 ) > 0;
    perOver_offset_morePCs(i) = ( sum(overThreshold_morePCs) / length(overThreshold_morePCs) ) * 100;

end

%% Binned counts of single-trial mean amplitudes - VO

curMvmtData_CMN = movementsStore_VO_CMN_noABS;
curMvmtData_offset = movementsStore_VO_offset_noABS;

% Lims and bins for 'mean' response binning
lims_thSel_mean = [-6 6];
binSize = 0.5;
bins_thSel_mean = lims_thSel_mean(1):binSize:lims_thSel_mean(2);
numBins_mean = length(bins_thSel_mean)-1;

totalNumPCs = size(VO.trialMvmtAll{1},3);

h_CMN_mean = nan(numRec_VO,numBins_mean,totalNumPCs);
h_off_mean = nan(numRec_VO,numBins_mean,totalNumPCs);
for i = 1:numRec_VO

    % CMN
    numTrialsCur = size(curMvmtData_CMN{i},1);
    optoTrials = VO.FRs_trialType_optoTrials{recFirstPositions_VO(i)}; % find opto trials
    % Mean
    % Get mean responses in the stimulus window for each trial
    meanResp_curRec = reshape( mean( curMvmtData_CMN{i}(:,CMNactWin_VO,:),2) , [numTrialsCur totalNumPCs] );
    % Bind data within limits (lims_thSel)
    meanResp_curRec_lims = meanResp_curRec;
    meanResp_curRec_lims(meanResp_curRec_lims < lims_thSel_mean(1)) = lims_thSel_mean(1);
    meanResp_curRec_lims(meanResp_curRec_lims > lims_thSel_mean(2)) = lims_thSel_mean(2);
    % Get counts of maximum amplitudes for each amplitude bin
    for j = 1:totalNumPCs
        h_CMN_mean(i,:,j) = histcounts(meanResp_curRec_lims(~optoTrials,j),bins_thSel_mean);
    end

    % Noise burst
    numTrialsCur = size(curMvmtData_offset{i},1);
    % Mean
    % Get maximum responses for each trial
    meanResp_curRec = reshape( mean( curMvmtData_offset{i}(:,CMNactWin_VO,:),2) , [numTrialsCur totalNumPCs] );
    % Bind data within limits (lims_thSel)
    meanResp_curRec_lims = meanResp_curRec;
    meanResp_curRec_lims(meanResp_curRec_lims < lims_thSel_mean(1)) = lims_thSel_mean(1);
    meanResp_curRec_lims(meanResp_curRec_lims > lims_thSel_mean(2)) = lims_thSel_mean(2);
    % Get counts of maximum amplitudes for each amplitude bin
    for j = 1:totalNumPCs
        h_off_mean(i,:,j) = histcounts(meanResp_curRec_lims(:,j),bins_thSel_mean);
    end

end

%% Noise burst VC spiking stuff 

fil_VC = VO.fracDepths >= fracDepthLims(1) & VO.fracDepths <= fracDepthLims(2) & ~avoidRec_VO_FM_perCell; % intitial SU filter

numUnits = length(fil_VC);

% Noise burst evoked p-values for all cells
pVal_offResp_CMNDur = nan(numUnits,1);
for i = 1:numUnits
    pVal_offResp_CMNDur(i) = signrank(VO.FRs_off_baseline_CMNlen{i},VO.FRs_off_act_CMNlen{i});
end
% FDR correct p values
pVal_offResp_FDR_CMNDur = nan(numUnits,1);
pVal_offResp_FDR_CMNDur(fil_VC & VO.probeNum == 1) = mafdr(pVal_offResp_CMNDur(fil_VC & VO.probeNum == 1),'BHFDR',true);
pVal_offResp_FDR_CMNDur(fil_VC & VO.probeNum == 2) = mafdr(pVal_offResp_CMNDur(fil_VC & VO.probeNum == 2),'BHFDR',true);

% Noise burst and non-noise burst responsive cells
noiseRespCell_VC_CMNdur = fil_VC & VO.probeNum == 2 & pVal_offResp_FDR_CMNDur < alphaLevel & VO.movieStatus < 3;
nonNoiseRespCell_VC_CMNdur = fil_VC & VO.probeNum == 2 & pVal_offResp_FDR_CMNDur >= alphaLevel & VO.movieStatus < 3;

% Index of trials to take
numOffTrials = cellfun(@(x) size(x,1), VO.FRs_off_baseline_s);
indexTake = arrayfun(@(n) true(n, 1), numOffTrials, 'UniformOutput', false);

% FR change for noise burst
FRch_base_to_noise_CMNlen = getFRchangeVals_log2(VO.FRs_off_baseline_CMNlen, VO.FRs_off_act_CMNlen, indexTake, FRchangeLims);

% FR change for noise burst - Z-scores
FRzScore_noiseResp_CMNlen = getFRzScore(VO.FRs_off_baseline_CMNlen, VO.FRs_off_act_CMNlen, indexTake);

% Index to trials with movement (all PCs), for each cell
trialsOverThresh_off_VO_perCell = repelem(trialsOverThresh_offset_VO_noABS, numCellPerRec_VO, 1);

% Index to trials w/ & w/o movement - numPCsToUse
index_noMvmt_off = cellfun(@(x,y) sum(x(:,:,1:numPCsToUse),3) == 0 & y , trialsOverThresh_off_VO_perCell , indexTake , 'UniformOutput', false );
index_mvmt_off = cellfun(@(x,y) sum(x(:,:,1:numPCsToUse),3) > 0  & y , trialsOverThresh_off_VO_perCell , indexTake , 'UniformOutput',  false );

% FR change for low and high movement trials
FRch_base_to_off_lowMvmt_CMNdur = getFRchangeVals_log2(VO.FRs_off_baseline_CMNlen, VO.FRs_off_act_CMNlen, index_noMvmt_off, [-inf inf]);
FRch_base_to_off_highMvmt_CMNdur = getFRchangeVals_log2(VO.FRs_off_baseline_CMNlen,VO.FRs_off_act_CMNlen,index_mvmt_off, [-inf inf]);

%% VO - spiking stuff

fil = VO.fracDepths >= fracDepthLims(1) & VO.fracDepths <= fracDepthLims(2) & ~avoidRec_VO_FM_perCell; % intitial SU filter

numUnits = length(fil);

% Visually evoked p-values for all cells
pVal_visResp = nan(numUnits,1);
for i = 1:numUnits
    curOptoTrials = VO.FRs_trialType_optoTrials{i};
    pVal_visResp(i) = signrank(VO.FRs_baseline{i}(~curOptoTrials), VO.FRs_Vonset{i}(~curOptoTrials));
end
% FDR correct p values - correct the VC and AC separately, and do not include cells out of the responsive span
pVal_visResp_FDR = nan(numUnits,1);
pVal_visResp_FDR(fil & VO.probeNum == 1) = mafdr(pVal_visResp(fil & VO.probeNum == 1),'BHFDR',true);
pVal_visResp_FDR(fil & VO.probeNum == 2) = mafdr(pVal_visResp(fil & VO.probeNum == 2),'BHFDR',true);

% Visually and non-visually responsive cells
visRespCell_AC = fil & VO.probeNum == 1 & pVal_visResp_FDR < alphaLevel & VO.movieStatus < 3;
nonVisRespCell_AC = fil & VO.probeNum == 1 & pVal_visResp_FDR >= alphaLevel & VO.movieStatus < 3;
visRespCell_VC = fil & VO.probeNum == 2 & pVal_visResp_FDR < alphaLevel & VO.movieStatus < 3;

% No opto trials
indexToNonOptoCells = cellfun(@(x) ~x, VO.FRs_trialType_optoTrials,  'UniformOutput',  false);
% Opto trials
indexToOptoCells = cellfun(@(x) x, VO.FRs_trialType_optoTrials,  'UniformOutput',  false);

% FR change for CMN
FRch_base_to_V = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, indexToNonOptoCells, FRchangeLims);
FRch_base_to_V_opto = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, indexToOptoCells, FRchangeLims);

% FR change - Z-scores
FRzScore_Vresp = getFRzScore(VO.FRs_baseline, VO.FRs_Vonset, indexToNonOptoCells);

% Index to trials with movement (all PCs), for each cell
trialsOverThresh_CMN_VO_perCell = repelem(trialsOverThresh_CMN_VO_noABS, numCellPerRec_VO, 1);

% Index to trials with no opto, and w/ & w/o movement - numPCsToUse
index_noMvmt = cellfun(@(x) sum(x(:,:,1:numPCsToUse),3) == 0  , trialsOverThresh_CMN_VO_perCell , 'UniformOutput',  false );
index_noOpto_noMvmt = cellfun(@(x,y) ~x & sum(y(:,:,1:numPCsToUse),3) == 0  , ...
    VO.FRs_trialType_optoTrials, trialsOverThresh_CMN_VO_perCell , ...
    'UniformOutput',  false );
index_noOpto_mvmt = cellfun(@(x,y) ~x & sum(y(:,:,1:numPCsToUse),3) > 0  , ...
    VO.FRs_trialType_optoTrials, trialsOverThresh_CMN_VO_perCell , ...
    'UniformOutput',  false );
% Index to trials with no opto, and w/ & w/o movement - numPCsToUse_higher
index_noMvmt_morePCs = cellfun(@(x) sum(x(:,:,1:numPCsToUse_higher),3) == 0  , trialsOverThresh_CMN_VO_perCell , 'UniformOutput',  false );
index_noOpto_noMvmt_morePCs = cellfun(@(x,y) ~x & sum(y(:,:,1:numPCsToUse_higher),3) == 0  , ...
    VO.FRs_trialType_optoTrials, trialsOverThresh_CMN_VO_perCell , ...
    'UniformOutput',  false );
index_noOpto_mvmt_morePCs = cellfun(@(x,y) ~x & sum(y(:,:,1:numPCsToUse_higher),3) > 0  , ...
    VO.FRs_trialType_optoTrials, trialsOverThresh_CMN_VO_perCell , ...
    'UniformOutput',  false );

% FR change for low and high movement trials
FRch_base_to_V_lowMvmt = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_noMvmt, FRchangeLims);
FRch_base_to_V_highMvmt = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_mvmt, FRchangeLims);
% FR change for low and high movement trials - more PCs
FRch_base_to_V_lowMvmt_morePCs = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_noMvmt_morePCs, FRchangeLims);
FRch_base_to_V_highMvmt_morePCs = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_mvmt_morePCs, FRchangeLims);

% Units with significant movement modulation 
[pVals_lowVsHighMvmt, pVals_lowVsHighMvmt_perm] = getLowVsHighMvmtPvals(VO.FRs_baseline, VO.FRs_Vonset, VO.FRs_trialType_optoTrials, index_noMvmt);
% Units with significant movement modulation - more PCs
[pVals_lowVsHighMvmt_morePCs, pVals_lowVsHighMvmt_morePCs_perm] = getLowVsHighMvmtPvals(VO.FRs_baseline, VO.FRs_Vonset, VO.FRs_trialType_optoTrials, index_noMvmt_morePCs);

% Mean FRs
funcMean = @(x,y) mean(x(y));
meanFR_noOpto = cellfun(funcMean, VO.FRs_Vonset, indexToNonOptoCells);
meanFR_opto = cellfun(funcMean, VO.FRs_Vonset, indexToOptoCells);
meanFR_noOpto_baseline = cellfun(funcMean,  VO.FRs_baseline, indexToNonOptoCells);
meanFR_opto_baseline = cellfun(funcMean,VO.FRs_baseline, indexToOptoCells);

% abs data
trialsOverThresh_CMN_VO_perCell_abs = repelem(trialsOverThresh_CMN_VO, numCellPerRec_VO, 1);
index_noMvmt_abs = cellfun(@(x) sum(x(:,:,1:numPCsToUse),3) == 0  , trialsOverThresh_CMN_VO_perCell_abs , 'UniformOutput',  false );
index_noOpto_noMvmt_abs = cellfun(@(x,y) ~x & sum(y(:,:,1:numPCsToUse),3) == 0  , ...
    VO.FRs_trialType_optoTrials, trialsOverThresh_CMN_VO_perCell_abs , ...
    'UniformOutput',  false );
index_noOpto_mvmt_abs = cellfun(@(x,y) ~x & sum(y(:,:,1:numPCsToUse),3) > 0  , ...
    VO.FRs_trialType_optoTrials, trialsOverThresh_CMN_VO_perCell_abs , ...
    'UniformOutput',  false );
FRch_base_to_V_lowMvmt_abs = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_noMvmt_abs, FRchangeLims);
FRch_base_to_V_highMvmt_abs = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_mvmt_abs, FRchangeLims);
[pVals_lowVsHighMvmt_abs, pVals_lowVsHighMvmt_abs_perm] = getLowVsHighMvmtPvals(VO.FRs_baseline, VO.FRs_Vonset, VO.FRs_trialType_optoTrials, index_noMvmt_abs);

%% Bins for spike counts above baseline plots 

% Spike timing
bins_spkCounts = plotWin_VO(1):0.1:plotWin_VO(2);
binCenters_spkCounts = (bins_spkCounts(1:end-1) + bins_spkCounts(2:end)) / 2;
binBaseInd = binCenters_spkCounts > baseWin_VO(1) & binCenters_spkCounts < baseWin_VO(2);
binActCC = binCenters_spkCounts > CMNwin_VO(1) & binCenters_spkCounts < CMNwin_VO(2);

%% Binned spike counts and deltas - V

binnedSpkCounts_bs_store_V = cell(numUnits,1);
binnedSpkCounts_bs_mvmt_V = nan(numUnits,length(binCenters_spkCounts));
binnedSpkCounts_bs_noMvmt_V = nan(numUnits,length(binCenters_spkCounts));
delta_V = nan(numUnits,1);
delta_V_noABS = nan(numUnits,1);
for curUnit = 1:numUnits

    numTrials = sum(~VO.FRs_trialType_optoTrials{curUnit});
    curRastTimes = VO.raster_times_V{curUnit}; % spike times
    curRastTrInds = VO.raster_trInd_V{curUnit}; % spike time trials indices
    binnedSpkCounts = nan(numTrials,length(binCenters_spkCounts));
    for j = 1:numTrials
        binnedSpkCounts(j,:) = histcounts(curRastTimes(curRastTrInds==j), bins_spkCounts);
    end
    binnedSpkCounts_bs = binnedSpkCounts - mean(binnedSpkCounts(:,binBaseInd),2);
    binnedSpkCounts_bs_store_V{curUnit} = binnedSpkCounts_bs;

    optoTrials = VO.FRs_trialType_optoTrials{curUnit};
    binnedSpkCounts_bs_mvmt_V(curUnit,:) = mean(binnedSpkCounts_bs(trialsOverThresh_CMN_VO_perCell{curUnit}(~optoTrials,:,1),:));
    binnedSpkCounts_bs_noMvmt_V(curUnit,:) = mean(binnedSpkCounts_bs(~trialsOverThresh_CMN_VO_perCell{curUnit}(~optoTrials,:,1),:));
    
    delta_V(curUnit) = sum(abs(binnedSpkCounts_bs_mvmt_V(curUnit,binActCC) - binnedSpkCounts_bs_noMvmt_V(curUnit,binActCC)));
    delta_V_noABS(curUnit) = sum(binnedSpkCounts_bs_mvmt_V(curUnit,binActCC) - binnedSpkCounts_bs_noMvmt_V(curUnit,binActCC));

end

%% Binned spike counts and deltas - noise bursts

binnedSpkCounts_bs_store_off = cell(numUnits,1);
binnedSpkCounts_bs_mvmt_off = nan(numUnits,length(binCenters_spkCounts));
binnedSpkCounts_bs_noMvmt_off = nan(numUnits,length(binCenters_spkCounts));
delta_off = nan(numUnits,1);
delta_off_noABS = nan(numUnits,1);
for curUnit = 1:numUnits

    numTrials = size(VO.trialMvmtAll_offset{curUnit},1);
    curRastTimes = VO.raster_times_off{curUnit}; % spike times
    curRastTrInds = VO.raster_trInd_off{curUnit}; % spike time trials indices
    binnedSpkCounts = nan(numTrials,length(binCenters_spkCounts));
    for j = 1:numTrials
        binnedSpkCounts(j,:) = histcounts(curRastTimes(curRastTrInds==j),bins_spkCounts);
    end
    binnedSpkCounts_bs = binnedSpkCounts - mean(binnedSpkCounts(:,binBaseInd),2);
    binnedSpkCounts_bs_store_off{curUnit} = binnedSpkCounts_bs;

    binnedSpkCounts_bs_mvmt_off(curUnit,:) = mean(binnedSpkCounts_bs(trialsOverThresh_off_VO_perCell{curUnit}(:,:,1),:));
    binnedSpkCounts_bs_noMvmt_off(curUnit,:) = mean(binnedSpkCounts_bs(~trialsOverThresh_off_VO_perCell{curUnit}(:,:,1),:));
    
    delta_off(curUnit) = sum(abs(binnedSpkCounts_bs_mvmt_off(curUnit,binActCC) - binnedSpkCounts_bs_noMvmt_off(curUnit,binActCC)));
    delta_off_noABS(curUnit) = sum(binnedSpkCounts_bs_mvmt_off(curUnit,binActCC) - binnedSpkCounts_bs_noMvmt_off(curUnit,binActCC));

end

%% GLM data

% %{

GLMfil = VO.fracDepths >= fracDepthLims(1) & VO.fracDepths <= fracDepthLims(2) & ~avoidRec_VO_FM_perCell & VO.probeNum == 1 & VO.movieStatus < 3; % intitial SU filter

warning('off', 'all');

modelSpec_full = 'numSpikes ~ VtrialYN * FM_1'; % model spec for full GLM

% Get GLM data
pVal_V = nan(numUnits,1);
pVal_FM = nan(numUnits,1);
pVal_VxFM = nan(numUnits,1);
FRchange_GLM = nan(numUnits,1);
beta_V_full = nan(numUnits,1);
beta_FM_full = nan(numUnits,1);
beta_VxFM_full = nan(numUnits,1);
fprintf('\nGLM.  Looping through %s cells.  Completed: ',num2str(sum(numUnits)));
for i = 1:numUnits

    curNoOptoInd = ~VO.FRs_trialType_optoTrials{i};
    VFRs = VO.FRs_Vonset{i}(curNoOptoInd);
    numTrig = length(VFRs);
    ranFRs = VO.FRs_baseline{i}(curNoOptoInd);
    FRchange_GLM(i) = mean(VFRs) / mean(ranFRs);
    numSponTrig = length(ranFRs);
    numSpikes = [ VFRs*0.5 ; ranFRs*0.5 ]; % convert back to spike counts
    VtrialYN = [ true(numTrig,1) ; false(numSponTrig,1) ];
    FM_1 = [ VO.trialMvmtMean_VO{i}(curNoOptoInd,1) ; VO.trialMvmtMean_VO_baseline{i}(curNoOptoInd,1) ];
    FM_1 = (FM_1 - mean(FM_1)) ./ std(FM_1); % z-score the FM

    % GLM with V and FM - one PC
    tbl = table(VtrialYN, FM_1, numSpikes); % convert input data to a table for the GLM
    glmResults = fitglm(tbl,modelSpec_full,'Distribution','Poisson','Link','log','CategoricalVars','VtrialYN','ResponseVar','numSpikes'); % run GLM
    pVal_V(i) = glmResults.Coefficients.pValue(2);
    pVal_FM(i) = glmResults.Coefficients.pValue(3);
    pVal_VxFM(i) = glmResults.Coefficients.pValue(4);
    beta_V_full(i) = glmResults.Coefficients.Estimate(2);
    beta_FM_full(i) = glmResults.Coefficients.Estimate(3);
    beta_VxFM_full(i) = glmResults.Coefficients.Estimate(4);

    % Counter
    counterStr(i);

end

warning('on', 'all');

% FDR correct p values
pVal_V_FDR = pVal_V;
pVal_V_FDR(GLMfil) = mafdr(pVal_V_FDR(GLMfil),'BHFDR',true);
pVal_FM_FDR = pVal_FM;
pVal_FM_FDR(GLMfil) = mafdr(pVal_FM_FDR(GLMfil),'BHFDR',true);
pVal_VxFM_FDR = pVal_VxFM;
pVal_VxFM_FDR(GLMfil) = mafdr(pVal_VxFM_FDR(GLMfil),'BHFDR',true);

% Visually responsive cells according to the GLM
visRespCells_GLM = GLMfil & pVal_V_FDR < alphaLevel & FRchange_GLM > 1;

% Beta values within limits
lims = [-1 1];
beta_V_full_lim = beta_V_full;
beta_V_full_lim(beta_V_full_lim<lims(1)) = lims(1);
beta_V_full_lim(beta_V_full_lim>lims(2)) = lims(2);
beta_FM_full_lim = beta_FM_full;
beta_FM_full_lim(beta_FM_full_lim<lims(1)) = lims(1);
beta_FM_full_lim(beta_FM_full_lim>lims(2)) = lims(2);
beta_VxFM_full_lim = beta_VxFM_full;
beta_VxFM_full_lim(beta_VxFM_full_lim<lims(1)) = lims(1);
beta_VxFM_full_lim(beta_VxFM_full_lim>lims(2)) = lims(2);

%}

%% Regression to the mean stats

% %{

fil_reg = VO.fracDepths >= fracDepthLims(1) & VO.fracDepths <= fracDepthLims(2) & ~avoidRec_VO_FM_perCell;

numIter = 1000;

% Random half and other half of visual only trials
numCell = length(VO.FRs_Vonset);
indsRanHalf = cell(numCell,numIter);
indsOtherHalf = cell(numCell,numIter);
fprintf('\nGetting indexes to random half and other half trials for %s iterations: ',num2str(numIter))
for j = 1:numIter
    for i = 1:numCell
        noOptoNumInd = find(indexToNonOptoCells{i}); % all visual trials from current cell without opto
        halfNumTrials = floor(length(noOptoNumInd)/2);
        curInds = noOptoNumInd(randperm(length(noOptoNumInd),halfNumTrials));  
        indsRanHalf{i,j} = false(length(indexToNonOptoCells{i}),1);
        indsRanHalf{i,j}(curInds) = true;
        indsOtherHalf{i,j} = indexToNonOptoCells{i} & ~indsRanHalf{i};
    end
    counterStr(j)
end

% FR change values (from baseline to Vonset), for all V trials, for the random half and the other half of trials - to create the null distribution
FRch_base_to_V_ranHalf = nan(numUnits,numIter);
FRch_base_to_V_otherHalf = nan(numUnits,numIter);
FRch_base_to_V_ranHalf_notLog = nan(numUnits,numIter);
FRch_base_to_V_otherHalf_notLog = nan(numUnits,numIter);
FRdiff_base_to_V_ranHalf = nan(numUnits,numIter);
FRdiff_base_to_V_otherHalf = nan(numUnits,numIter);
fprintf('\nGetting FR changes for random half and other half trials for %s iterations: ',num2str(numIter))
for i = 1:numIter
    % FR change
    FRch_base_to_V_ranHalf(:,i) = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, indsRanHalf(:,i), [-inf inf]);
    FRch_base_to_V_otherHalf(:,i) = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, indsOtherHalf(:,i), [-inf inf]);
    % FR change - not log2
    FRch_base_to_V_ranHalf_notLog(:,i) = getFRchangeVals(VO.FRs_baseline, VO.FRs_Vonset, indsRanHalf(:,i));
    FRch_base_to_V_otherHalf_notLog(:,i) = getFRchangeVals(VO.FRs_baseline, VO.FRs_Vonset, indsOtherHalf(:,i));
    % FR diff
    FRdiff_base_to_V_ranHalf(:,i) = getFRdiff(VO.FRs_baseline, VO.FRs_Vonset, indsRanHalf(:,i));
    FRdiff_base_to_V_otherHalf(:,i) = getFRdiff(VO.FRs_baseline, VO.FRs_Vonset, indsOtherHalf(:,i));
    counterStr(i)
end

% Visually evoked p-values for random half of trials
pVal_visResp_ranHalf = nan(numUnits,numIter);
fprintf('\nGetting p values for random half trials for %s iterations: ',num2str(numIter))
for j = 1:numIter
    for i = 1:numUnits
        pVal_visResp_ranHalf(i,j) = signrank(VO.FRs_baseline{i}(indsRanHalf{i,j}), VO.FRs_Vonset{i}(indsRanHalf{i,j}));
    end
    counterStr(j)
end

% FDR correct p values - correcting the VC and AC separately, and not including cells out of the responsive span when correcting
pVal_visResp_ranHalf_FDR = nan(numUnits,numIter);
fprintf('\nFDR correcting p values for random half trials for %s iterations: ',num2str(numIter))
for i = 1:numIter
    pVal_visResp_ranHalf_FDR(fil_reg & VO.probeNum == 1,i) = mafdr(pVal_visResp_ranHalf(fil_reg & VO.probeNum == 1,i),'BHFDR',true);
    pVal_visResp_ranHalf_FDR(fil_reg & VO.probeNum == 2,i) = mafdr(pVal_visResp_ranHalf(fil_reg & VO.probeNum == 2,i),'BHFDR',true);
    counterStr(i)
end

% Visually responsive cells using a random half of V trials
visRespCell_AC_ranHalf = repmat(fil_reg & VO.probeNum == 1,1,numIter) & pVal_visResp_ranHalf_FDR < alphaLevel ;

% Difference values for increased FR cells
diffFR_ranOther_inc = nan(numIter,1);
diffFR_sub_ranOther_inc = nan(numIter,1);
diffFR_ranOther_inc_notLog = nan(numIter,1);
for i = 1:numIter
    incFR = FRch_base_to_V_ranHalf(:,i) > 0;
    % FR change
    group1 = FRch_base_to_V_ranHalf(visRespCell_AC_ranHalf(:,i) & incFR,i);
    group2 = FRch_base_to_V_otherHalf(visRespCell_AC_ranHalf(:,i) & incFR,i);
    isInf = isinf(group1) | isinf(group2);
    diffFR_ranOther_inc(i) = median(group2(~isInf) - group1(~isInf), 'omitnan');
    % FR diff
    group1 = FRdiff_base_to_V_ranHalf(visRespCell_AC_ranHalf(:,i) & incFR,i);
    group2 = FRdiff_base_to_V_otherHalf(visRespCell_AC_ranHalf(:,i) & incFR,i);
    diffFR_sub_ranOther_inc(i) = median(group2 - group1, 'omitnan');
        
    % FR change - not log
    group1 = FRch_base_to_V_ranHalf_notLog(visRespCell_AC_ranHalf(:,i) & incFR,i);
    group2 = FRch_base_to_V_otherHalf_notLog(visRespCell_AC_ranHalf(:,i) & incFR,i);
    isInf = isinf(group1) | isinf(group2);
    diffFR_ranOther_inc_notLog(i) = median(group2(~isInf) - group1(~isInf), 'omitnan');
    
end

%}

%% Longer baseline

% %{

fprintf('\nGetting FM movement threshold crossing stuff for longer baseline...')

dontTakeAfterTrial = 2; % don't use baseline periods within this duration after a trial to avoid after-effects from previous trials

trialsOverThresh_CMN_longer = cell(numRec_VO,1);
trialsOverThresh_off_longer = cell(numRec_VO,1);
perOver_CMN_longer = nan(numRec_VO,1);
perOver_off_longer = nan(numRec_VO,1);
for i = 1:numRec_VO

    % CMN
    optoTrials = VO.FRs_trialType_optoTrials{recFirstPositions_VO(i)}; % find opto trials
    % Make the first trial, which may not have enough FM data preceding the
    % analysis window, an opto trial, so it doesn't get
    % included in all the stats below
    optoTrials(1) = true;
    curRecMvmts = VO.trialMvmtAll{recFirstPositions_VO(i)}; % CMN movement data for current recording
    preTimes = VO.durUntilPrevOffTrig_V{recFirstPositions_VO(i)};
    movements_zScore_temp = nan(size(curRecMvmts,1),size(curRecMvmts,2));
    for j = 1:size(curRecMvmts,1)
        startPoint = -preTimes(j) + dontTakeAfterTrial;
        baseInd = tAx_VO > max([tAx_VO(1) startPoint]) & tAx_VO < 0;
        curRecMvmts_bs = curRecMvmts(j,:,1) - mean(curRecMvmts(j,baseInd,1));
        movements_zScore_temp(j,:) = curRecMvmts_bs / std(curRecMvmts_bs(baseInd));
    end
    overTh = mean(movements_zScore_temp(:,CMNactWin_VO,:),2) > zScoreThresh;
    trialsOverThresh_CMN_longer{i} = overTh;
    perOver_CMN_longer(i) = ( sum(overTh(~optoTrials)) / sum(~optoTrials) ) * 100;

    % Noise bursts
    curRecMvmts = VO.trialMvmtAll_offset{recFirstPositions_VO(i)}; % CMN movement data for current recording
    preTimes = VO.durUntilPrevOffTrig_A{recFirstPositions_VO(i)};
    movements_zScore_temp = nan(size(curRecMvmts,1),size(curRecMvmts,2));
    for j = 1:size(curRecMvmts,1)
        startPoint = -preTimes(j) + dontTakeAfterTrial;
        baseInd = tAx_VO > max([tAx_VO(1) startPoint]) & tAx_VO < 0;
        curRecMvmts_bs = curRecMvmts(j,:,1) - mean(curRecMvmts(j,baseInd,1));
        movements_zScore_temp(j,:) = curRecMvmts_bs / std(curRecMvmts_bs(baseInd));
    end
    overTh = mean(movements_zScore_temp(:,CMNactWin_VO,:),2) > zScoreThresh;
    trialsOverThresh_off_longer{i} = overTh;
    perOver_off_longer(i) = ( sum(overTh) / length(overTh) ) * 100;

end

% Spiking stuff for movement and no movement trials

trialsOverThresh_longerB_perCell = repelem(trialsOverThresh_CMN_longer,numCellPerRec_VO,1);
index_noMvmt_longerB = cellfun(@(x) sum(x(:,:,1:numPCsToUse),3) == 0  , trialsOverThresh_longerB_perCell , 'UniformOutput',  false );
index_noOpto_noMvmt_longerB = cellfun(@(x,y) ~x & sum(y(:,:,1:numPCsToUse),3) == 0  , ...
    VO.FRs_trialType_optoTrials, trialsOverThresh_longerB_perCell , ...
    'UniformOutput',  false );
index_noOpto_mvmt_longerB = cellfun(@(x,y) ~x & sum(y(:,:,1:numPCsToUse),3) > 0  , ...
    VO.FRs_trialType_optoTrials, trialsOverThresh_longerB_perCell , ...
    'UniformOutput',  false );
[pVals_lowVsHighMvmt_longerB, pVals_lowVsHighMvmt_longerB_perm] = getLowVsHighMvmtPvals(VO.FRs_baseline, VO.FRs_Vonset, VO.FRs_trialType_optoTrials, index_noMvmt_longerB);
FRch_base_to_V_lowMvmt_longerB = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_noMvmt_longerB, FRchangeLims);
FRch_base_to_V_highMvmt_longerB = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_mvmt_longerB, FRchangeLims);

fprintf('\n...complete')

%}

%% FM ROIs

FRchangeLims_ROI = [-inf inf];
% FRchangeLims_ROI = FRchangeLims;

fprintf('\nGetting FM movement threshold crossing stuff for ROIs...')

trialsOverThresh_CMN_ROI = cell(numRec_VO,4);
trialsOverThresh_off_ROI = cell(numRec_VO,4);
perOver_CMN_ROI = nan(numRec_VO,4);
perOver_off_ROI = nan(numRec_VO,4);
% abs
trialsOverThresh_CMN_ROI_abs = cell(numRec_VO,4);
trialsOverThresh_off_ROI_abs = cell(numRec_VO,4);
perOver_CMN_ROI_abs = nan(numRec_VO,4);
perOver_off_ROI_abs = nan(numRec_VO,4);
for i = 1:numRec_VO

    optoTrials = VO.FRs_trialType_optoTrials{recFirstPositions_VO(i)}; % find opto trials

    for j = 1:4 % Loop through the 4 FM ROIs

        % CMN
        curRecMvmts = VO.trialMvmtROI{recFirstPositions_VO(i)}(:,:,j); % CMN movement data for current recording and ROI
        curRecMvmts_bs = curRecMvmts - mean(curRecMvmts(:,baseInd_VO,:),2); % baseline subtract
        curRecMvmts_zScore = curRecMvmts_bs ./ std(curRecMvmts_bs(:,baseInd_VO,:),[],2);
        % Threshold
        overTh = mean(curRecMvmts_zScore (:,CMNactWin_VO,:),2) > zScoreThresh;
        trialsOverThresh_CMN_ROI{i,j} = overTh;
        perOver_CMN_ROI(i,j) = ( sum(overTh(~optoTrials)) / sum(~optoTrials) ) * 100;
        % abs
        curRecMvmts_zScore_abs = abs(curRecMvmts_zScore);
        overTh = mean(curRecMvmts_zScore_abs(:,CMNactWin_VO,:),2) > zScoreThresh;
        trialsOverThresh_CMN_ROI_abs{i,j} = overTh;
        perOver_CMN_ROI_abs(i,j) = ( sum(overTh(~optoTrials)) / sum(~optoTrials) ) * 100;

        % Noise bursts
        curRecMvmts = VO.trialMvmtROI_offset{recFirstPositions_VO(i)}(:,:,j); 
        curRecMvmts_bs = curRecMvmts - mean(curRecMvmts(:,baseInd_VO,:),2);
        curRecMvmts_zScore = curRecMvmts_bs ./ std(curRecMvmts_bs(:,baseInd_VO,:),[],2);
        % Threshold
        overTh = mean(curRecMvmts_zScore (:,CMNactWin_VO,:),2) > zScoreThresh;
        trialsOverThresh_off_ROI{i,j} = overTh;
        perOver_off_ROI(i,j) = ( sum(overTh) / length(overTh) ) * 100;
        % abs
        curRecMvmts_zScore_abs = abs(curRecMvmts_zScore);
        overTh = mean(curRecMvmts_zScore_abs(:,CMNactWin_VO,:),2) > zScoreThresh;
        trialsOverThresh_off_ROI_abs{i,j} = overTh;
        perOver_off_ROI_abs(i,j) = ( sum(overTh) / length(overTh) ) * 100;

    end

end

pVals_lowVsHighMvmt_ROI = nan(sum(numCellPerRec_VO),4);
pVals_lowVsHighMvmt_ROI_perm = nan(sum(numCellPerRec_VO),4);
FRch_base_to_V_lowMvmt_ROI = nan(sum(numCellPerRec_VO),4);
FRch_base_to_V_highMvmt_ROI = nan(sum(numCellPerRec_VO),4);
% abs
pVals_lowVsHighMvmt_ROI_abs = nan(sum(numCellPerRec_VO),4);
pVals_lowVsHighMvmt_ROI_abs_perm = nan(sum(numCellPerRec_VO),4);
FRch_base_to_V_lowMvmt_ROI_abs = nan(sum(numCellPerRec_VO),4);
FRch_base_to_V_highMvmt_ROI_abs = nan(sum(numCellPerRec_VO),4);
for i = 1:4 % loop through each FM ROI
    trialsOverThresh_ROI_perCell = repelem(trialsOverThresh_CMN_ROI(:,i),numCellPerRec_VO,1);
    index_noMvmt_ROI = cellfun(@(x) sum(x(:,:,1:numPCsToUse),3) == 0  , trialsOverThresh_ROI_perCell , 'UniformOutput',  false );
    index_noOpto_noMvmt_ROI = cellfun(@(x,y) ~x & sum(y(:,:,1:numPCsToUse),3) == 0  , ...
        VO.FRs_trialType_optoTrials, trialsOverThresh_ROI_perCell , ...
        'UniformOutput',  false );
    index_noOpto_mvmt_ROI = cellfun(@(x,y) ~x & sum(y(:,:,1:numPCsToUse),3) > 0  , ...
        VO.FRs_trialType_optoTrials, trialsOverThresh_ROI_perCell , ...
        'UniformOutput',  false );
    [pVals_lowVsHighMvmt_ROI(:,i), pVals_lowVsHighMvmt_ROI_perm(:,i)] = getLowVsHighMvmtPvals(VO.FRs_baseline, VO.FRs_Vonset, VO.FRs_trialType_optoTrials, index_noMvmt_ROI);
    FRch_base_to_V_lowMvmt_ROI(:,i) = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_noMvmt_ROI, FRchangeLims_ROI);
    FRch_base_to_V_highMvmt_ROI(:,i) = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_mvmt_ROI, FRchangeLims_ROI);
    % abs
    trialsOverThresh_ROI_perCell = repelem(trialsOverThresh_CMN_ROI_abs(:,i),numCellPerRec_VO,1);
    index_noMvmt_ROI_abs = cellfun(@(x) sum(x(:,:,1:numPCsToUse),3) == 0  , trialsOverThresh_ROI_perCell , 'UniformOutput',  false );
    index_noOpto_noMvmt_ROI_abs = cellfun(@(x,y) ~x & sum(y(:,:,1:numPCsToUse),3) == 0  , ...
        VO.FRs_trialType_optoTrials, trialsOverThresh_ROI_perCell , ...
        'UniformOutput',  false );
    index_noOpto_mvmt_ROI_abs = cellfun(@(x,y) ~x & sum(y(:,:,1:numPCsToUse),3) > 0  , ...
        VO.FRs_trialType_optoTrials, trialsOverThresh_ROI_perCell , ...
        'UniformOutput',  false );
    [pVals_lowVsHighMvmt_ROI_abs(:,i), pVals_lowVsHighMvmt_ROI_abs_perm(:,i)] = getLowVsHighMvmtPvals(VO.FRs_baseline, VO.FRs_Vonset, VO.FRs_trialType_optoTrials, index_noMvmt_ROI_abs);
    FRch_base_to_V_lowMvmt_ROI_abs(:,i) = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_noMvmt_ROI_abs, FRchangeLims_ROI);
    FRch_base_to_V_highMvmt_ROI_abs(:,i) = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_mvmt_ROI_abs, FRchangeLims_ROI);
end

fprintf('\n...complete')

%% Ball movement

fprintf('\nGetting ball movement threshold crossing stuff...')

FRchangeLims_ball = [-inf inf];
% FRchangeLims_ball = FRchangeLims;

trialsOverThresh_CMN_ball = cell(numRec_VO,1);
trialsOverThresh_off_ball = cell(numRec_VO,1);
perOver_CMN_ball = nan(numRec_VO,1);
perOver_off_ball = nan(numRec_VO,1);
% abs
trialsOverThresh_CMN_ball_abs = cell(numRec_VO,1);
trialsOverThresh_off_ball_abs = cell(numRec_VO,1);
perOver_CMN_ball_abs = nan(numRec_VO,1);
perOver_off_ball_abs = nan(numRec_VO,1);
for i = 1:numRec_VO

    optoTrials = VO.FRs_trialType_optoTrials{recFirstPositions_VO(i)}; % find opto trials

    % CMN
    curRecMvmts = VO.trialMvmtBall{recFirstPositions_VO(i)}; % CMN movement data for current recording and ROI
    curRecMvmts_bs = curRecMvmts - mean(curRecMvmts(:,baseInd_VO,:),2); % baseline subtract
    curRecMvmts_zScore = curRecMvmts_bs ./ std(curRecMvmts_bs(:,baseInd_VO,:),[],2);
    % Threshold
    overTh = mean(curRecMvmts_zScore (:,CMNactWin_VO,:),2) > zScoreThresh;
    trialsOverThresh_CMN_ball{i} = overTh;
    perOver_CMN_ball(i) = ( sum(overTh(~optoTrials)) / sum(~optoTrials) ) * 100;
    % abs
    curRecMvmts_zScore_abs = abs(curRecMvmts_zScore);
    overTh = mean(curRecMvmts_zScore_abs(:,CMNactWin_VO,:),2) > zScoreThresh;
    trialsOverThresh_CMN_ball_abs{i} = overTh;
    perOver_CMN_ball_abs(i) = ( sum(overTh(~optoTrials)) / sum(~optoTrials) ) * 100;

    % Offset
    curRecMvmts = VO.trialMvmtBallOff{recFirstPositions_VO(i)}; % CMN movement data for current recording and ROI
    curRecMvmts_bs = curRecMvmts - mean(curRecMvmts(:,baseInd_VO,:),2); % baseline subtract
    curRecMvmts_zScore = curRecMvmts_bs ./ std(curRecMvmts_bs(:,baseInd_VO,:),[],2);
    % Threshold
    overTh = mean(curRecMvmts_zScore (:,CMNactWin_VO,:),2) > zScoreThresh;
    trialsOverThresh_off_ball{i} = overTh;
    perOver_off_ball(i) = ( sum(overTh) / length(overTh) ) * 100;
    % abs
    curRecMvmts_zScore_abs = abs(curRecMvmts_zScore);
    overTh = mean(curRecMvmts_zScore_abs(:,CMNactWin_VO,:),2) > zScoreThresh;
    trialsOverThresh_off_ball_abs{i} = overTh;
    perOver_off_ball_abs(i) = ( sum(overTh) / length(overTh) ) * 100;

end

trialsOverThresh_ball_perCell = repelem(trialsOverThresh_CMN_ball,numCellPerRec_VO,1);
index_noMvmt_ball = cellfun(@(x) sum(x(:,:,1:numPCsToUse),3) == 0  , trialsOverThresh_ball_perCell , 'UniformOutput',  false );
index_noOpto_noMvmt_ball = cellfun(@(x,y) ~x & sum(y(:,:,1:numPCsToUse),3) == 0  , ...
    VO.FRs_trialType_optoTrials, trialsOverThresh_ball_perCell , ...
    'UniformOutput',  false );
index_noOpto_mvmt_ball = cellfun(@(x,y) ~x & sum(y(:,:,1:numPCsToUse),3) > 0  , ...
    VO.FRs_trialType_optoTrials, trialsOverThresh_ball_perCell , ...
    'UniformOutput',  false );
[pVals_lowVsHighMvmt_ball, pVals_lowVsHighMvmt_ball_perm] = getLowVsHighMvmtPvals(VO.FRs_baseline,VO.FRs_Vonset,VO.FRs_trialType_optoTrials,index_noMvmt_ball);
FRch_base_to_V_lowMvmt_ball = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_noMvmt_ball, FRchangeLims_ball);
FRch_base_to_V_highMvmt_ball = getFRchangeVals_log2(VO.FRs_baseline, VO.FRs_Vonset, index_noOpto_mvmt_ball, FRchangeLims_ball);
% abs
trialsOverThresh_ball_perCell = repelem(trialsOverThresh_CMN_ball_abs,numCellPerRec_VO,1);
index_noMvmt_ball_abs = cellfun(@(x) sum(x(:,:,1:numPCsToUse),3) == 0  , trialsOverThresh_ball_perCell , 'UniformOutput',  false );
index_noOpto_noMvmt_ball_abs = cellfun(@(x,y) ~x & sum(y(:,:,1:numPCsToUse),3) == 0  , ...
    VO.FRs_trialType_optoTrials, trialsOverThresh_ball_perCell , ...
    'UniformOutput',  false );
index_noOpto_mvmt_ball_abs = cellfun(@(x,y) ~x & sum(y(:,:,1:numPCsToUse),3) > 0  , ...
    VO.FRs_trialType_optoTrials, trialsOverThresh_ball_perCell , ...
    'UniformOutput',  false );
[pVals_lowVsHighMvmt_ball_abs, pVals_lowVsHighMvmt_ball_abs_perm] = getLowVsHighMvmtPvals(VO.FRs_baseline,VO.FRs_Vonset,VO.FRs_trialType_optoTrials,index_noMvmt_ball_abs);
FRch_base_to_V_lowMvmt_ball_abs = getFRchangeVals_log2(VO.FRs_baseline,VO.FRs_Vonset,index_noOpto_noMvmt_ball_abs,FRchangeLims_ball);
FRch_base_to_V_highMvmt_ball_abs = getFRchangeVals_log2(VO.FRs_baseline,VO.FRs_Vonset,index_noOpto_mvmt_ball_abs,FRchangeLims_ball);

fprintf('\n...complete')

%% Save worksapce

% workspaceSaveDateTime = sprintf('allVars_%s.mat',datetime('now','Format','yyyy-MM-dd_HH-mm-ss'));
% saveFol = '\\heward.cin.ucsf.edu\work\Tim\Projects\veryNew\visualNeuronsAC\visualAndMovementPaper\workSpaceSave';
% save(fullfile(saveFol,workspaceSaveDateTime),'-v7.3');
% 
% tEnd = toc;
% fprintf('\nThat took %s minutes to run',num2str(round(tEnd/60,1)));


