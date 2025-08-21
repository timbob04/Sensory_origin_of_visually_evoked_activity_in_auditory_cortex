clear all; close all; clc; restoredefaultpath

cd('\\heward.cin.ucsf.edu\Work\Tim\Projects\veryNew\visualNeuronsAC\code\VO_CMNonly\collectStats')

%% Inputs

DOCID = '1BcTbmKkqf761zizHoId0gt9dLbF5Kc6HlF2vKy7BzWM'; % google sheet with details for all VO recordings
DOCID_bigLB = '1lJlNOdE3_iISDkPPPDakdCS528kprU_qTSoXKAm-FhM'; % google sheet with lab book
dataFol = '\\heward.cin.ucsf.edu\Data\Tim\NP'; % where the NP data live
saveFol = '\\heward.cin.ucsf.edu\Work\Tim\Projects\veryNew\visualNeuronsAC\code\VO_CMNonly\summaryData';
localDataFol = 'E:\Data'; % where the local data lives, to get TDT data

labBook = gsheet2table(DOCID,'0'); % read the google sheet
labBook_big = gsheet2table(DOCID_bigLB,'0'); % read the google sheet

% Some timing info for PSTHs and other things
TPWthresh = 0.00045; % for splitting BS and NS units (seconds)
CMNdur = .5; % duration of CMN
optoPreTime = .5; % duration the opto precedes the CMN
timeWin = [-6 2]; % time window to collect PSTH and raster data
PSTHbinSize = 0.05; % bin size for PSTH generation (in sec)
timeWin_Vonset = [ 0 CMNdur ]; % time window for the visual response
timeWin_baseline = [ -CMNdur 0 ]; % [ -1 -optoPreTime ] time window for baseline.  Use [-1 -optoPreTime] if I want to see how the data looks when avoiding the opto-period
lenOffStim = 1.05; % duration of the offet stimulus (trigger on to trigger off)
timeWin_baseline_off_s = [-.15 0];
timeWin_baseline_off_l = [-1.05 0];
timeWin_off_s = [0 .15];
timeWin_off_l = [0 1.05];

% SU waveform stats inputs
numEachWay = [1 2 5]; % get the slope at 0.5 ms using this many samples each way

% Some timing info for movement-based stuff
frameSP = 0.0332; % frame sampling period (1/fps)

% Here control which recordings to analyze, using info in the VO lab book (DOCID)
recsToAnalyze = strcmp(labBook.Ephys,'1') ... % has ephys recording
    & strcmp(labBook.ExptType,'VO_CMN_flash_NP') ... % VO_CMN_flash experiment type
    & strcmp(labBook.AC_sorted,'1') ... % AC spikes sorted
    & strcmp(labBook.VC_sorted,'1') ... % VC spikes sorted
    & strcmp(labBook.Face_movie,'1') ... % contains face movie
    & strcmp(labBook.AnalyzeV,'1'); % analyze Y/N. Sometimes 'no' if something went wrong

% Number of FM PCs to look at
numOfPCsToCheck = 10;

% Inputs for designating a recording site as primary- or secondary-like
baseWin_T = [-.2 0];
actWin_T = [0 .1];
lenStim = .1;
latTh = 2.5;
attenTake = 1:7;
binTimeLat = 0.002;
binTimeLatUS = 0.0001;

% ISI violation time threshold
ISIviolTime = 0.001;

%% Some starter stats

numFramesToPlot = floor(diff(timeWin) / frameSP); % total number of FM frames to store for each trial
binCenters_faceMvmts = timeWin(1)+frameSP/2 : frameSP : timeWin(1) + diff(timeWin) - (frameSP/2); % bin centers for the face movement data

numFramesForMean = floor(CMNdur / frameSP);

%% Get the stats for each recording

recsToAnalyze_ind = find(recsToAnalyze); % index to the recordings to analyze, in the VO lab book

numRec = sum(recsToAnalyze); % number of recordings to process

% Initialize variables for stats collection
numUnit = nan(numRec*2,1);
V1_silencing_manDec = nan(numRec*2,1); % manual decision on whether there was silencing in the VC
analyzeVC = nan(numRec*2,1); % analyze opto data for summary plots (Y/N) - maybe something went wrong with the opto, or some timing was change, but I still want to analyze the V responses
signsOfBadExpression = nan(numRec*2,1);
recNumber = nan(numRec*2,1);
mouseID = cell(numRec*2,1);
probeNum = nan(numRec*2,1);
Ai32Mouse = false(numRec*2,1);
wildTypeMouse = false(numRec*2,1);
num_Vtrig = nan(numRec*2,1);
num_optoTrials = nan(numRec*2,1);
% Metadata
unitIDs = cell(numRec*2,1);
BSunits = cell(numRec*2,1);
NSunits = cell(numRec*2,1);
medUnits = cell(numRec*2,1);
fracDepths = cell(numRec*2,1);
channelsCol = cell(numRec*2,1);
fracDepthDataCol = cell(numRec*2,1);
fracDepth_recID = nan(numRec*2,1);
fracDepth_probeCol = nan(numRec*2,1);
% PSTHs and rasters
PSTHout_V = cell(numRec*2,1);
raster_times_V = cell(numRec*2,1);
raster_trInd_V = cell(numRec*2,1);
PSTHout_VO = cell(numRec*2,1);
raster_times_VO = cell(numRec*2,1);
raster_trInd_VO = cell(numRec*2,1);
PSTHout_off = cell(numRec*2,1);
raster_times_off = cell(numRec*2,1);
raster_trInd_off = cell(numRec*2,1);
% FRs
FRs_baseline = cell(numRec*2,1);
FRs_Vonset = cell(numRec*2,1);
FRs_trialType_optoTrials = cell(numRec*2,1);
% Movement stuff
trialMvmtAll = cell(numRec*2,1);
trialMvmtROI = cell(numRec*2,1);
movieStatus = nan(numRec*2,1);
ROIdataExist = nan(numRec*2,1);
trialMvmtMean_VO = cell(numRec*2,1);
trialMvmtMean_VO_baseline = cell(numRec*2,1);
% Offset stuff
trialMvmtAll_offset = cell(numRec*2,1);
trialMvmtROI_offset = cell(numRec*2,1);
FRs_off_baseline_s = cell(numRec*2,1);
FRs_off_baseline_l = cell(numRec*2,1);
FRs_off_baseline_CMNlen = cell(numRec*2,1);
FRs_off_act_s = cell(numRec*2,1);
FRs_off_act_l = cell(numRec*2,1);
FRs_off_act_CMNlen = cell(numRec*2,1);
offSpikingResponses = cell(numRec*2,1);
% Eigen values
eigenValues_varExp = cell(numRec,1);
eigenValues_recID = nan(numRec,1);
% SU waveeform data
TPWs = cell(numRec,1);
troughHeightRatios = cell(numRec,1);
endSlopes = cell(numRec,1);
waveforms_mean = cell(numRec,1);
waveforms_SD = cell(numRec,1);
ISIviolPer = cell(numRec,1);
% MU tone stuff for primary/secondary designation
latencies_T = cell(numRec*2,1);
pVals_T = cell(numRec*2,1);
pVals_twoSided_T = cell(numRec*2,1);
overHalfHeight_T = cell(numRec*2,1);
varFrac_T = cell(numRec*2,1);
depthWithin_T = cell(numRec*2,1);
% Opto voltage
maxOptoVol = nan(numRec*2,1);
% Store frame times and FM PC1 (whole trace)
FM_PC1_store = cell(numRec,1);
FM_ROI_store = cell(numRec,1);
frameTimes_store = cell(numRec,1);
% Store V and offset trigger on times
trigOn_store_VO = cell(numRec,1);
trigOn_store_off = cell(numRec,1);
% Duration until previous trigger off
durUntilPrevOffTrig_col_V = cell(numRec*2,1);
durUntilPrevOffTrig_col_A = cell(numRec*2,1);
% Ball movement
trialMvmtBall = cell(numRec*2,1);
trialMvmtBallOff = cell(numRec*2,1);
minTimeDiff_V = cell(numRec*2,1);
minTimeDiff_A = cell(numRec*2,1);
treadmill_vel_store = cell(numRec,1);
treadmill_time_store = cell(numRec,1);

fprintf('\nCollecting AVO stats for %s recordings',num2str(numRec*2));
tic
k = 1;
recCounter = 1;
tuningPass = false(1,numRec*2);
for i = 1:numRec % loop through all recordings to analyze

    curInd = recsToAnalyze_ind(i); % index to row in lab book for current recording

    for curProbe = 1:2 % loop through probes 1 (VC) and 2 (AC)

        %% Some basic recording info (date, rec start time, etc)

        curDate = labBook.ExptDate{curInd}; % date of rec
        rec = labBook.Multiblock{curInd}; % recording time stamp
        blk = labBook.Block{curInd}; % block time
        V1_silencing_manDec(k) = strcmp(labBook.V1_silencing_manDec{curInd},'1'); % Is there V1 silencing?  Manual decision written in the lab book
        analyzeVC(k) = strcmp(labBook.AnalyzeVC{curInd},'1'); % can be set to zero for multiple reasons, including there being no recognizable cortex boundaries
        signsOfBadExpression(k) = strcmp(labBook.SignsOfBadExpression{curInd},'1');
        probeNum(k) = curProbe;
        recNumber(k) = str2double(labBook.RecID{curInd});
        mouseID{k} = labBook.MouseID{curInd};
        Ai32Mouse(k) = strcmp(labBook.Strain{curInd},'Ai32/PV-Cre');
        wildTypeMouse(k) = strcmp(labBook.Strain{curInd},'C57 WT');

        %% Load channel map

        chanMapDirAndName = fullfile(dataFol,curDate,sprintf('%s_g0_t0.imec%s.ap_kilosortChanMap.mat',rec,num2str(curProbe-1)));
        curChanMap = load(chanMapDirAndName);
        numChan = length(curChanMap.connected)-1;

        %% Spiking data

        spkData = load( fullfile(dataFol,curDate, sprintf('SUs_checked_%s_%s_to_%s_probe%s_chanFixed.mat',rec,blk,blk,num2str(curProbe)) ) ); % load the spiking data
        spkTimes = spkData.spikeTimes; % unpack the spike times
        channels = spkData.channel; % unpack the SU channels
        unitID = spkData.unitID;
        numUnits = length(spkTimes);
        medUnit = spkData.mediumUnits; % the units which were marked 'medium' in Phy or when doing my final SU plot check
        BS = spkData.troughtToPeakWidthTime > TPWthresh;
        NS = spkData.troughtToPeakWidthTime <= TPWthresh;
        % Store
        medUnits{k} = ismember(unitID,medUnit);
        numUnit(k) = numUnits;
        unitIDs{k} = unitID;
        BSunits{k} = BS';
        NSunits{k} = NS';
        channelsCol{k} = channels;
        % Get spike waveform properties
        TPWs{k} = spkData.troughtToPeakWidthTime'; % trough to peak widths
        troughHeightRatios{k} = transpose( abs(spkData.peakSampleAndAmp(2,:)) ./ abs(spkData.troughSampleAndAmp(2,:) ) ) ; % the peak amplitude divided by the trough amplitude
        slopes = nan(numUnits,length(numEachWay)); % for collecting the slope difference in amplitude values
        for curUnit = 1:numUnits
            if curUnit == 1
                Fs = spkData.Fs; SP = 1/Fs;
                numSampsPointFive = round(0.0005/SP);
                lenWaveform = length(spkData.waveforms{curUnit}.meanWaveform);
                curRecWaveforms_mean = nan(numUnits,lenWaveform);
                curRecWaveforms_SD = nan(numUnits,lenWaveform);
                waveFormTax = SP:SP:SP*lenWaveform;
            end
            curWaveform = spkData.waveforms{curUnit}.meanWaveform;
            curWaveformSD = spkData.waveforms{curUnit}.SDwaveform;
            peakPoint = spkData.troughSampleAndAmp(1,curUnit);
            for ew = 1:length(numEachWay)
                ampDiff = curWaveform(peakPoint+numSampsPointFive+numEachWay(ew)) - curWaveform(peakPoint+numSampsPointFive-numEachWay(ew));
                timeDiff = ((numEachWay(ew)*2)*SP) ;
                slopes(curUnit,ew) = ampDiff / timeDiff ;
            end
            curRecWaveforms_mean(curUnit,:) = curWaveform;
            curRecWaveforms_SD(curUnit,:) =  curWaveformSD;
        end
        endSlopes{k} = slopes; % slope difference in amplitude at 0.5 ms
        waveforms_mean{k} = curRecWaveforms_mean;
        waveforms_SD{k} = curRecWaveforms_SD;
        ISIviolPer_cur = nan(numUnits,1);
        for curUnit = 1:numUnits
            diffs = diff(spkData.spikeTimes{curUnit});
            ISIviolPer_cur(curUnit) = (sum(diffs<ISIviolTime) / length(diffs))*100;
        end
        ISIviolPer{k} = ISIviolPer_cur;

        %% Trigger data

        if curProbe == 1 % only load trigger stuff if on probe 1. Don't need to load again for probe 2.
            trda = load(fullfile(dataFol,curDate,sprintf('trigData_%s.mat',rec))); % load the trigger times
            % Block on and off times
            [secOn,secOff] = getBlockTime_stimMac(trda.blockOn_time,trda.blockOff_time,rec,blk,blk,false);
            indTrigWithinBlockTimes = trda.trigOn_time > secOn & trda.trigOn_time < secOff;
            % VO trigger on and off times
            trigOn_VO = trda.trigOn_time(indTrigWithinBlockTimes & trda.trigOn_numTrig == 2 ); % VO trigger on times
            trigOff_VO = trigOn_VO + CMNdur;
            numVtrig = length(trigOn_VO);
            % Opto trials
            optoTrials = sum(trda.optoOn_time' > trigOn_VO - optoPreTime - .2 & trda.optoOn_time' < trigOn_VO) == 1; % the V trials with an opto stimulus
            % Offset on and off times
            trigOn_offset = trda.trigOn_time(indTrigWithinBlockTimes & trda.trigOn_numTrig == 4 );
            trigOff_offset = trigOn_offset + lenOffStim;
            % Store trigger times
            trigOn_store_VO{i} = trigOn_VO;
            trigOn_store_off{i} = trigOn_offset;
            % All trigger offs
            allTrigOn = [ trigOn_VO trigOn_offset];
            allTrigOff = [ trigOff_VO trigOff_offset];
        end
        num_Vtrig(k) = numVtrig;
        num_optoTrials(k) = sum(optoTrials);

        %% Opto voltage level

        if curProbe == 1
            TDTname = labBook_big.Multiblock_TDT{find(strcmp(labBook_big.ExptDate,curDate) & strcmp(labBook_big.Multiblock_NP,rec),1,'first')};
            TDTfileName = fullfile(dataFol,curDate,TDTname);
            Lramp = TDTbin2mat(TDTfileName,'T1',secOn+100,'T2',secOn+300, 'TYPE', {'streams'}, 'STORE', 'Lrap');
            ligg = TDTbin2mat(TDTfileName,'T1',secOn+100,'T2',secOn+300, 'TYPE', {'streams'}, 'STORE', 'ligg');
            l_trigOn = detectRisingEdge(ligg.streams.ligg.data,.5);
            l_trigOff = detectRisingEdge(-ligg.streams.ligg.data,.5);
        end
        maxOptoVol(k) = max(Lramp.streams.Lrap.data(l_trigOn(end):l_trigOff(end)));

        %% Get the fractional depth for each unit

        fracDepths{k} = nan(numUnits,1);
        fracDepthData = nan(2,numChan);
        fracDepthFile = fullfile(dataFol,curDate,sprintf('fracDepth_%s_probe%s.mat',rec,num2str(curProbe)));
        if exist(fracDepthFile,'file') > 0
            load(fracDepthFile)
            fracDepthData = getFracDepth_NP1(sh,dp,numChan);
            [~,indMx] = max(channels' == fracDepthData(2,:)');
            fracDepths{k} = fracDepthData(1,indMx)';
        end
        fracDepthDataCol{k} = fracDepthData;
        fracDepth_recID(k) = recNumber(k);
        fracDepth_probeCol(k) = curProbe;

        %% Load face movement data

        if curProbe == 1 % only load trigger stuff if on probe 1. Don't need to load again for probe 2.
            ROIsvd = [];
            movieStatus_temp = [];
            ROIdataExist_temp = [];
            facemapFile = fullfile(dataFol,curDate,sprintf('eyecam_%s_%s_proc.mat',rec,blk));
            facemapFileROI = fullfile(dataFol,curDate,sprintf('eyecam_%s_%s_proc_ROI.mat',rec,blk));
            frameTimeFile = fullfile(dataFol,curDate,sprintf('frameTimes_%s_%s.mat',rec,blk));
            if exist(facemapFile,'file') ~= 0 && exist(frameTimeFile,'file') ~= 0
                fprintf('\nLoading data for face movement analysis...')
                load(facemapFile,'motSVD_0'); % load the video motion data
                % If the frame triggers exist, load them. Otherwise make them
                if exist(frameTimeFile,'file') ~= 0
                    load(frameTimeFile)
                    frameTimes = frameTimesOut;
                end
                % Store PC1 and frame times
                FM_PC1_store{i} = motSVD_0(:,1);
                frameTimes_store{i} = frameTimes;
                % Check that the number of video frames and the number of frame
                % triggers is the same
                if size(motSVD_0,1) ~= length(frameTimes) 
                    if abs(diff([size(motSVD_0,1) length(frameTimes)])) == 1 % if they differ by one frame
                        minLen = min([size(motSVD_0,1) length(frameTimes)]);
                        motSVD_0 = motSVD_0(1:minLen,:);
                        frameTimes = frameTimes(1:minLen);
                        movieStatus_temp = 2;
                    else % if they differ by more than one frame
                        warning('The number to frame trigger times and movie SVD points differs by more than 1')
                        movieStatus_temp = 3;
                    end
                else % if they are the same
                    movieStatus_temp = 1;
                end
                fprintf('complete')
                % Get the eigen values
                load(facemapFile,'motSv'); % load the video motion data
                eigenValues = motSv.^2;
                eigenValues_fraction = eigenValues / sum(eigenValues);
                eigenValues_varExp{recCounter} = eigenValues_fraction * 100;
                eigenValues_recID(recCounter) = recNumber(k);
            else
                movieStatus_temp = 4; % missing facemap or frame trigger file
            end
            % ROIs
            if exist(facemapFileROI,'file') ~= 0 && exist(frameTimeFile,'file') ~= 0
                ROIdataExist_temp = true;
                load(facemapFileROI,'motSVD_1','motSVD_2','motSVD_3','motSVD_4'); % load the video motion data
                if size(motSVD_0,1) ~= length(frameTimes)
                    if abs(diff([size(motSVD_0,1) length(frameTimes)])) == 1
                        motSVD_1 = motSVD_1(1:minLen,:);
                        motSVD_2 = motSVD_2(1:minLen,:);
                        motSVD_3 = motSVD_3(1:minLen,:);
                        motSVD_4 = motSVD_4(1:minLen,:);
                    else
                        motSVD_1 = nan(length(frameTimes),4);
                        motSVD_2 = nan(length(frameTimes),4);
                        motSVD_3 = nan(length(frameTimes),4);
                        motSVD_4 = nan(length(frameTimes),4);
                    end
                end
                ROIsvd = [ motSVD_1(:,1) motSVD_2(:,1) motSVD_3(:,1) motSVD_4(:,1) ];
                FM_ROI_store{i} = ROIsvd;
            else
                ROIdataExist_temp = false;
            end
        end
        movieStatus(k) = movieStatus_temp;
        ROIdataExist(k) = ROIdataExist_temp;

        %% Get the CMN trial-aligned face movements from PC1 - for plotting

        % Get the movements for the entire trial duration
        if curProbe == 1 % only load trigger stuff if on probe 1. Don't need to load again for probe 2.
            trialMvmtAll_temp = nan(numVtrig,numFramesToPlot,numOfPCsToCheck); % for storing the single trial face movements
            trialMvmtROI_temp = nan(numVtrig,numFramesToPlot,size(ROIsvd,2));
            durUntilPrevOffTrig_V = nan(numVtrig,1);
            if movieStatus(k) < 3
                for jj = 1:numOfPCsToCheck
                    vidMot = motSVD_0(:,jj)'; % get PC1 of the video motion
                    for j = 1:numVtrig
                        indStart = find(frameTimes > trigOn_VO(j) + timeWin(1),1,'first');
                        indEnd = indStart + numFramesToPlot - 1;
                        if indStart ~= 1 && indEnd <= length(vidMot)
                            trialMvmtAll_temp(j,:,jj) = vidMot(indStart:indEnd); % get the AV trial movements
                            if jj == 1
                                trigsBefore = allTrigOff < trigOn_VO(j);
                                if sum(trigsBefore) == 0
                                    durUntilPrevOffTrig_V(j) = -inf;
                                else
                                    durUntilPrevOffTrig_V(j) = trigOn_VO(j) - max(allTrigOff(trigsBefore));
                                end
                            end
                        elseif jj == 1
                            fprintf('\nTrial: %s.  CMN - starting or ending trial outside of frame trig bounds',num2str(j))
                        end
                    end
                end
                if ROIdataExist(k)
                    for jj = 1:size(ROIsvd,2)
                        vidMot = ROIsvd(:,jj)'; % get PC1 of the current ROI
                        for j = 1:numVtrig
                            indStart = find(frameTimes > trigOn_VO(j) + timeWin(1),1,'first');
                            indEnd = indStart + numFramesToPlot - 1;
                            if indStart ~= 1 && indEnd <= length(vidMot)
                                trialMvmtROI_temp(j,:,jj) = vidMot(indStart:indEnd);
                            elseif jj == 1
                                fprintf('\nTrial: %s.  CMN ROI - Starting or ending trial outside of frame trig bounds',num2str(j))
                            end
                        end
                    end
                end
            end
        end
        trialMvmtAll{k} = trialMvmtAll_temp;
        trialMvmtROI{k} = trialMvmtROI_temp;
        durUntilPrevOffTrig_col_V{k} = durUntilPrevOffTrig_V;

        %% Get the single-trial FRs for each window (also get the info for identifying what each trial type is)

        % Intialize the vectors for collecting the FRs
        FRs_baseline_temp = cell(numUnits,1); % For collecting the baseline FRs for the current recording's units
        FRs_Vonset_temp = cell(numUnits,1); % For collecting the visual onset FRs for the current recording's units
        % Define the trials to collect FRs from
        curTrigsOn = trigOn_VO;
        for j = 1:numUnits
            % Single-trial FRs for the baseline window
            timeWin_temp = timeWin_baseline;
            FRs_baseline_temp{j} = getSingleTrialFRs_NP(spkTimes{j},curTrigsOn,timeWin_temp);
            % Single-trial FRs for the visual onset window
            timeWin_temp = timeWin_Vonset;
            FRs_Vonset_temp{j} = getSingleTrialFRs_NP(spkTimes{j},curTrigsOn,timeWin_temp);
        end
        FRs_baseline{k} = FRs_baseline_temp;
        FRs_Vonset{k} = FRs_Vonset_temp;
        FRs_trialType_optoTrials{k} = repmat({optoTrials'},numUnits,1);

        %% Get the PSTHs and rasters for different trial types (V, VO, etc)

        % Initialize cell arrays for collecting PSTHs and rasters
        PSTHout_V_temp = cell(numUnits,1);
        PSTHout_VO_temp = cell(numUnits,1);
        raster_times_V_temp = cell(numUnits,1);
        raster_trInd_V_temp = cell(numUnits,1);
        raster_times_VO_temp = cell(numUnits,1);
        raster_trInd_VO_temp = cell(numUnits,1);
        for j = 1:numUnits
            % V triggers
            trigOnOffs = [ trigOn_VO(~optoTrials)' trigOff_VO(~optoTrials)'];
            % V PSTH
            [PSTHout_V_temp{j},binCenters_VO,~] = PSTHbasic_fixedBinSize_hColOut(spkTimes{j},trigOnOffs,-timeWin(1),PSTHbinSize,timeWin(2)-CMNdur,CMNdur);
            % V raster
            [rasterSpikesTime,raster_trInd_V_temp{j}] = getRaster(trigOnOffs(:,1)+timeWin(1),trigOnOffs(:,1)+timeWin(2),spkTimes{j});
            raster_times_V_temp{j} = rasterSpikesTime + timeWin(1);
            % VO triggers
            trigOnOffs = [ trigOn_VO(optoTrials)' trigOff_VO(optoTrials)'];
            % VO PSTH
            [PSTHout_VO_temp{j},~,~] = PSTHbasic_fixedBinSize_hColOut(spkTimes{j},trigOnOffs,-timeWin(1),PSTHbinSize,timeWin(2)-CMNdur,CMNdur);
            % VO raster
            [rasterSpikesTime,raster_trInd_VO_temp{j}] = getRaster(trigOnOffs(:,1)+timeWin(1),trigOnOffs(:,1)+timeWin(2),spkTimes{j});
            raster_times_VO_temp{j} = rasterSpikesTime + timeWin(1);
        end
        PSTHout_V{k} = cat(1,PSTHout_V_temp{:});
        raster_times_V{k} = raster_times_V_temp;
        raster_trInd_V{k} = raster_trInd_V_temp;
        PSTHout_VO{k} = cat(1,PSTHout_VO_temp{:});
        raster_times_VO{k} = raster_times_VO_temp;
        raster_trInd_VO{k} = raster_trInd_VO_temp;

        %% Mean face movement for 500 ms windows

        % Get the mean face movement for all visual triggers (0 to 0.5 sec), and the baselines
        curTrigOn = trigOn_VO;
        if curProbe == 1 % only load trigger stuff if on probe 1. Don't need to load again for probe 2.
            trialMvmtAll_temp_VO = nan(length(curTrigOn),numOfPCsToCheck); % for storing the single trial face movement
            trialMvmtAll_temp_VO_baseline = nan(length(curTrigOn),numOfPCsToCheck); % for storing the single trial baseline face movement
            if movieStatus(k) < 3
                for jj = 1:numOfPCsToCheck
                    vidMot = motSVD_0(:,jj)'; % get PC1 of the video motion
                    for j = 1:length(curTrigOn)
                        % V window
                        indStart = find(frameTimes > curTrigOn(j),1,'first');
                        indEnd = indStart + numFramesForMean - 1;
                        trialMvmtAll_temp_VO(j,jj) = mean(vidMot(indStart:indEnd)); % get the AV trial movements
                        % Baseline
                        indStart = find(frameTimes > curTrigOn(j)+timeWin_baseline(1),1,'first');
                        indEnd = indStart + numFramesForMean - 1;
                        trialMvmtAll_temp_VO_baseline(j,jj) = mean(vidMot(indStart:indEnd)); % get the AV trial movements
                    end
                end
            end
        end
        trialMvmtMean_VO{k} = trialMvmtAll_temp_VO;
        trialMvmtMean_VO_baseline{k} = trialMvmtAll_temp_VO_baseline;

        %% Get the offset trial-aligned face movements from PC1 - for plotting

        trigOnOffs = [trigOn_offset' trigOff_offset']; % offset trigger on and off times
        % Get the movements for the entire trial duration
        numOffTrig = size(trigOnOffs,1);
        if curProbe == 1 % only load trigger stuff if on probe 1. Don't need to load again for probe 2.
            trialMvmtAll_temp_offset = nan(numOffTrig,numFramesToPlot,numOfPCsToCheck); % for storing the single trial face movements
            trialMvmtROI_temp_offset = nan(numOffTrig,numFramesToPlot,size(ROIsvd,2));
            durUntilPrevOffTrig_A = nan(numOffTrig,1);
            if movieStatus(k) < 3
                for jj = 1:numOfPCsToCheck
                    vidMot = motSVD_0(:,jj)'; % get PC1 of the video motion
                    for j = 1:numOffTrig
                        indStart = find(frameTimes > trigOnOffs(j,1) + timeWin(1),1,'first');
                        indEnd = indStart + numFramesToPlot - 1;
                        if indStart ~= 1 && indEnd <= length(vidMot)
                            trialMvmtAll_temp_offset(j,:,jj) = vidMot(indStart:indEnd); % get the AV trial movements
                            if jj == 1
                                trigsBefore = allTrigOff < trigOnOffs(j,1);
                                if sum(trigsBefore) == 0
                                    durUntilPrevOffTrig_A(j) = -inf;
                                else
                                    durUntilPrevOffTrig_A(j) = trigOnOffs(j,1) - max(allTrigOff(trigsBefore));
                                end
                            end
                        end
                    end
                end
                % ROI data
                if ROIdataExist(k)
                    for jj = 1:size(ROIsvd,2)
                        vidMot = ROIsvd(:,jj)'; % get PC1 of the current ROI
                        for j = 1:numOffTrig
                            indStart = find(frameTimes > trigOnOffs(j,1) + timeWin(1),1,'first');
                            indEnd = indStart + numFramesToPlot - 1;
                            trialMvmtROI_temp_offset(j,:,jj) = vidMot(indStart:indEnd); % get the AV trial movements
                        end
                    end
                end


            end
        end
        trialMvmtAll_offset{k} = trialMvmtAll_temp_offset;
        trialMvmtROI_offset{k} = trialMvmtROI_temp_offset;
        durUntilPrevOffTrig_col_A{k} = durUntilPrevOffTrig_A;

        %% Get the single-trial FRs for each window - offset

        % Intialize the vectors for collecting the FRs
        FRs_baseline_temp_s = cell(numUnits,1); % For collecting the baseline FRs for the current recording's units
        FRs_off_temp_s = cell(numUnits,1); % For collecting the visual onset FRs for the current recording's units
        FRs_baseline_temp_l = cell(numUnits,1); % For collecting the baseline FRs for the current recording's units
        FRs_off_temp_l = cell(numUnits,1); % For collecting the visual onset FRs for the current recording's units
        FRs_baseline_temp_CMNlen = cell(numUnits,1); % For collecting the baseline FRs for the current recording's units
        FRs_off_temp_CMNlen = cell(numUnits,1); % For collecting the visual onset FRs for the current recording's units
        % Define the trials to collect FRs from
        curTrigsOn = trigOn_offset;
        for j = 1:numUnits
            %
            timeWin_temp = timeWin_baseline_off_s;
            FRs_baseline_temp_s{j} = getSingleTrialFRs_NP(spkTimes{j},curTrigsOn,timeWin_temp);
            % 
            timeWin_temp = timeWin_off_s;
            FRs_off_temp_s{j} = getSingleTrialFRs_NP(spkTimes{j},curTrigsOn,timeWin_temp);
            %
            timeWin_temp = timeWin_baseline_off_l;
            FRs_baseline_temp_l{j} = getSingleTrialFRs_NP(spkTimes{j},curTrigsOn,timeWin_temp);
            % 
            timeWin_temp = timeWin_off_l;
            FRs_off_temp_l{j} = getSingleTrialFRs_NP(spkTimes{j},curTrigsOn,timeWin_temp);
            %
            timeWin_temp = timeWin_baseline;
            FRs_baseline_temp_CMNlen{j} = getSingleTrialFRs_NP(spkTimes{j},curTrigsOn,timeWin_temp);
            % 
            timeWin_temp = timeWin_Vonset;
            FRs_off_temp_CMNlen{j} = getSingleTrialFRs_NP(spkTimes{j},curTrigsOn,timeWin_temp);
        end
        FRs_off_baseline_s{k} = FRs_baseline_temp_s;
        FRs_off_baseline_l{k} = FRs_baseline_temp_l;
        FRs_off_baseline_CMNlen{k} = FRs_baseline_temp_CMNlen;
        FRs_off_act_s{k} = FRs_off_temp_s;
        FRs_off_act_l{k} = FRs_off_temp_l;
        FRs_off_act_CMNlen{k} = FRs_off_temp_CMNlen;

        %% Get the PSTHs and rasters for offset trials

        % Initialize cell arrays for collecting PSTHs and rasters
        PSTHout_off_temp = cell(numUnits,1);
        raster_times_off_temp = cell(numUnits,1);
        raster_trInd_off_temp = cell(numUnits,1);
        for j = 1:numUnits
            % Off PSTH
            trigOnOffs = [trigOn_offset' trigOff_offset'];
            [PSTHout_off_temp{j},binCenters_off,~] = PSTHbasic_fixedBinSize_hColOut(spkTimes{j},trigOnOffs,-timeWin(1),PSTHbinSize,timeWin(2)-lenOffStim,lenOffStim);
            % Off raster
            [rasterSpikesTime,raster_trInd_off_temp{j}] = getRaster(trigOnOffs(:,1)+timeWin(1),trigOnOffs(:,1)+timeWin(2),spkTimes{j});
            raster_times_off_temp{j} = rasterSpikesTime + timeWin(1);
        end
        PSTHout_off{k} = cat(1,PSTHout_off_temp{:});
        raster_times_off{k} = raster_times_off_temp;
        raster_trInd_off{k} = raster_trInd_off_temp;

        %% Primary/secondary AC stats

        % For each channel, get:
        % 1) Tone onset latency
        % 2 P value to significantly elevated or not
        % 2) Tone FR BF to other amplitude ratio
        % 3) Tone FTC variance significance

        curRecFol = fullfile(localDataFol,curDate);
        tuningInd = strcmp(labBook_big.ExptDate,curDate) & strcmp(labBook_big.Multiblock_NP,rec) & ...
            strcmp(labBook_big.ExptType,'Tuning'); % position in lab book for the tuning protocol

        if exist(curRecFol,'dir') > 0 && sum(tuningInd) == 1 % if the folder and tuning protocol exist
            MU = load(fullfile(curRecFol,sprintf('MUspikeTimes_%s_probe%s.mat',rec,num2str(curProbe)))); % MU data
            tuningBlock = labBook_big.Block{tuningInd}; % block name for the tuning protoco.
            [secOn_T,secOff_T] = getBlockTime_stimMac(trda.blockOn_time,trda.blockOff_time,rec,tuningBlock,tuningBlock,false);
            indTrigWithinBlockTimes = trda.trigOn_time > secOn_T & trda.trigOn_time < secOff_T;
            load(fullfile(localDataFol,curDate,sprintf('exptParams_%s.mat',tuningBlock)),'frequencies')
            load(fullfile(localDataFol,curDate,sprintf('exptParams_%s.mat',tuningBlock)),'attenuations')
            trigOnOffs = [trda.trigOn_time(indTrigWithinBlockTimes)' trda.trigOn_time(indTrigWithinBlockTimes)'+lenStim];

            % Only take trials with a certain attenuation/s
            attenInd = ismember(attenuations,attenTake);
            frequencies_attenTake = frequencies(attenInd);
            trigOnOffs_attenTake = trigOnOffs(attenInd,:);

            curLats = nan(length(MU.spikeTimes),1);
            cPvals = nan(length(MU.spikeTimes),1);
            cPvals_twoSided = nan(length(MU.spikeTimes),1);
            overHalfHeight = false(length(MU.spikeTimes),1);
            varFrac = false(length(MU.spikeTimes),1);
            for curMU = 1:length(MU.spikeTimes)
                % Get PSTH - filtered and then z-scored
                curTimes = MU.spikeTimes{curMU};
                [curPSTH,binCenters_tone,~] = PSTHbasic_fixedBinSize_hColOut(curTimes,trigOnOffs_attenTake,0.2,binTimeLat,0.2,lenStim);
                % Filter PSTH
                curPSTH_fil = sgolayfilt(curPSTH,3,5);
                % Z-score PSTH
                baseInd = binCenters_tone > baseWin_T(1) & binCenters_tone < baseWin_T(2);
                zScored = ( curPSTH_fil - mean(curPSTH_fil(baseInd)) ) / std(curPSTH_fil(baseInd));
                % Up-sample filtered, z-scored PSTH for detecting the onset latency
                numPointsUS = length(binCenters_tone(1):binTimeLatUS:binCenters_tone(end));
                binCenters_interp1 = linspace(binCenters_tone(1),binCenters_tone(end),numPointsUS);
                zScored_interp1 = interp1(binCenters_tone,zScored,binCenters_interp1);
                actInd_US = binCenters_interp1 > actWin_T(1) & binCenters_interp1 < actWin_T(2);
                actTime_US = binCenters_interp1(actInd_US);
                % Get onset latency
                curLat = actTime_US(find(zScored_interp1(actInd_US)>=latTh,1,'first'));
                if ~isempty(curLat)
                    curLats(curMU) = curLat;
                end
                % Determine tone response p-value using spike counts
                spkCo_base = sum(curTimes >= trigOnOffs_attenTake(:,1)' - lenStim & curTimes < trigOnOffs_attenTake(:,1)',1);
                spkCo_act = sum(curTimes >= trigOnOffs_attenTake(:,1)' & curTimes < trigOnOffs_attenTake(:,2)',1);
                cPvals(curMU) = signrank(spkCo_base,spkCo_act,'tail','left');
                cPvals_twoSided(curMU) = signrank(spkCo_base,spkCo_act,'tail','both');

                spkCo_base_all = sum(curTimes >= trigOnOffs_attenTake(:,1)' - lenStim & curTimes < trigOnOffs_attenTake(:,1)');
                spkCo_act_all = sum(curTimes >= trigOnOffs_attenTake(:,1)' & curTimes < trigOnOffs_attenTake(:,2)');

                % Baseline subtracted spike counts for each freqency
                numFreq = length(unique(frequencies_attenTake));
                spkCo_base = zeros(1,numFreq);
                spkCo_act = zeros(1,numFreq);
                for curFreq = 1:numFreq
                    curFreqInd = frequencies_attenTake == curFreq;
                    spkCo_base(curFreq) = sum(spkCo_base_all(curFreqInd));
                    spkCo_act(curFreq) = sum(spkCo_act_all(curFreqInd));
                end
                spkCo_bs = spkCo_act - spkCo_base;
                
                % % Tuning BF half height test
                % [mx,mxInd] = max(spkCo_bs);
                % indNonBF = true(1,numFreq); indNonBF(mxInd(1)) = false;
                % overHalfHeight(curMU) = (mx/2) > mean(spkCo_bs(indNonBF));

                % Tuning variance test
                numIter = 1000;
                varReal = var(spkCo_bs);
                varShuffled = nan(numIter,1);
                numTrials = length(frequencies_attenTake);
                for curIter = 1:numIter
                    shuffledFreq = frequencies_attenTake(randperm(numTrials,numTrials));
                    spkCo_base_shuff = nan(1,numFreq);
                    spkCo_act_shuff = nan(1,numFreq);
                    for curFreq = 1:numFreq
                        curFreqInd = shuffledFreq == curFreq;
                        spkCo_base_shuff(curFreq) = sum(spkCo_base_all(curFreqInd));
                        spkCo_act_shuff(curFreq) = sum(spkCo_act_all(curFreqInd));
                    end
                    spkCo_bs_shuff = spkCo_act_shuff - spkCo_base_shuff;
                    varShuffled(curIter) = var(spkCo_bs_shuff);
                end
                varFrac(curMU) = sum(varShuffled > varReal) / numIter;

            end

            latencies_T{k} = curLats;
            pVals_T{k} = cPvals;
            pVals_twoSided_T{k} = cPvals_twoSided;
            % overHalfHeight_T{k} = overHalfHeight;
            varFrac_T{k} = varFrac;
            
            if exist(fracDepthFile,'file') > 0
                curDepthWithin = false(length(MU.spikeTimes),1);
                curDepthWithin(dp:sh) = true;
                depthWithin_T{k} = curDepthWithin;
            end

        end

        %% Ball movement

        fprintf('\nGetting trial-aligned ball movement...')

        % Load data
        ballMovementData = TDTbin2mat(TDTfileName, 'TYPE', {'scalars'}, 'STORE', 'UDP1');
        ballTime = ballMovementData.scalars.UDP1.ts;
        ballData = ballMovementData.scalars.UDP1.data;
        % Get the ball time and velocity
        cvx = ballData(1,:);
        cvy = ballData(1,:);
        dT = ballTime;
        dChannel = ballData(3,:);
        [ tsamp, v, ~ ] = mouse_velocity_twoOptiMice( cvx, cvy, dT, dChannel, frameSP );
        mx_val = prctile( v, 99 );
        mn_val = prctile( v, 0.1 );
        v( v > mx_val ) = mx_val;
        v( v < mn_val ) = mn_val;
        tsamp_c = tsamp(1:end-1) + median(diff(tsamp))/2; % time bin centers
        treadmill_vel = medfilt1(v,3);

        % Block on and off times for TDT
        blokData = TDTbin2mat(TDTfileName, 'TYPE', {'streams'}, 'STORE', 'Blok');
        startTimeTDT = blokData.streams.Blok.startTime;
        FsTDT = blokData.streams.Blok.fs;
        SPTDT = 1/FsTDT;
        blokTimes_TDT = (detectRisingEdge(blokData.streams.Blok.data,0.5)*SPTDT) + startTimeTDT;
        blocTimesWithOneAfter = [diff(blokTimes_TDT) < 0.02 false];
        blokOff = sum( [ blocTimesWithOneAfter ; [ false blocTimesWithOneAfter(1:end-1) ] ] ) == 1;
        blokOnTimes_TDT = blokTimes_TDT(~blokOff);
        blokOffTimes_TDT = blokTimes_TDT(blocTimesWithOneAfter);
        % Difference in time between NP and TDT recording start times
        recTimeSeconds = (str2double(rec(1:2)) * 3600) + (str2double(rec(4:5)) * 60) + str2double(rec(7:8));
        TDTrecTimeSeconds = (str2double(TDTname(1:2)) * 3600) + (str2double(TDTname(4:5)) * 60) + str2double(TDTname(7:8));
        timeDiff = TDTrecTimeSeconds - recTimeSeconds;
        % Trigger on and off time for VO protocol
        VOblokOn = blokOnTimes_TDT(abs((blokOnTimes_TDT + timeDiff) - secOn) == min(abs((blokOnTimes_TDT + timeDiff) - secOn))); % the TDT block on time closet to the NP VO block on time
        VOblokOff = blokOffTimes_TDT(abs((blokOffTimes_TDT + timeDiff) - secOff) == min(abs((blokOffTimes_TDT + timeDiff) - secOff))); % the TDT block off time closet to the NP VO block off time
        % Trigger on times
        trigData = TDTbin2mat(TDTfileName, 'TYPE', {'streams'}, 'STORE', 'Trig','T1',VOblokOn,'T2',VOblokOff);
        trigOnsTDT_VO = (detectRisingEdge(trigData.streams.Trig.data,0.5) * SPTDT) + VOblokOn; % the trigger on times, in TDT time
        % Find the trigger on times for the CMN and offset stimuli
        lastInSequence = [diff(trigOnsTDT_VO) > 1 true]; % Which triggers don't have a trigger immediately after them
        temp = sum( [ [ false ~lastInSequence(1:end-1) ] ; [ false false ~lastInSequence(1:end-2) ] ; ...
            [ false false false ~lastInSequence(1:end-3)  ] ] ) == 3; % find which triggers (lastInSequence) have three triggers immedidately before them
        offsetTrigInd_VO = false(1,length(lastInSequence));
        offsetTrigInd_VO(temp) = true; % the triggers with three triggers immediately before them are noise burst triggers
        trigOn_V_TDT = trigOnsTDT_VO(lastInSequence & ~offsetTrigInd_VO); % the TDT CMN trigger ons, in TDT time
        trigOn_off_TDT = trigOnsTDT_VO(lastInSequence & offsetTrigInd_VO); % the TDT noise bursts trigger ons, in TDT time

        if curProbe == 1
            treadmill_vel_store{i} = treadmill_vel;
            treadmill_time_store{i} = tsamp_c + (secOn - VOblokOn);
        end
        
        % Loop through each trigger on from the NP data.  Then find the
        % closest TDT trigger on, in actual time.  Then use this trigger
        % (in TDT time) to get the trial-aligned movement data
        minTimeDiff_V_temp = nan(numVtrig,1);
        minTimeDiff_A_temp = nan(length(trigOn_offset),1);
        if curProbe == 1
            % Collect trial-aligned ball movement for CMN trials
            trialMvmtBall_temp = nan(numVtrig,numFramesToPlot);
            timeFromBlockOn_TDT = trigOn_V_TDT - VOblokOn; % TDT trigger ons relative to TDT block on
            for j = 1:numVtrig
                % Get the TDT trigger closest to the current NP trigger (in actual time)
                timeFromBlockOn_NP = trigOn_VO(j) - secOn; % current NP trigger on (V trig) relative to the NP block on
                minTimeDiff_V_temp(j) = min(abs(timeFromBlockOn_TDT - timeFromBlockOn_NP));
                getTrig = abs(timeFromBlockOn_TDT - timeFromBlockOn_NP) == min(abs(timeFromBlockOn_TDT - timeFromBlockOn_NP)); % get the TDT trigger on time that has same duration between this trigger and the TDT block on, compared to the duration between the NP trigger and block on
                curTrigOn = trigOn_V_TDT(getTrig); % use this trigger, in TDT time
                % Store the trial-aligned ball movement
                indStart = find(tsamp_c > curTrigOn + timeWin(1),1,'first');
                indEnd = indStart + numFramesToPlot - 1;
                trialMvmtBall_temp(j,:) = treadmill_vel(indStart:indEnd);
            end
            % For noise burst trials
            trialMvmtBallOff_temp = nan(length(trigOn_offset),numFramesToPlot);
            timeFromBlockOn_TDT = trigOn_off_TDT - VOblokOn;
            for j = 1:length(trigOn_offset)
                % Get the TDT trigger closest to the current NP trigger (in actual time)
                timeFromBlockOn_NP = trigOn_offset(j) - secOn;
                minTimeDiff_A_temp(j) = min(abs(timeFromBlockOn_TDT - timeFromBlockOn_NP));
                getTrig = abs(timeFromBlockOn_TDT - timeFromBlockOn_NP) == min(abs(timeFromBlockOn_TDT - timeFromBlockOn_NP));
                curTrigOn = trigOn_off_TDT(getTrig);
                % Store the trial-aligned ball movement
                indStart = find(tsamp_c > curTrigOn + timeWin(1),1,'first');
                indEnd = indStart + numFramesToPlot - 1;
                trialMvmtBallOff_temp(j,:) = treadmill_vel(indStart:indEnd);
            end
        end
        trialMvmtBall{k} = trialMvmtBall_temp;
        trialMvmtBallOff{k} = trialMvmtBallOff_temp;
        minTimeDiff_V{k} = minTimeDiff_V_temp;
        minTimeDiff_A{k} = minTimeDiff_A_temp;
        fprintf('complete')

        %% end of loop

        fprintf('\nCompleted %s of %s recordings',num2str(k),num2str(numRec*2))

        k = k+1;
        if curProbe == 1
            recCounter = recCounter+1;
        end

    end

end

%% Store data and save

VO.V1_silencing_manDec = repelem(V1_silencing_manDec,numUnit);
VO.analyzeVC = repelem(analyzeVC,numUnit);
VO.signsOfBadExpression = repelem(signsOfBadExpression,numUnit);
VO.recNumber = repelem(recNumber,numUnit);
VO.mouseID = repelem(mouseID,numUnit);
VO.Ai32Mouse = repelem(Ai32Mouse,numUnit);
VO.wildTypeMouse = repelem(wildTypeMouse,numUnit);
VO.probeNum = repelem(probeNum,numUnit);
VO.num_Vtrig = repelem(num_Vtrig,numUnit);
VO.num_optoTrials = repelem(num_optoTrials,numUnit);
% Metadata
VO.unitIDs = cat(1,unitIDs{:});
VO.BSunits = cat(1,BSunits{:});
VO.NSunits = cat(1,NSunits{:});
VO.medUnits = cat(1,medUnits{:});
VO.fracDepths = cat(1,fracDepths{:});
VO.channels = cat(1,channelsCol{:});
VO.TPWs = cat(1,TPWs{:});
VO.troughHeightRatios = cat(1,troughHeightRatios{:});
VO.endSlopes = cat(1,endSlopes{:});
VO.numEachWay = numEachWay;
VO.waveforms_mean = cat(1,waveforms_mean{:});
VO.waveforms_SD = cat(1,waveforms_SD{:});
VO.ISIviolPer = cat(1,ISIviolPer{:});
VO.waveforms_tAx = waveFormTax;
% PSTHs and rasters
VO.PSTHout_V = cat(1,PSTHout_V{:});
VO.raster_times_V = cat(1,raster_times_V{:});
VO.raster_trInd_V = cat(1,raster_trInd_V{:});
VO.PSTHout_VO = cat(1,PSTHout_VO{:});
VO.raster_times_VO = cat(1,raster_times_VO{:});
VO.raster_trInd_VO = cat(1,raster_trInd_VO{:});
VO.PSTHout_off = cat(1,PSTHout_off{:});
VO.raster_times_off = cat(1,raster_times_off{:});
VO.raster_trInd_off = cat(1,raster_trInd_off{:});
% FRs (and opto trial index)
VO.FRs_baseline = cat(1,FRs_baseline{:});
VO.FRs_Vonset = cat(1,FRs_Vonset{:});
VO.FRs_trialType_optoTrials = cat(1,FRs_trialType_optoTrials{:});
% FRs - offset
VO.FRs_off_baseline_s = cat(1,FRs_off_baseline_s{:});
VO.FRs_off_baseline_l = cat(1,FRs_off_baseline_l{:});
VO.FRs_off_baseline_CMNlen = cat(1,FRs_off_baseline_CMNlen{:});
VO.FRs_off_act_s = cat(1,FRs_off_act_s{:});
VO.FRs_off_act_l = cat(1,FRs_off_act_l{:});
VO.FRs_off_act_CMNlen = cat(1,FRs_off_act_CMNlen{:});
% Movement stuff
VO.trialMvmtAll = repelem(trialMvmtAll,numUnit);
VO.trialMvmtROI = repelem(trialMvmtROI,numUnit);
VO.movieStatus = repelem(movieStatus,numUnit);
VO.ROIdataExist = repelem(ROIdataExist,numUnit);
VO.trialMvmtMean_VO = repelem(trialMvmtMean_VO,numUnit);
VO.trialMvmtMean_VO_baseline = repelem(trialMvmtMean_VO_baseline,numUnit);
VO.durUntilPrevOffTrig_V = repelem(durUntilPrevOffTrig_col_V,numUnit);
VO.durUntilPrevOffTrig_A = repelem(durUntilPrevOffTrig_col_A,numUnit);
% Bin centers
VO.binCenters_faceMvmts = binCenters_faceMvmts;
VO.binCenters_VO = binCenters_VO;
VO.binCenters_off = binCenters_off;
% FracDepth details for each recording
VO.fracDepthRecData.depths = fracDepthDataCol;
VO.fracDepthRecData.recID = fracDepth_recID;
VO.fracDepthRecData.probe = fracDepth_probeCol;
% Offset stuff
VO.trialMvmtAll_offset = repelem(trialMvmtAll_offset,numUnit);
VO.trialMvmtROI_offset = repelem(trialMvmtROI_offset,numUnit);
VO.offSpikingResponses = cat(1,offSpikingResponses{:});
% Windows and other timing things
VO.timeWin = timeWin;
VO.PSTHbinSize = PSTHbinSize;
VO.timeWin_Vonset = timeWin_Vonset;
VO.timeWin_baseline = timeWin_baseline;
VO.TPWthresh = TPWthresh;
VO.CMNdur = CMNdur;
VO.lenOffStim = lenOffStim;
VO.frameSP = frameSP;
VO.timeWin_baseline_off_s = timeWin_baseline_off_s;
VO.timeWin_baseline_off_l = timeWin_baseline_off_l;
VO.timeWin_off_s = timeWin_off_s;
VO.timeWin_off_l = timeWin_off_l;
% Eigen value stuff
VO.eigen.eigenValues_varExp = eigenValues_varExp;
VO.eigen.eigenValues_recID = eigenValues_recID;
% MU tone stuff for primary/seconary designation
VO.toneLatencies_MU = latencies_T;
VO.tonePvals_MU = pVals_T;
VO.tonePvalsTwoSided_MU = pVals_twoSided_T;
% VO.toneOverhalfHeight_MU = overHalfHeight_T;
VO.toneVarFrac_MU = varFrac_T;
VO.toneDepthWithin_MU = depthWithin_T;
VO.toneLatSDthreshold = latTh;
% Opto voltage
VO.maxOptoVol = maxOptoVol;
% Store frame times and FM PC1 (whole trace)
VO.store_FM_PC1 = FM_PC1_store;
VO.store_FM_ROI = FM_ROI_store;
VO.store_frameTimes = frameTimes_store;
% Store V and offset trigger on times
VO.store_trigOn_VO = trigOn_store_VO;
VO.store_trigOn_off = trigOn_store_off;
% Ball movement
VO.trialMvmtBall = repelem(trialMvmtBall,numUnit);
VO.trialMvmtBallOff = repelem(trialMvmtBallOff,numUnit);
VO.minTimeDiff_V = repelem(minTimeDiff_V,numUnit);
VO.minTimeDiff_A = repelem(minTimeDiff_A,numUnit);
VO.store_ball_vel = treadmill_vel_store;
VO.store_ball_time = treadmill_time_store;

saveName = fullfile(saveFol,sprintf('VOstats_NP_%s.mat',datetime('now','Format','yyyy-MM-dd_HH-mm-ss')));

save(saveName,'-struct','VO','-v7.3');






