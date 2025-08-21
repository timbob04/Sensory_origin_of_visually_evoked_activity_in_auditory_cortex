function [pVals_Mann, pVals_perm] = getLowVsHighMvmtPvals(FRs_one, FRs_two, trialsToAvoid, referenceGroup)

% numIter = 10000;

referenceGroup = cellfun(@(x,y) ~x & y,  trialsToAvoid, referenceGroup,'UniformOutput',false);
targetGroup    = cellfun(@(x,y) ~x & ~y, trialsToAvoid, referenceGroup,'UniformOutput',false);

fprintf('\nGetting FR change rank-sum & permutation p-vals for mvmt vs no mvmt trials, %s cells. Completed: ',num2str(length(FRs_one)))

pVals_Mann  = nan(length(FRs_one),1);
pVals_perm  = nan(length(FRs_one),1);

for i = 1:length(FRs_one)

    % trial-level FR changes
    dFR = FRs_two{i} - FRs_one{i};

    dFR_noMvmt = dFR(referenceGroup{i});
    dFR_mvmt   = dFR(targetGroup{i});

    % Skip if either group empty
    if isempty(dFR_noMvmt) || isempty(dFR_mvmt)
        counterStr(i);
        continue
    end

    % Mann–Whitney (two-sided)
    pVals_Mann(i) = ranksum(dFR_mvmt, dFR_noMvmt, 'tail','both');

%     % Observed median difference (keep sign convention consistent)
%     obs = median(dFR_noMvmt) - median(dFR_mvmt);
% 
%     % Label-permutation null
%     pooledData = [dFR_noMvmt; dFR_mvmt];
%     nNoMvmt    = numel(dFR_noMvmt);
%     nTotal     = numel(pooledData);
% 
%     null = nan(numIter,1);
%     for j = 1:numIter
%         indNoMvmt = false(nTotal,1);
%         indNoMvmt(randperm(nTotal, nNoMvmt)) = true;
% 
%         dFR_noMvmt_iter = pooledData(indNoMvmt);
%         dFR_mvmt_iter   = pooledData(~indNoMvmt);
% 
%         null(j) = median(dFR_noMvmt_iter) - median(dFR_mvmt_iter);
%     end
% 
%     % Two-sided permutation p with +1 correction
%     pVals_perm(i) = (sum(abs(null) >= abs(obs)) + 1) / (numIter + 1);

    counterStr(i);
    
end
