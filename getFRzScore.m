function FRzScore = getFRzScore(FRs_one,FRs_two,inds)

FRzScore = nan(length(FRs_one),1);
for i = 1:length(FRs_one)
    FRzScore(i) = ( mean(FRs_two{i}(inds{i})) - mean(FRs_one{i}(inds{i})) ) / std(FRs_one{i}(inds{i})) ;
end

