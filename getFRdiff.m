function FRdiff = getFRdiff(FRs_one,FRs_two,inds)

FRdiff = nan(length(FRs_one),1);
for i = 1:length(FRs_one)
    FRdiff(i) = mean(FRs_two{i}(inds{i})) - mean(FRs_one{i}(inds{i}));
end

