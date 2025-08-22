function FRchange = getFRchangeVals(FRs_one,FRs_two,inds)

FRchange = nan(length(FRs_one),1);
for i = 1:length(FRs_one)
    FRchange(i) = mean(FRs_two{i}(inds{i})) / mean(FRs_one{i}(inds{i}));
end