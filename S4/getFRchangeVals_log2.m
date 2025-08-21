function FRchange_log2 = getFRchangeVals_log2(FRs_one,FRs_two,inds,lims)

FRchange = nan(length(FRs_one),1);
for i = 1:length(FRs_one)
    FRchange(i) = mean(FRs_two{i}(inds{i})) / mean(FRs_one{i}(inds{i}));
end

FRchange_log2 = log2(FRchange);

FRchange_log2(FRchange_log2<lims(1)) = lims(1);
FRchange_log2(FRchange_log2>lims(2)) = lims(2);