function zScore_scatSize = getZScoreScatterSizes(zScores,zScoreScatterSizeRange,minMax)

zScoreAdj = (zScores - minMax(1)) / diff(minMax);
zScore_scatSize = (zScoreAdj * diff(zScoreScatterSizeRange)) + zScoreScatterSizeRange(1);
