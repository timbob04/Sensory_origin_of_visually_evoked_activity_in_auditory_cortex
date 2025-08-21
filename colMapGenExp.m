function colMapOut = colMapGenExp(lowCol,midCol,topCol,numCol)

slope = 0.4;

x = 1:numCol;
y = 2.^(slope*x);
y_norm = y / diff([y(1) y(end)]);
y_norm_bs = y_norm - y_norm(1);

lowCols = nan(numCol,3);
for i = 1:3
    minVal = midCol(i);
    maxVal = lowCol(i);
    lowCols(:,i) = y_norm_bs * diff([minVal maxVal]) + minVal;
end
    
topCols = nan(numCol,3);
for i = 1:3
    minVal = midCol(i);
    maxVal = topCol(i);
    topCols(:,i) = y_norm_bs * diff([minVal maxVal]) + minVal;
end

colMapOut = [ flipud(lowCols) ; topCols(2:end,:) ];


    
    

