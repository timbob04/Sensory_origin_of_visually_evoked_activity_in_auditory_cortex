function [sortedImsc,orderVecOut] = plotImagesc(imscData,indUseVec,orderVec,ColAxisLims,tAx,plotBetween,xTicks,colMap)

[sortedImsc,sortNum] = sortrows([orderVec(indUseVec) imscData(indUseVec,:)]); % sorted imagesc data
sortedImsc = sortedImsc(:,2:end);

curCells = find(indUseVec);
orderVecOut = curCells(sortNum);

imagesc(sortedImsc)
colormap(colMap)
caxis(ColAxisLims)

xStart = find( abs(tAx - plotBetween(1)) == min(abs(tAx - plotBetween(1))) ,1,'first');
xEnd = find( abs(tAx - plotBetween(2)) == min(abs(tAx - plotBetween(2))) ,1,'first');

xTickPos = nan(1,length(xTicks));
for i = 1:length(xTicks)
    xTickPos(i) = find( abs(tAx - xTicks(i)) == min(abs(tAx - xTicks(i))) ,1,'first');
end

set(gca,'XLim',[xStart xEnd],'xtick',xTickPos,'xticklabel',xTicks,'tickdir','out','linewidth',1,'box','off','fontsize',14,'YDir','normal')