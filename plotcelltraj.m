function cdist = plotcelltraj(cposition,xrange,yrange,color)
dataplot = [];
for j = 1:size(cposition,2)
    clear dataplot
    dataplot = bsxfun(@minus,cposition{j},cposition{j}(1,:));
    cdist{j} = dataplot;
    plot(dataplot(:,1),dataplot(:,2),'color',color, 'LineWidth', 1)
    hold on
    xlim(xrange)
    ylim(yrange)
end
end