function trajsum(pm,savename,vrt)
%trajectorysum

if nargin <3
    vrt = 1
end
    

%plot cell trajectory in graph
for k = 1:vrt
    figure(9)
    clf
    limx = 60/k;
    for i = 1:size(pm.distdata,2)
        for j = 1:size(pm.distdata{i},2)
            plot(pm.distdata{i}{j}(:,1),pm.distdata{i}{j}(:,2),'color',pm.color, 'LineWidth', 1)
            hold on
        end
    end
    ylim([-limx limx])
    xlim([-250 0])
    xlabel('X-displacement (um)');
    ylabel('Y-displacement (um)');
    
    figure(10)
    clf
    
    for i = 1:size(pm.dispdata,2)
        for j = 1:size(pm.dispdata{i},2)
            plot([0,pm.dispdata{i}{j}(:,1)],[0,pm.dispdata{i}{j}(:,2)],'color',pm.color, 'LineWidth', 1)
            hold on
        end
    end
    
    ylim([-limx limx])
    xlim([-250 0])
    
    xlabel('X-displacement (um)');
    ylabel('Y-displacement (um)');
    fldname = [savename(1:max(findstr(savename,'\'))),'sum\']
    mkdir(fldname)
    saveas(figure(9),[fldname,'trajectorysum_',num2str(limx),'.jpg'])
    saveas(figure(10),[fldname,'\trajectorysum_disp_',num2str(limx),'.jpg'])
end
disp('finish lah')
end