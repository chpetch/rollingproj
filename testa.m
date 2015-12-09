auc_diabetic.ch = [];
auc_diabetic.norm = [];
auc_healthy.ch = [];
auc_healthy.norm = [];

%patient#9 video 9
figure(2)
for i = 1:1:4
    auc_ch_9 = [];
    auc_9 = [];
    delcellarray = [];
    if i == 1
        load('C:\Users\Chayakorn\OneDrive for Business\drhou\HOUHANWEI TTSH\TTSH clinical testing (April 13 onwards)\04172015 WBCs-0009 n 0010\0009\60x - 1million mL\Time Series-7.mat')
        delcellarray = [20];
        result_edited = delcell(result,delcellarray);
        result = result_edited;
    elseif i == 2
        load('C:\Users\Chayakorn\OneDrive for Business\drhou\HOUHANWEI TTSH\TTSH clinical testing (April 13 onwards)\04172015 WBCs-0009 n 0010\0009\60x - 1million mL\Time Series-8.mat')
        delcellarray = [];
        result_edited = delcell(result,delcellarray);
        result = result_edited;
    elseif i == 3
        load('C:\Users\Chayakorn\OneDrive for Business\drhou\HOUHANWEI TTSH\TTSH clinical testing (April 13 onwards)\04172015 WBCs-0009 n 0010\0009\60x - 1million mL\Time Series-9.mat')
        delcellarray = [];
        result_edited = delcell(result,delcellarray);
        result = result_edited;
    elseif i == 4
        load('C:\Users\Chayakorn\OneDrive for Business\drhou\HOUHANWEI TTSH\TTSH clinical testing (April 13 onwards)\04172015 WBCs-0009 n 0010\0009\60x - 1million mL\Time Series-10.mat')
        delcellarray = [13,14,34];
        result_edited = delcell(result,delcellarray);
        result = result_edited;
    end
    for j = 1:size(result.cpos_new,2)
        clear dataplot
        dataplot = bsxfun(@minus,result.cpos_new{j},result.cpos_new{j}(1,:));
        cdist{j} = dataplot;
        plot(abs(dataplot(:,1)),dataplot(:,2),'color','k', 'LineWidth', 1)
        hold on
        %xlim([-300,0])
        ylim([-60,60])
        auc_ch_9(j)=linerotation(dataplot(:,1),dataplot(:,2),0)/abs(dataplot(end,1));
        %title(['Video: ',num2str(i),' Cell: ',num2str(j)])
        auc_9(j) = abs(trapz(dataplot(:,1),dataplot(:,2))/abs(dataplot(end,1)));
    end
    auc_ch_9 = auc_ch_9';
    auc_9 = auc_9';
    disp(['diabetic',num2str(i),'_',num2str(length(auc_9))])
    auc_diabetic.ch = [auc_diabetic.ch;auc_ch_9];
    auc_diabetic.norm = [auc_diabetic.norm;auc_9];
end

clear delcellarray
%figure(3)
for i =1:2
    %patient#28 video 7
    auc_ch_28 = [];
    auc_28 = [];
    delcellarray = [];
    if i == 1
        load('C:\Users\Chayakorn\OneDrive for Business\drhou\09042015 WBCs-0028 (ctrl)\60x 1 million per mL\Time Series-9.mat')
        delcellarray = [4,12,20,30];
        result_edited = delcell(result,delcellarray);
        result = result_edited;
    elseif i == 2
        load('C:\Users\Chayakorn\OneDrive for Business\drhou\09042015 WBCs-0028 (ctrl)\60x 1 million per mL\Time Series-10.mat')
        delcellarray = [6,18,22,24,30];
        result_edited = delcell(result,delcellarray);
        result = result_edited;
    end
    for j = 1:size(result.cpos_new,2)
        clear dataplot
        dataplot = bsxfun(@minus,result.cpos_new{j},result.cpos_new{j}(1,:));
        cdist{j} = dataplot;
        plot(dataplot(:,1),dataplot(:,2),'color','b', 'LineWidth', 1)
        hold on
%         xlim([-460,460])
%         ylim([-60,60])
        auc_ch_28(j)=linerotation(dataplot(:,1),dataplot(:,2),0)/abs(dataplot(end,1));
% auc_ch_28(j)=linerotation(dataplot(:,1),dataplot(:,2),1)/abs(dataplot(end,1));
% title(['Video: ',num2str(i),' Cell: ',num2str(j)])        
% xlim([0,500])
% ylim([-60,60])

        auc_28(j) = abs(trapz(dataplot(:,1),dataplot(:,2))/abs(dataplot(end,1)));
    end
    auc_ch_28 = auc_ch_28';
    auc_28 = auc_28';
    disp(['healthy',num2str(i),'_',num2str(length(auc_28))])
    auc_healthy.ch = [auc_healthy.ch;auc_ch_28];
    auc_healthy.norm = [auc_healthy.norm;auc_28];
end