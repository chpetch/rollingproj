%histogram sum
function  output = histsum(histdata,step,lng,savename,color,md,pc)

if nargin < 6 
    md = 1;
else nargin < 7
    pc = 0;
end

%histdata = [histdata zeros(1,37)];

[histFreq, histXout] = hist(histdata,[0:step:lng]);
output.value = histFreq;
output.freq = histXout;
if pc == 0
    figure(6)
    clf
    if md == 1 || md == 3
        bar(histXout, histFreq/sum(histFreq)*100,1,color);
        hold on
    end
    if md == 2 || md == 3
        PD = fitdist(histdata', 'normal');
        plot(histXout, pdf(PD, histXout)*100,[color ':'],'LineWidth', 3);
    end
    xlim([0,lng])
    xlabel('Rolling cell velocity (um/s)');
    ylabel('Frequency (Normalized)');
end

sv = step/2;
for i = 1 : length(histXout)
    if i == 1
        bandwidth{i} = ['0 - ',num2str(sv)];
    elseif i == length(histXout)
        bandwidth{i} = [num2str(sv),' - ','more than ',num2str(lng)];
    else
        bandwidth{i} = [num2str(sv),' - ',num2str(sv+step)];
        sv = sv + step;
    end
end
tb = {'Bin Center (X)','Frequency (Y)', 'Percentage'};
tbb = [histXout',histFreq',(histFreq/sum(histFreq)*100)'];
tb = [tb;num2cell(tbb)];
tb = [['Binwidth' ; bandwidth'] tb];

if md == 1
    mkdir([savename])
    xlswrite([savename,'\histdatasum.xls'], tb, 'Histogram', 'A1');
    saveas(figure(6),[savename,'\histogramsum.jpg'])
elseif md == 2
    mkdir([savename])
    saveas(figure(6),[savename,'\curvesum.jpg'])
elseif md == 3
    mkdir([savename])
    xlswrite([savename,'\histdatasum.xls'], tb, 'Histogram', 'A1');
    saveas(figure(6),[savename,'histogramcurve.jpg'])    
end

end