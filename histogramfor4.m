clf
clear
clc
numguass = [1,1,1,1];
ccolor = ['brmg'];
shape = ['----'];
tt = {'Collected data from Heathy neutrophil'; ...
      'Collected data from TNF-Alpha treated neutrophil'; ... 
      'Collected data from Glucose treated neutrophil'; ...
      'Collected data from PMA treated neutrophil'}
file = {'G:\drhou\HANWEI\new\03262015 60um H 400um W 1cm channel length w e-SEL+BSa blocking 0.5%BSA\healthy DFF_O2\Time Series.mat'; ...
    'G:\drhou\HANWEI\new\03202015 TNFa 10ngml\Time Series.mat'; ...
    'G:\drhou\HANWEI\new\03272015 60um H 400um W 1cm channel length w e-SEL+BSa blocking 0.5%BSA\30mM glucose DFF_O2\Time Series.mat'; ...
    'G:\drhou\HANWEI\new\03262015 60um H 400um W 1cm channel length w e-SEL+BSa blocking 0.5%BSA\PMA DFF_O2\Time Series-9.mat'
    };
for k = 1 : length(file)
    histdata = [];
    load(file{k});
    savename = fpname(1:max(findstr(fpname,'\')))
    if k == 2
        buffer = pm.vidvel;
        load('G:\drhou\HANWEI\new\03272015 60um H 400um W 1cm channel length w e-SEL+BSa blocking 0.5%BSA\TNFa 10ngml DFF_O2\Time Series-9.mat');
        pm.vidvel = [pm.vidvel buffer];
    end
    for i = 1:size(pm.vidvel,2)
        for j = 1:size(pm.vidvel{i},2)
            histdata = [histdata mean(pm.vidvel{i}{j})];
        end
    end
    
    hplot.raw{k} = histdata;
    
    [histFreq, histXout] = hist(histdata,[0:pm.step:pm.lng]);
    
%     %plot histogram in the same figure
%     figure(1)
%     subplot(2,length(file),k)
%     hb = bar(histXout, histFreq/sum(histFreq)*100,1,pm.color);
%     hold on
%     xlim([0,pm.lng])
%     if k ~=4
%         ylim([0,20])
%     end
%     xlabel('Rolling cell velocity (um/s)');
%     ylabel('Frequency (Normalized)');
    
    %plot distribution curve separately
    for q = 1:2
        figure(q+1)
        if q == 1
            clf(figure(2))
            hb = bar(histXout, histFreq/sum(histFreq)*100,1,ccolor(k));
            hold on
            
            PD = fitdist(histdata', 'normal');
            gcurve = gmdistribution.fit(histdata',numguass(k));
            gc{k}.curve = gcurve;
            %title(tt{k})
        end
        plot(histXout, pdf(gcurve, histXout')*100,[ccolor(k) '-' shape(k)]);
        hold on
        xlim([0,pm.lng])
        xlabel('Rolling cell velocity (um/s)');
        ylabel('Frequency (Normalized)');
    end

    pause
%     %create side by side histogram after collect all data
%     if k == length(file)
%         hplot.freq = histFreq;
%         hplot.value = [gdata.value/sum(gdata.value)*100;tdata.value/sum(tdata.value)*100;hdata.value/sum(hdata.value)*100;pdata.value/sum(pdata.value)*100]';
%         figure(4)
%         hb = bar(hplot.freq,hplot.value);
%         xlim([-2,pm.lng])
%         ylim([0,30])
%         xlabel('Rolling cell velocity (um/s)');
%         ylabel('Frequency (Normalized)');
%         hold on
%         
%         for pp = 1:4
%             set(hb(pp),'FaceColor',ccolor(pp))
%             plot(hplot.freq, pdf(gc{pp}.curve, hplot.freq')*100,[ccolor(pp) '--'],'LineWidth', 2);
%             hold on
%         end
%     end
end
