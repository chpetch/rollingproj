clf
clear
clc
numguass = [2,3,1,1];
ccolor = ['mrbg'];
for k = 1 : 4
    histdata = [];
    if k == 1
        load('C:\Users\nwbi\Google Drive\NTU\Progress\drhou\HANWEI\new\03272015 60um H 400um W 1cm channel length w e-SEL+BSa blocking 0.5%BSA\30mM glucose DFF_O2\Time Series.mat')
        pm.color = 'm';
        savename = fpname(1:max(findstr(fpname,'\')))
        gdata = histsum(pm.vidvel,pm.step,pm.lng,savename,pm.color,3);
        %trajsum(pm,savename,3)
    elseif k == 2
        load('C:\Users\nwbi\Google Drive\NTU\Progress\drhou\HOUHANWEI TTSH\03202015 TNFa 10ngml\Time Series.mat')
        adddata = pm.vidvel;
        load('C:\Users\nwbi\Google Drive\NTU\Progress\drhou\HANWEI\new\03272015 60um H 400um W 1cm channel length w e-SEL+BSa blocking 0.5%BSA\TNFa 10ngml DFF_O2\Time Series-9.mat')
        pm.vidvel = [pm.vidvel adddata]
        pm.color = 'r';
        savename = fpname(1:max(findstr(fpname,'\')))
        tdata = histsum(pm.vidvel,pm.step,pm.lng,savename,pm.color,3)
        %trajsum(pm,savename,3)
    elseif k == 3
        load('C:\Users\nwbi\Google Drive\NTU\Progress\drhou\HANWEI\new\03262015 60um H 400um W 1cm channel length w e-SEL+BSa blocking 0.5%BSA\healthy DFF_O2\Time Series.mat')
        % load('C:\Users\nwbi\Google Drive\NTU\Progress\drhou\HANWEI\new\03272015 60um H 400um W 1cm channel length w e-SEL+BSa blocking 0.5%BSA\healthy DFF_O2\Time Series-7.mat')
        
        pm.color = 'b';
        savename = fpname(1:max(findstr(fpname,'\')))
        hdata = histsum(pm.vidvel,pm.step,pm.lng,savename,pm.color,3)
        %trajsum(pm,savename,3)
    elseif k==4
        load('C:\Users\nwbi\Google Drive\NTU\Progress\drhou\HANWEI\new\03262015 60um H 400um W 1cm channel length w e-SEL+BSa blocking 0.5%BSA\PMA DFF_O2\Time Series-8.mat')
        pm.color = 'g';
        savename = fpname(1:max(findstr(fpname,'\')))
        pdata = histsum(pm.vidvel,pm.step,pm.lng,savename,pm.color,3)
    end
    for i = 1:size(pm.vidvel,2)
        for j = 1:size(pm.vidvel{i},2)
            histdata = [histdata mean(pm.vidvel{i}{j})];
        end
    end
    
    hplot.raw{k} = histdata;
    
    [histFreq, histXout] = hist(histdata,[0:pm.step:pm.lng]);
    
    %plot histogram in the same figure
    figure(1)
    subplot(2,2,k)
    hb = bar(histXout, histFreq/sum(histFreq)*100,1,pm.color);
    hold on
    xlim([0,pm.lng])
    if k ~=4
        ylim([0,20])
    end
    xlabel('Rolling cell velocity (um/s)');
    ylabel('Frequency (Normalized)');
    
    %plot distribution curve separately
    for q = 1:2
        figure(q+1)
        if q == 1
            clf(figure(2))
            hb = bar(histXout, histFreq/sum(histFreq)*100,1,pm.color);
            hold on
            
            PD = fitdist(histdata', 'normal');
            gcurve = gmdistribution.fit(histdata',numguass(k));
            gc{k}.curve = gcurve;            
        end
        plot(histXout, pdf(gcurve, histXout')*100,[ccolor(k) '--'],'LineWidth', 2);
        hold on
        xlim([0,pm.lng])
        xlabel('Rolling cell velocity (um/s)');
        ylabel('Frequency (Normalized)');
        if q == 2 && k == 4
            legend('Glucose','TNF-Alpha','Healthy','PMA')
        end
    end
    pause
    %create side by side histogram after collect all data
    if k == 4
        hplot.freq = gdata.freq;
        hplot.value = [gdata.value/sum(gdata.value)*100;tdata.value/sum(tdata.value)*100;hdata.value/sum(hdata.value)*100;pdata.value/sum(pdata.value)*100]';
        figure(4)
        hb = bar(hplot.freq,hplot.value);
        xlim([-2,pm.lng])
        ylim([0,30])
        xlabel('Rolling cell velocity (um/s)');
        ylabel('Frequency (Normalized)');
        hold on
        
        for pp = 1:4
            set(hb(pp),'FaceColor',ccolor(pp))
            plot(hplot.freq, pdf(gc{pp}.curve, hplot.freq')*100,[ccolor(pp) '--'],'LineWidth', 2);
            hold on
        end
        legend('Glucose','TNF-Alpha','Healthy','PMA')
    end
end
