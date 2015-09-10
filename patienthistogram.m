clf
clear
clc
numguass = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
ccolor = ['mmmgggbbccrrrryryrykk'];
shape = ['+o*.xsd^v><+o*.xsd--'];
tt = {'Collected data from Patient 005'; ...
      'Collected data from Patient 006'; ... 
      'Collected data from Patient 007'; ...
      'Collected data from Patient 008'; ...
      'Collected data from Patient 009'; ...
      'Collected data from Patient 010'; ...
      'Collected data from Patient 011'; ...
      'Collected data from Patient 012'; ...
      'Collected data from Patient 013'; ...
      'Collected data from Patient 014'; ...
      'Collected data from Patient 015'; ...
      'Collected data from Patient 016'; ...
      'Collected data from Patient 017'; ...
      'Collected data from Patient 018'; ...
      'Collected data from Patient 019'; ...
      'Collected data from Patient 019 fingerprick'; ...
      'Collected data from Patient 020'; ...
      'Collected data from Patient 020 fingerprick'; ...
      'Collected data from Patient 021'; ...
      'Collected data from CW'; ...
      'Collected data from XW'}
file = {'C:\Users\Chayakorn\OneDrive for Business\drhou\HOUHANWEI TTSH\TTSH clinical testing (April 13 onwards)\04132015 WBCs-0005 n 0006\0005\20x-1million mL\Time Series-9.mat'; ...
    'C:\Users\Chayakorn\OneDrive for Business\drhou\HOUHANWEI TTSH\TTSH clinical testing (April 13 onwards)\04132015 WBCs-0005 n 0006\0006\20x 1 million mL\Time Series.mat'; ...
    'C:\Users\Chayakorn\OneDrive for Business\drhou\HOUHANWEI TTSH\TTSH clinical testing (April 13 onwards)\04162015 WBCs-0007 n 0008\0007\20x - 1million mL\Time Series.mat'; ...
    'C:\Users\Chayakorn\OneDrive for Business\drhou\HOUHANWEI TTSH\TTSH clinical testing (April 13 onwards)\04162015 WBCs-0007 n 0008\0008\20x- 1million mL\Time Series.mat'; ...
    'C:\Users\Chayakorn\OneDrive for Business\drhou\HOUHANWEI TTSH\TTSH clinical testing (April 13 onwards)\04172015 WBCs-0009 n 0010\0009\20x - 1million mL\Time Series.mat'; ...
    'C:\Users\Chayakorn\OneDrive for Business\drhou\HOUHANWEI TTSH\TTSH clinical testing (April 13 onwards)\04172015 WBCs-0009 n 0010\0010\20x - 1million mL\Time Series.mat'; ...
    'C:\Users\Chayakorn\OneDrive for Business\drhou\HOUHANWEI TTSH\TTSH clinical testing (April 13 onwards)\04272015 WBCs- 0011\1 million per mL 20x\Time Series.mat'; ...
    'C:\Users\Chayakorn\OneDrive for Business\drhou\HOUHANWEI TTSH\TTSH clinical testing (April 13 onwards)\04272015 WBCs- 0012 n 0013\0012\1 million per mL 20x\Time Series.mat'; ...
    'C:\Users\Chayakorn\OneDrive for Business\drhou\05192015 WBCs-0013 ctrl\20x\Time Series-6.mat'; ...
    'C:\Users\Chayakorn\OneDrive for Business\drhou\05212015 WBCs- 0014\20x\Time Series.mat'; ...
    'C:\Users\Chayakorn\OneDrive for Business\drhou\05272015 WBCs- 0015\20x- 1 million per mL\Time Series.mat'; ...
    'C:\Users\Chayakorn\OneDrive for Business\drhou\05292015 WBCs-0016 and 0017\0016\20x 1 million per mL\Time Series.mat'; ...
    'C:\Users\Chayakorn\OneDrive for Business\drhou\05292015 WBCs-0016 and 0017\0017\20x-1 million per mL\Time Series.mat'; ...
    'C:\Users\Chayakorn\OneDrive for Business\drhou\06032015 WBCs-0018\20x-1 million per mL\Time Series.mat'; ...
    'C:\Users\Chayakorn\OneDrive for Business\drhou\06262015 WBCs-0019 n 0020 (ctrls)\0019\20x-1 million per ml\Time Series.mat'; ...
    'C:\Users\Chayakorn\OneDrive for Business\drhou\06262015 WBCs-0019 n 0020 (ctrls)\0019-fingerprick\20x-1 million per mL\Time Series-8.mat'; ...
    'C:\Users\Chayakorn\OneDrive for Business\drhou\06262015 WBCs-0019 n 0020 (ctrls)\0020\20x- 1 million per mL\Time Series.mat'; ...
    'C:\Users\Chayakorn\OneDrive for Business\drhou\06262015 WBCs-0019 n 0020 (ctrls)\0020-fingerprick\20x-1 million per mL\Time Series.mat'; ...
    'C:\Users\Chayakorn\OneDrive for Business\drhou\06302015 WBCs-0021 (ctrl)\20x-1 million per mL\Time Series.mat'; ...
    'C:\Users\Chayakorn\OneDrive for Business\drhou\HOUHANWEI TTSH\TTSH clinical testing (April 13 onwards)\04212015 CW blood\20x\Time Series.mat'; ...
    'C:\Users\Chayakorn\OneDrive for Business\drhou\HOUHANWEI TTSH\TTSH clinical testing (April 13 onwards)\04222015 XW blood\20x\Time Series.mat'};
for k = 1 : length(file)
    histdata = [];
    load(file{k});
    savename = fpname(1:max(findstr(fpname,'\')))
    
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
            title(tt{k})
        end
        plot(histXout, pdf(gcurve, histXout')*100,[ccolor(k) '-' shape(k)]);
        hold on
        xlim([0,pm.lng])
        xlabel('Rolling cell velocity (um/s)');
        ylabel('Frequency (Normalized)');
    end
    legend( 'Patient 005', ...
            'Patient 006', ...
            'Patient 007', ...
            'Patient 008', ...
            'Patient 009', ...
            'Patient 010', ...
            'Patient 011', ...
            'Patient 012', ...
            'Patient 013', ...
            'Patient 014', ...
            'Patient 015', ...
            'Patient 016', ...
            'Patient 017', ...
            'Patient 018', ...
            'Patient 019', ...
            'Patient 019 fingerprick', ...
            'Patient 020', ...
            'Patient 020 fingerprick', ...
            'Patient 021', ...
            'CW', ...
            'XW' )
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
