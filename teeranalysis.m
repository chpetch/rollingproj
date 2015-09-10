
freqused = [1000 100000];
i = 1;
color = {'r' 'b'};

conf = {'EC' 'ED' 'FA' 'FB'}
for crun = 1:length(conf)
    for i = 1:2
        folder = ['C:\Users\nwbi\Google Drive\NTU\Progress\drhou\TEER\742015\Cells\D Normal _C','\'];
        data = xlsread([folder,conf{crun},'_f',num2str(freqused(i)),'_zmag']);
        phase = xlsread([folder,conf{crun},'_f',num2str(freqused(i)),'_zphase']);
        savename = [folder,conf{crun},'_f',num2str(freqused(i)),'_graph'];
        for k = 1:2
            figure (k)
            if k == 1
                suptitle(['Collected data across ',conf{crun},' using ',num2str(freqused(i)),' Hz'])
            end
            Z = data(:,2).*exp(j*phase(:,2));
            R = real(Z);
            X = imag(Z);
            hold on
            subplot(3,1,1)
            plot(data(:,1),data(:,2),color{i})
            ylabel('Impedance |Z| (Ohm)')
            xlabel('Time (s)')
            hold on
            subplot(3,1,2)
            plot(data(:,1),R,color{i})
            ylabel('Resistance (Ohm)')
            xlabel('Time (s)')
            hold on
            subplot(3,1,3)
            plot(data(:,1),X,color{i})
            ylabel('Reactance')
            xlabel('Time (s)')
            if k == 2
                hold on
                if i == 2
                    subplot(3,1,1)
                    ylim([0 8000])
                    subplot(3,1,2)
                    ylim([0 8000])
                    suptitle(['Comparison of collected data across ',conf{crun},' using ',num2str(freqused(1)),' Hz',' and ',num2str(freqused(2)),' Hz'])
                    saveas(figure(k),[savename,'_sum.jpg'])
                    clf
                    %legend(['Measured values from ',num2str(freqused(1)),'hz'],['Measured values from ',num2str(freqused(2)),'hz'],'Location','southoutside')
                end
            else
                saveas(figure(k),[savename,'.jpg'])
                clf
            end
        end
        
    end
end