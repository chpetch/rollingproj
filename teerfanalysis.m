clc
clear
name = { ['C:\Users\nwbi\Google Drive\NTU\Progress\drhou\TEER\14042015\'];
    ['C:\Users\nwbi\Google Drive\NTU\Progress\drhou\TEER\14042015\treated\'];
    ['C:\Users\nwbi\Google Drive\NTU\Progress\drhou\TEER\14042015\treated\afterawhile\']}
chip = {'D'};
conf = {'FB' 'FA'};
color = {'r','b','m'}
i = 1;
h = 2;
for g = 1 : 3
    for i = 1 : 1
        data = xlsread([name{g},'Chip',chip{i},'_',conf{h},'_100_100k_mag']);
        phase = xlsread([name{g},'Chip',chip{i},'_',conf{h},'_100_100k_phase']);
        
        Z = data(:,2).*exp(j*phase(:,2));
        R = real(Z);
        X = imag(Z);
        pdata(:,2) = log10(data(:,2));
        
        figure(1)
        semilogx(data(:,1),pdata(:,2),color{g})
        hold on
        [v.min,v.idx]=min(pdata(:,2));
        %semilogx(data(v.idx,1),pdata(v.idx,2),'ro')
        hold on
        figure(2)
        semilogx(data(:,1),R,color{g})
        hold on
        %semilogx(data(v.idx,1),R(v.idx),'ro')
        figure(3)
        semilogx(data(:,1),X,color{g})
        hold on
        %semilogx(data(v.idx,1),X(v.idx),'ro')
        cZ{g} = Z;
    end
end

data_c = xlsread([name{1},'ChipF_',conf{h},'_100_100k_mag']);
phase_c = xlsread([name{1},'ChipF_',conf{h},'_100_100k_phase']);
Z = data_c(:,2).*exp(j*phase_c(:,2));
R = real(Z);
X = imag(Z);
data_c(:,2) = log10(data_c(:,2));
control_Z = Z;

figure(1)
semilogx(data_c(:,1),data_c(:,2),'k')
ylabel('Impedance |Z| (Ohm) in Log scale')
xlabel('Frequency (Hz)')
xlim([100 100000])
t= title(['Impedance comparison between chip ',chip{i},' and control (chip F)'])
set(t, 'FontSize', 20)
legend('1st attempt','2nd attempt (after treatment)','3rd attempt','Control','location','eastoutside')
hold on
figure(2)
semilogx(data_c(:,1),R,'k')
ylabel('Resistance (Ohm)')
xlabel('Frequency (Hz)')
xlim([100 100000])
t = title(['Resistance comparison between chip ',chip{i},' and control (chip F)'])
set(t, 'FontSize', 20)
legend('1st attempt','2nd attempt (after treatment)','3rd attempt','Control','location','eastoutside')
hold on
figure(3)
semilogx(data_c(:,1),X,'k')
ylabel('Reactance')
xlabel('Frequency (Hz)')
xlim([100 100000])
t = title(['Reactance comparison between chip ',chip{i},' and control (chip F)'])
set(t, 'FontSize', 20)
refline(0,1)
legend('1st attempt','2nd attempt (after treatment)','3rd attempt','Control','location','eastoutside')
hold on

%resistance and reactance

im = imag(control_Z);
re = real(control_Z);
CResistance = (im.^2 + re.^2)./re;
CReactance = (im.^2 + re.^2)./im;

for h = 1:3
    im = imag(cZ{h});
    re = real(cZ{h});
    Resistance = (im.^2 + re.^2)./re;
    Reactance = (im.^2 + re.^2)./im;
    figure(5)
    semilogx(data_c(:,1),(Resistance-CResistance),color{h})
    diff.resis{h} = Resistance-CResistance;
    hold on
    figure(6)
    semilogx(data_c(:,1),Reactance-CReactance,color{h})
    diff.react{h} = Reactance-CReactance;
    hold on
end

figure(5)
ylabel('Resistance')
xlabel('Frequency (Hz)')
t = title(['Resistance Difference between chip ',chip{i},' and control (chip F)'])
set(t, 'FontSize', 20)
refline(0,1)
legend('1st attempt','2nd attempt (after treatment)','3rd attempt','location','eastoutside')
figure(6)
ylabel('Reactance')
ylim([-0.5*10^5 0.5*10^5]);
xlabel('Frequency (Hz)')
t = title(['Reactance Difference between chip ',chip{i},' and control (chip F)'])
set(t, 'FontSize', 20)
refline(0,1)
legend('1st attempt','2nd attempt (after treatment)','3rd attempt','location','eastoutside')

% figure(5)
% semilogx(data_c(:,1),log10(CResistance),'k')
% hold on
% figure(6)
% semilogx(data_c(:,1),CReactance,'k')
% hold on