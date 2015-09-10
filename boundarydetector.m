function [pk,lc] = boundarydetector(vpic)

for j = 1:1:size(vpic,1)
    hc(j) = sum(vpic(j,:))/(512*65535);
    vc(j) = sum(vpic(:,j))/(512*65535);
end
hc = hc - min(hc(:));
hc = hc/max(hc);
vc = vc - min(vc(:));
vc = vc/max(vc);

figure(9)
subplot(1,2,1)
plot(hc)
hold on
subplot(1,2,2)
plot(vc)
hold on

[pk.hc,lc.hc] = findpeaks(hc,'NPeaks',2,'MinPeakDistance',100,'SortStr','descend');
[pk.vc,lc.vc] = findpeaks(vc,'NPeaks',1,'MinPeakDistance',100,'SortStr','descend');

if abs(lc.hc(1) - lc.hc(2)) < 200
    figure(10)
    s(1) = subplot(1,2,1)
    plot(hc)
    hold on
    s(2) = subplot(1,2,2)
    plot(vc)
    hold on
    title(s(1),'Horizontal collection')
    title(s(2),'Vertical collection')
    lc.hc(1) = input('upbd : ');
    lc.hc(2) = input('lwbd : ');
end

end