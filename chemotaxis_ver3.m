%chemotaxis version 3
clear
clc
clf

filetype{1} = '*.tif';
filetype{2} = '*.avi';

%associated variables
pm.tbf = 0.5;
pm.contrast = [0.6 1];
pm.nFrames = 61;
pm.pkheight = 0.15;
pm.tlr = 15;
pm.d_threshold = 25;
pm.mft = 10;
pm.pkdist = 5;
pm.input_binsize = 9;
pm.convexarea = [25,40];
pm.umperpx = 0;
pm.step = 0.5;
pm.lng = 20;
pm.color = 'b';
pm.ccf = 5;
pm.smc = 30;
dt = clock;

% pm.fldname = [num2str(dt(3)),'_',num2str(dt(2)),'_',num2str(dt(1)),'_H',num2str(dt(4)),'M',num2str(dt(5))];
%
% [pm.FileName,PathName] = uigetfile(filetype','Select the video file');
% %
% handles.pathname = [PathName,pm.FileName];
% %handles.pathname = 'C:\Users\nwbi\Google Drive\NTU\Progress\drhou\HOUHANWEI TTSH\TTSH clinical testing (April 13 onwards)\04132015 WBCs-0005 n 0006\0005\20x-1million mL\Time Series-9.tif'
% pm.FileName = pm.FileName(1:find(pm.FileName=='.')-1);
% fld=dir(handles.pathname(1:max(findstr(handles.pathname,'\'))));
% listoffilename = {fld.name};
% k = 1;

close all
clf
[pm.FileName,PathName] = uigetfile(filetype','Select the video file');
handles.pathname = [PathName,pm.FileName];
bw = input('black(0) or white(1) : ');
for i = 1:1:2
       
    vpic = imread(handles.pathname,i);
    
    if bw == 1
        vpic = imcomplement(vpic);
    end
    
    [pk,lc] = boundarydetector(vpic);
    pk.hc = sort(pk.hc);
    lc.hc = sort(lc.hc);    
    croppic = vpic(lc.hc(1)+5:lc.hc(2)-5,lc.vc:lc.vc+(end/2));
    scl = 400/(lc.hc(2)-lc.hc(1));
    croppic = imsharpen(croppic);
    
    hgamma = ...
        vision.GammaCorrector(2.2,'Correction','De-gamma');
    y = step(hgamma, croppic);
    
    ori = croppic;
    croppic = y;
    
    thpic = imtophat(croppic, strel('disk', 5));
    logpic = imagefilt(thpic,'log',9,3);
    
    I = logpic;
    thresh = multithresh(I,2);
    seg_I = imquantize(I,thresh);
    RGB = label2rgb(seg_I);
    BW = (seg_I == 3);
    [BW_out,properties] = filterRegions(BW,10);
    center = cat(1, properties.Centroid);
    
    figure(1)
    subplot(2,1,i)
    imshow(ori)
    hold on
    plot(center(:,1), center(:,2), 'r*');
    title(['Detected cells in ',pm.FileName])
    
    figure(2)
    subplot(2,1,i)
    color = 'r';
    [histFreq, histXout] = hist(center(:,1),[0:2:size(seg_I,2)]);
    output.value{i} = histFreq';
    output.freq{i} = histXout';
    output.raw{i} = center(:,1);
    hff = (histFreq/sum(histFreq)*100);
    bar(histXout, histFreq/sum(histFreq)*100,1,color)
    hold on
    PD = fitdist(center(:,1), 'normal');
    plot(histXout, pdf(PD, histXout)*100,['b' '-'],'LineWidth', 3);
    title(['Cell distribution in ',pm.FileName])
    
end