%determine circularity

clear
clc
%clf

filetype{1} = '*.tif';
filetype{2} = '*.avi';

%associated variables
pm.tbf = 0.5;
pm.nFrames = 61;
pm.d_threshold = 100;
%minimum number of frames that cell has to appear
pm.mft = 10;
pm.input_binsize = 9;
pm.step = 1;
pm.lng = 100;
pm.color = 'g';
%number of connected frames
pm.ccf = 5;
pm.smc = 30;
fldop = 0;
pm.bdd = 0;
pm.th_area = 900;
dt = clock;

pm.fldname = [num2str(dt(3)),'_',num2str(dt(2)),'_',num2str(dt(1)),'_H',num2str(dt(4)),'M',num2str(dt(5))];

[FileName,PathName] = uigetfile(filetype','Select the video file');

handles.pathname = [PathName,FileName];

fld=dir([handles.pathname(1:max(findstr(handles.pathname,'\'))),FileName(1:end-4),'\']);
listoffilename = {fld.name};
k = 1;
if fldop == 1;
    for i = 1:length(listoffilename)
        if findstr(listoffilename{i},'Time Series') == 1
            if findstr(listoffilename{i},'.tif') > 0
                listoftif{k} = listoffilename{i};
                
                k=k+1;
            end
        end
    end
    pm.nof = length(listoftif);
else
    pm.nof = 1;
end

%prepare save path
savename = handles.pathname(1:end-4);
fldname = pm.fldname;
fldname = [savename(1:max(findstr(savename,'\'))),FileName(1:end-4),'\'];
mkdir(fldname)
ba = 1;

for i = 1:60:pm.nFrames
    i
    pic(i).cdata = imread(handles.pathname,i);
    procpic = pic(i).cdata;
    
    %reverse color
    procpic = imcomplement(procpic);
    procpic = imsharpen(procpic);
    %store original iamge
    ori = procpic;
    %background subtraction
    background = imopen(procpic,strel('disk',7));
    I2 = procpic - background;
    %make circular objects stand out? do we need it?
    I2=imtophat(I2, strel('disk', 20));
    I3 = imadjust(I2);
    level = graythresh(I3);
    bw = im2bw(I3,level);
    %fill in the holes
    bw = imfill(bw,'holes');
    %delete obj smaller than the second agruement
    bw = bwareaopen(bw, 150);
    %delete border obj
    bw = imclearborder(bw, 4);
    bw = bwareaopen(bw, 150);
    %Hessian ridge detection (will be used as segmentation function)
    Hf=hessianridgefilter(bw,0);
    %Sobel edge detection (will be used as segmentation function)
    I_sobel = imagefilt(I3,'sobel',1,2);
    
    %determine foreground mask
    I = imcomplement(I3);
    %default value = 8
    se = strel('disk', 8);
    Io = imopen(I, se);
    Ie = imerode(I, se);
    Iobr = imreconstruct(Ie, I);
    Iobrd = imdilate(Iobr, se);
    Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
    Iobrcbr = imcomplement(Iobrcbr);
    level = graythresh(Iobrcbr);
    fgm = imregionalmin(Iobrcbr);
    se2 = strel('disk', 5);
    fgm3 = imerode(fgm, se2);
    fgm4 = bwareaopen(fgm3, 20);
    fgm4 = imclearborder(fgm4, 18);
    %fgm4 = imerode(bw,se2);
    %create foreground mask by merging 2 masks
    fgm4 = bwareaopen(fgm4, 20);
    %fgm4 = imdilate(fgm4,ones(2,2));
    fgm4 = fgm4&bw;
    fgm4 = bwareaopen(fgm4, 20);
    
    %determine background mask
    %D = ~bwdist(im2bw(Iobrcbr, graythresh(Iobrcbr)));
    D = bwdist(bw);
    DL = watershed(D);
    bgm = DL == 0;
    
    Hf = (Hf.^2);
    mHf = (Hf./max(Hf(:)))*65535;
    %mHf = imcomplement(mHf);
    
    gradmag2 = imimposemin(I_sobel, bgm|fgm4);
    %figure
    %imshow(gradmag2)
    L = watershed(gradmag2);
    
    I4 = I;
    I4(imdilate(L == 0, ones(3, 3)) | bgm | fgm4) = 65535;
    Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(1,2,1)
    imshow(ori)
    hold on
    subplot(1,2,2)
    imshow(I)
    himage = imshow(Lrgb);
    himage.AlphaData = 0.3;
    title('Segmented cells on original image')
    
    %extract infomation from bw image
    info = regionprops(bwconncomp(L,8),'Centroid','Eccentricity','MajorAxisLength','MinorAxisLength','Area','Perimeter');
    %delete background and large size elements in image
    area = cat(1,info.Area);
    delbg = find(area>pm.th_area);
    info(delbg) = [];    
    area = cat(1,info.Area);
    center = cat(1, info.Centroid);
    perimeter = cat(1, info.Perimeter);
    csi = (area*(4*pi))./((perimeter+pi).^2); 
        
    %show fancy plot blah blah
    for k = 1 : size(info,1)
        text(center(k,1),center(k,2), ...
            sprintf('%d: %1.3f',k, csi(k)), ...
            'Color','r');
    end
    
    %delete the cells that we don't want
    delcell = input('Which cell do you want to delete?');
    if (isempty(delcell)==0)
        clear csi
        info(delcell) = [];
        area = cat(1,info.Area);
        center = cat(1, info.Centroid);
        perimeter = cat(1, info.Perimeter);
        csi = (area*(4*pi))./((perimeter+pi).^2);
        %plot again
        clf
        subplot(1,2,1)
        imshow(ori)
        subplot(1,2,2)
        imshow(I)
        himage = imshow(Lrgb);
        himage.AlphaData = 0.3;
        title('Segmented cells on original image')
        
        for k = 1 : size(info,1)
        text(center(k,1),center(k,2), ...
            sprintf('%d: %1.3f',k, csi(k)), ...
            'Color','r');
        end
    end

    %aspect ratio    
    %old definition of circularity
    mn = cat(1,info.MinorAxisLength);
    mj = cat(1,info.MajorAxisLength);
    ar = mj./mn;
    %store necessary information
    store.features{ba} = [fliplr(center)];
    store.csi{ba} = csi;
    store.aspectratio{ba} = ar;
    ba = ba +1 ;
    pause
end

if isempty(store.csi{1})==0
    xlswrite([fldname,'\ccr.xls'], store.csi{1}, 'circularity of first frame');
end
if isempty(store.csi{2})==0
    xlswrite([fldname,'\ccr.xls'], store.csi{2}, 'circularity of last frame');
end
if isempty(store.aspectratio{1})==0
    xlswrite([fldname,'\ccr.xls'], store.aspectratio{1}, 'aspectratio of first frame');
end
if isempty(store.aspectratio{2})==0
    xlswrite([fldname,'\ccr.xls'], store.aspectratio{2}, 'aspectratio of last frame');
end
