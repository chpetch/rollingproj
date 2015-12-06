%determine circularity

clear
clc
%clf

filetype{1} = '*.tif';
filetype{2} = '*.avi';

%associated variables
pm.tbf = 0.5;
pm.nFrames = 61;
pm.d_threshold = 25;
pm.d_threshold_bwd = 0.5;
pm.d_ydis = 1;
%minimum number of frames that cell has to appear
pm.mft = 10;
pm.input_binsize = 9;
pm.step = 1;
pm.lng = 100;
pm.umperpx = 1;
pm.color = 'g';
%number of connected frames
pm.ccf = 5;
pm.smc = 30;
fldop = 0;
pm.bdd = 0;
pm.th_area = 900;
pm.removeslowspeed = 0;
dt = clock;

pm.fldname = [num2str(dt(3)),'_',num2str(dt(2)),'_',num2str(dt(1)),'_H',num2str(dt(4)),'M',num2str(dt(5))];

[FileName,PathName] = uigetfile(filetype','Select the video file');

handles.pathname = [PathName,FileName];
pm.FileName = FileName;
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

for i = 1:1:pm.nFrames
    i
    pic(i).cdata = imread(handles.pathname,i);
    procpic = pic(i).cdata;
    
    %reverse color
    procpic = imcomplement(procpic);
    procpic = imsharpen(procpic);
    %store original iamge
    ori = procpic;
    %background subtraction
    background = imopen(procpic,strel('disk',5));
    I2 = procpic - background;
    %make circular objects stand out? do we need it?
    I2=imtophat(I2, strel('disk', 20));
    
    if i == 1
       I3 = imadjust(I2);
       contrastrange = stretchlim(I2);
       level = graythresh(I3);
    else
       I3 = imadjust(I2,contrastrange); 
    end
    
    bw = im2bw(I3,level);
    %fill in the holes
    bw = imfill(bw,'holes');
    %delete obj smaller than the second agruement
    bw = bwareaopen(bw, 100);
    %delete border obj
    bw = imclearborder(bw, 4);
    bw = bwareaopen(bw, 100);

    %Hessian ridge detection (will be used as segmentation function)
    Hf=hessianridgefilter(bw,0);
    %Sobel edge detection (will be used as segmentation function)
    I_sobel = imagefilt(I3,'sobel',1,2);
    
    %determine foreground mask
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
    %D = bwdist(~bw2);
    D = bwdist(bw);
    DL = watershed(D,8);
    bgm = DL == 0;
    
    %mHf = imcomplement(mHf);
    
    gradmag2 = imimposemin(I_sobel, bgm|fgm4);
    %figure
    %imshow(gradmag2)
    L = watershed(gradmag2,8);
    
    I4 = I;
    I4(imdilate(L == 0, ones(3, 3)) | bgm | fgm) = 65535;
    Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
    %figure('units','normalized','outerposition',[0 0 1 1])
    subplot(1,3,1)
    imshow(bw)
    hold on
    subplot(1,3,2)
    imshow(Lrgb)
    title('Segmented cells on original image')
    subplot(1,3,3)
    imshow(bgm|fgm4)
    
    %extract infomation from bw image
    info = regionprops(bw,'Centroid','Eccentricity','MajorAxisLength','MinorAxisLength','Area','Perimeter');
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
    pause(0.25)
end
%match cell of frame t and frame t+1
run = 1;
d.index = [];
d.amp = [];
z.index = [];
mf = store.features;

sidis = matchingframes(mf);

[t_cell,sframe] = calculatevelocity(mf,sidis,pm);

[result,ntf,sframe] = constrainsforcells(t_cell,sframe,pm);

if isempty(t_cell) ~= 1
    %image projection
    for i = 1:size(result.cpos_new,2)
        clear kkk
        k = 1;
        for j = result.cse_new(2,i):1:result.cse_new(3,i)
            procpic = double(pic(j).cdata);
            procpic = procpic-min(procpic(:));
            procpic = procpic/max(procpic(:));
            kkk(:,:,k) = procpic;
            k = k+1;
        end
        procpic = std(kkk,0,3);
        ncpic{i} = procpic/max(procpic(:));
    end
    
    savename =[handles.pathname(1:end-4)];
    
    for i = 1:size(result.cpos_new,2)
        for j = 1 : size(result.cpos_new{i},1)-1
            cvel{i}(j) = (sqrt(sum((result.cpos_new{i}(j,:) - result.cpos_new{i}(j+1,:)).^2))*pm.umperpx)/pm.tbf;
            cvelxy{i}(j,:) = (result.cpos_new{i}(j,:) - result.cpos_new{i}(j+1,:)).*(pm.umperpx/pm.tbf);
        end
    end
    
    %show each cell trajectory on picture
    figure(2)
    %mkdir(savename)
    mkdir([fldname,'\',pm.FileName])
    htd = [];
    for j = 1:size(result.cse_new,2)
        figure(2)
        clf
        imshow(ncpic{j})
        hold on
        plot(result.cpos_new{j}(:,1),result.cpos_new{j}(:,2),'r', 'LineWidth', 2)
        legend(['Traveling Path of Cell ', num2str(j) ,' | Starting frame: ' , num2str(sframe(j)),' , Last frame: ' , num2str(sframe(j)+ntf(j)-1)])
        pause(0.5)
        saveas(figure(2),[fldname,'\',pm.FileName,'\tj_cell',num2str(j),'.jpg'])
        dataplot = bsxfun(@minus,result.cpos_new{j},result.cpos_new{j}(1,:));
        %auc(j)=linerotation(dataplot(:,1),dataplot(:,2),1)/abs(dataplot(end,1));
        %trace cell motion
        %cmpic = tracecellmotion(pic,result.cpos_new{j},sframe(j),1,[fldname,pm.FileName,'\tj_cell',num2str(j)]);
    end
    %auc=auc';
    
    %delete cell
    result.cpos_new(htd) = [];
    clear cvel cvelxy
    frun =1;
    %plot cell trajectory in graph
    figure(3)
    clf
    for j = 1:size(result.cpos_new,2)
        clear dataplot
        dataplot = bsxfun(@minus,result.cpos_new{j},result.cpos_new{j}(1,:));
        cdist{j} = dataplot;
        plot(dataplot(:,1),dataplot(:,2),'color','r', 'LineWidth', 1)
        hold on
        xlim([-300,0])
        ylim([-60,60])
    end
    cdist = plotcelltraj(result.cpos_new,[-300,0],[-60,60],'r');
    pm.distdata{frun} = cdist;
    
    figure(7)
    clf
    for j = 1:size(result.cpos_new,2)
        clear bcdisp
        bcdisp = result.cpos_new{j}(end,:) - result.cpos_new{j}(1,:);
        cdisp{j} = bcdisp;
        plot([0,bcdisp(1)],[0,bcdisp(2)],'color',pm.color, 'LineWidth', 1)
        hold on
        bcdisp = sqrt(sum(bcdisp.^2));
        bcdisp = (bcdisp/size(result.cpos_new{j},1))*(pm.umperpx/pm.tbf);
        cvel_disp{j} = bcdisp;
    end
    
    pm.dispdata{frun} = cdisp;
    
    fpname = handles.pathname;
    
    for i = 1:size(result.cpos_new,2)
        for j = 1 : size(result.cpos_new{i},1)-1
            cvel{i}(j) = (sqrt(sum((result.cpos_new{i}(j,:) - result.cpos_new{i}(j+1,:)).^2))*pm.umperpx)/pm.tbf;
            cvelxy{i}(j,:) = (result.cpos_new{i}(j,:) - result.cpos_new{i}(j+1,:)).*(pm.umperpx/pm.tbf);
        end
    end
    
    result.vidvel{frun} = cvel;
    
    save([handles.pathname(1:end-4),'.mat'],'fpname','pm','result')
else
    disp('no moving cell in this file')
    figure(6)
    clf
    imshow(npic)
    pause
end
