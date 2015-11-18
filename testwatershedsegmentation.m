%determine circularity

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
pm.tlr = 5;
pm.d_threshold = 100;
%minimum number of frames that cell has to appear
pm.mft = 10;
pm.pkdist = 5;
pm.input_binsize = 9;
pm.convexarea = [25,40];
pm.umperpx = 0;
pm.step = 1;
pm.lng = 100;
pm.color = 'g';
%number of connected frames
pm.ccf = 5;
pm.smc = 30;
fldop = 0;
pm.bdd = 0;
pm.th_area = 700;
dt = clock;

pm.fldname = [num2str(dt(3)),'_',num2str(dt(2)),'_',num2str(dt(1)),'_H',num2str(dt(4)),'M',num2str(dt(5))];

[FileName,PathName] = uigetfile(filetype','Select the video file');

handles.pathname = [PathName,FileName];

fld=dir(handles.pathname(1:max(findstr(handles.pathname,'\'))));
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
fldname = [savename(1:max(findstr(savename,'\'))),fldname];
mkdir(fldname)
ba = 1;

for i = 1:60:pm.nFrames
    i
    pic(i).cdata = imread(handles.pathname,i);
    procpic = pic(i).cdata;
    %     hgamma = ...
    %             vision.GammaCorrector(2.2,'Correction','De-gamma');
    %     y = step(hgamma, procpic);
    %     procpic = y;
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
    fgm4 = fgm4&bw;
    fgm4 = bwareaopen(fgm3, 20);
    
    %determine background mask
    D = ~bwdist(im2bw(Iobrcbr, graythresh(Iobrcbr)));
    %D = ~bwdist(bw);

    DL = watershed(D);
    bgm = DL == 0;
    
    %mHf = uint16((Hf.^2)*65535);
    gradmag2 = imimposemin(Hf, fgm4|bgm);
    L = watershed(gradmag2);
       
    I4 = I;
    I4(imdilate(L == 0, ones(3, 3)) | bgm | fgm4) = 255;
    Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
    figure
    imshow(I)
    hold on
    himage = imshow(Lrgb);
    himage.AlphaData = 0.3;
    title('Lrgb superimposed transparently on original image')
    
    %extract infomation from bw image
    info = regionprops(L,'Centroid','Eccentricity','MajorAxisLength','MinorAxisLength','Area');
    %delete background and large size elements in image
    area = cat(1,info.Area);
    info(find(area>pm.th_area)) = [];
    mn = cat(1,info.MinorAxisLength);
    mj = cat(1,info.MajorAxisLength);
    center = cat(1, info.Centroid);
    circularity = mj./mn;
    %show fancy plot blah blah
    for k = 1 : size(info,1)-1
        text(center(k,1),center(k,2), ...
            sprintf('%1.3f', circularity(k)), ...
            'Color','r');
    end
    
    %store necessary information
    store.features{ba} = [fliplr(center)];
    store.circularity{ba} = mj./mn;
    ba = ba +1 ;
    pause
end

%xlswrite([fldname,'\ccr.xls'], store.circularity{1}, 'circularity of first frame');
%xlswrite([fldname,'\ccr.xls'], store.circularity{2}, 'circularity of last frame');

%match cell of frame t and frame t+1
% run = 1;
% d.index = [];
% d.amp = [];
% z.index = [];
% mf = store.features;
%
% sidis = matchingframes(mf);
%
% [t_cell,sframe] = calculatevelocity(mf,sidis,pm);
%
% for i = 1:size(t_cell,2)
%     ntf(i) = size(t_cell{i},1);
% end
%
% %cell number (1) starting frame (2) - end frame (3)
% result.cpos_original = t_cell;
% result.cse_original = [[1:size(t_cell,2)];sframe;sframe + ntf - ones(size(ntf,1),1)];
%
% %find 2 detected objs that seem to be the same (disconnected cells)
% cccs = [];
% for cc = 1 : size(result.cse_original,2)
%     result.cse_original(3,cc);
%     for ccc = cc+1 : size(result.cse_original,2)
%         if (result.cse_original(2,ccc) - result.cse_original(3,cc)) > 0 ...
%                 && (result.cse_original(2,ccc) - result.cse_original(3,cc)) < pm.ccf
%             if abs(result.cpos_original{1,cc}(end,2) - result.cpos_original{1,ccc}(1,2)) < 1
%                 if abs(result.cpos_original{1,cc}(end,1) - result.cpos_original{1,ccc}(1,1)) < pm.d_threshold
%                     cccs = [cccs;cc,ccc];
%                     break
%                 end
%             end
%         end
%     end
% end
%
% %connect seem to be connected cells
% if isempty(cccs) == 0
%     for i = size(cccs,1):-1:1
%         t_cell{cccs(i,1)} = [t_cell{cccs(i,1)};t_cell{cccs(i,2)}];
%     end
%
%     t_cell(cccs(:,2)) = [];
%     sframe(cccs(:,2)) = [];
%
%     ntf = [];
%     for i = 1:size(t_cell,2)
%         ntf(i) = size(t_cell{i},1);
%     end
% end
% result.cpos_original = t_cell;
% result.cse_original = [[1:size(t_cell,2)];sframe;sframe + ntf - ones(size(ntf,1),1)];
%
% if isempty(t_cell) ~= 1
%
%     %delete cells that move faster than pm.d_threshold
%     for i = 1:size(t_cell,2)
%         for j = 1 : size(t_cell{i},1)-1
%             acvel{i}(j) = (sqrt(sum((t_cell{i}(j,:) - t_cell{i}(j+1,:)).^2)));
%             if acvel{i}(j) > pm.d_threshold
%                 %acvel{i}(j) = [];
%                 %t_cell{i}(j+1:end,:) = [];
%                 disp('cell moving faster than threshold is detected')
%                 break
%             end
%         end
%         ntf(i) = size(t_cell{i},1);
%     end
%
%     %delete cell that appear less than x (pm.mft) frames
%     t_cell = t_cell(find(ntf>pm.mft));
%     sframe = sframe(find(ntf>pm.mft));
%     ntf = [];
%     for i = 1:size(t_cell,2)
%         ntf(i) = size(t_cell{i},1);
%     end
%
%     %cell number (1) starting frame (2) - end frame (3)
%     result.cse_new = [[1:size(t_cell,2)];sframe;sframe + ntf - ones(size(ntf,1),1)];
%     result.cpos_new = t_cell;
%
%     %image projection
%         for i = 1:size(result.cpos_new,2)
%             k = 1;
%             for j = result.cse_new(2,i):1:result.cse_new(3,i)
%                 procpic = double(pic(j).cdata);
%                 procpic = procpic-min(procpic(:));
%                 procpic = procpic/max(procpic(:));
%                 kkk(:,:,k) = procpic;
%                 k = k+1;
%             end
%             procpic = std(kkk,0,3);
%             ncpic{i} = procpic/max(procpic(:));
%         end
% end