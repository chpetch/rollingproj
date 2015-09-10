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
pm.mft = 10;
pm.pkdist = 5;
pm.input_binsize = 9;
pm.convexarea = [25,40];
pm.umperpx = 0;
pm.step = 1;
pm.lng = 100;
pm.color = 'g';
fldop = 0;
pm.bdd = 0;

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

for frun = 1:pm.nof
    if fldop == 1
        handles.pathname = [handles.pathname(1:max(findstr(handles.pathname,'\'))),listoftif{frun}];
        disp(['File ',num2str(frun) ,' Name: ',listoftif{frun}])
    end
    close all
    for i = 1:1:pm.nFrames
        clf
        pic(i).cdata = imread(handles.pathname,i);
    end
    
    clf
    clearvars -except pic pm.upbd pm.lwbd pm.tbf pm.nFrames zproject npic pm handles fldop listoftif frun
    
    
    for i = 1:1:pm.nFrames
        i
        tic
        procpic = pic(i).cdata;
        procpic = imcomplement(procpic);
        I2 = imtophat(procpic, strel('disk', 5));
        
        I2 = I2*(65535/double(max(I2(:))));
        I3 = I2;
        I3 = medfilt2(I2,[3,3]);
        I3 = imfill(I3,'holes');
        bI3 = im2bw(I3,graythresh(I3));
        bI3 = bwareaopen(bI3,100,8);
        tobj = regionprops(bI3,'centroid','EquivDiameter','ConvexArea','MajorAxisLength','MinorAxisLength');
        mj = cat(1, tobj.MajorAxisLength);
        mn = cat(1, tobj.MinorAxisLength);
        cva = cat(1, tobj.ConvexArea);
        ct = cat(1, tobj.Centroid);
        if find(mn < 8) > 0
            delin = find(mn<8);
            mn(delin) = [];
            mj(delin) = [];
            cva(delin) = [];
            ct(delin,:) = [];
        end
        obj.convexarea{i} = cva;
        obj.center{i} = ct;
        obj.mj{i} = mj;
        obj.mn{i} = mn;
        obj.mjmn{i} = mj./mn;
        rconvex = max(cva) - min(cva);
        mp = procpic;
        mp(~bI3) = 0;
        histpic = imhist(mp,65536);
        histpic(1) = 0;
        dispim2 = imadjust(procpic,[min(find(histpic~=0)) max(find(histpic~=0))]/65536,[]);
        
        sc = 25;
        msc = 15;
        for j = 1:size(obj.center{i},1)
            if obj.center{i}(j,1)+sc > size(procpic,1)
                croparea = [obj.center{i}(j,2)-sc obj.center{i}(j,2)+sc ; obj.center{i}(j,1)-sc 512];      
            else
                croparea = [obj.center{i}(j,2)-sc obj.center{i}(j,2)+sc ; obj.center{i}(j,1)-sc obj.center{i}(j,1)+sc ];
            end
            detectedcells{j} = procpic(obj.center{i}(j,2)-sc:obj.center{i}(j,2)+sc,obj.center{i}(j,1)-sc:obj.center{i}(j,1)+sc);
            
            mask = zeros(size(detectedcells{j}));
            mask(msc:end-msc,msc:end-msc) = 1;
            bw{i,j} = activecontour(detectedcells{j},mask,5000,'Chan-Vese');
            
            figure(3)
            subplot(1,size(obj.center{i},1),j)
            bw{i,j} = imfill(bw{i,j},'holes');
            bw{i,j} = bwareaopen(bw{i,j},20,8);
            mp = detectedcells{j};
            mp(~bw{i,j}) = 0;
            histpic = imhist(mp,65536);
            histpic(1) = 0;
            dispim = imadjust(mp,[min(find(histpic~=0)) max(find(histpic~=0))]/65536,[]);
            bbwbd = bwboundaries(bw{i,j});
            bwbd{i,j} = fliplr(bbwbd{1});
            imshow(dispim)
        end
        dispim = imadjust(mp,[min(find(histpic~=0)) max(find(histpic~=0))]/65536,[]);
        figure(1)
        imshow(procpic)
        hold on
        for j = 1:size(obj.center{i},1)
            plot(obj.center{i}(j,1)-sc+bwbd{i,j}(:,1),obj.center{i}(j,2)-sc+bwbd{i,j}(:,2),'r')
        end
        plot(ct(:,1), ct(:,2), 'b*');
        pause(0.5)
    end
    
end
% I3 = procpic;
% %I3 = medfilt2(procpic,[3,3]);
% I2 = imtophat(I3, strel('disk', 5));
% I2 = I2*(65535/double(max(I2(:))));
% I2 = imfill(I2,'holes');
% subplot(1,2,2)
% bI2 = im2bw(I2,graythresh(I3));
% obj2 = regionprops(bI2,'centroid','EquivDiameter','ConvexArea','MajorAxisLength','MinorAxisLength');
% imshow(bI2)
%
% B2 = imbothat(I3, strel('disk', 5));
% B2 = B2*(65535/double(max(B2(:))));
% imshow(B2)
% hp = imhist(procpic);
% [~,idx] = max(hp);
%
% imshow(procpic)

% i = 1;
% procpic = pic(i).cdata;
% procpic = imcomplement(procpic);
%
% subplot(2,2,1)
% I1 = imtophat(procpic, strel('disk', 5));
% I1 = I1*(65535/double(max(I1(:))));
% I1 = medfilt2(I1);
% I1 = imtophat(I1, strel('disk', 5));
% imshow(I1)
%
% subplot(2,2,2)
% I2 = imtophat(procpic, strel('disk', 10));
% I2 = I2*(65535/double(max(I2(:))));
% I2 = medfilt2(I2);
% I2 = imtophat(I2, strel('disk', 10));
% imshow(I2)
%
% subplot(2,2,3)
% I3 = imtophat(procpic, strel('disk', 15));
% I3 = I3*(65535/double(max(I3(:))));
% I3 = medfilt2(I3);
% I3 = imtophat(I3, strel('disk', 15));
% imshow(I3)
%
% subplot(2,2,4)
% I4 = imtophat(procpic, strel('disk', 20));
% I4 = I4*(65535/double(max(I4(:))));
% I4 = medfilt2(I4);
% I4 = imtophat(I4, strel('disk', 20));
% imshow(I4)