clear
clc
clf

filetype{1} = '*.tif';
filetype{2} = '*.avi';

%associated variables
pm.tbf = 0.5;
pm.contrast = [0.6 1];
pm.nFrames = 51;
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
pm.bdd = 1;

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
        imshow(pic(i).cdata)
        npic = double(pic(i).cdata);
        npic = npic-min(npic(:));
        npic = npic/max(npic(:));
        kkk(:,:,i) = npic;
    end
end

clf
for i = 1:1:pm.nFrames
    i
    procpic = pic(i).cdata;
    if pm.bdd == 1
        procpic(1:pm.upbd+pm.tlr,:) = 0;
        procpic(pm.lwbd-pm.tlr:end,:) = 0;
    end
    
    %watershed segmentation
    I2 = imtophat(procpic, strel('disk', 10));
    
    procpic = procpic*(65535/double(max(procpic(:))));
    procpic = edge(procpic);
    procpic = imdilate(procpic,strel('disk', 2));
    procpic = imfill(procpic,'holes');
    imshow(procpic)

    [center, radii, metric] = imfindcircles(procpic,[10,50]);
    hold on
    viscircles(center, radii,'EdgeColor','b');
    
    pause(0.5)
    centers{i} = fliplr(center);
    
end
