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

for i = 1:1:pm.nFrames
    i
    pic(i).cdata = imread(handles.pathname,i);
    procpic = pic(i).cdata;
    procpic = imcomplement(procpic);
    
    ori = procpic;
    
    background = imopen(procpic,strel('disk',7));
    imshow(background)
    
    I2 = procpic - background;
    
    I3 = imadjust(I2);
        
    level = graythresh(I3);
    bw = im2bw(I3,level);
    bw = imfill(bw,'holes');
    bw = bwareaopen(bw, 150);
    imshow(bw)
    
    pause
    
    
end
