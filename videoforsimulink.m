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

pm.fldname = [num2str(dt(3)),'_',num2str(dt(2)),'_',num2str(dt(1)),'_H',num2str(dt(4)),'M',num2str(dt(5))];

[pm.FileName,PathName] = uigetfile(filetype','Select the video file');
%
handles.pathname = [PathName,pm.FileName];
%handles.pathname = 'C:\Users\nwbi\Google Drive\NTU\Progress\drhou\HOUHANWEI TTSH\TTSH clinical testing (April 13 onwards)\04132015 WBCs-0005 n 0006\0005\20x-1million mL\Time Series-9.tif'
pm.FileName = pm.FileName(1:find(pm.FileName=='.')-1);
fld=dir(handles.pathname(1:max(findstr(handles.pathname,'\'))));
listoffilename = {fld.name};
k = 1;

close all
for i = 1:1:pm.nFrames
    clf
    pic(i).cdata = imread(handles.pathname,i);
end