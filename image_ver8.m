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
%maximum distance 
pm.d_threshold = 25;
%minimum number of frames that cell has to appear
pm.mft = 10;
pm.pkdist = 5;
pm.input_binsize = 9;
pm.convexarea = [25,40];
pm.umperpx = 0;
pm.step = 0.5;
pm.lng = 20;
pm.color = 'b';
%number of connected frames
pm.ccf = 5;
pm.smc = 30;
pm.removeslowspeed = 0;
fldop = 1;
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
        pm.FileName = handles.pathname(max(find(handles.pathname == '\'))+1:max(find(handles.pathname=='.'))-1);
        disp(['File ',num2str(frun) ,' Name: ',listoftif{frun}])
    end
    close all
    for i = 1:1:pm.nFrames
        clf
        pic(i).cdata = imread(handles.pathname,i);
        npic = double(pic(i).cdata);
        npic = npic-min(npic(:));
        npic = npic/max(npic(:));
        kkk(:,:,i) = npic;
        if i == 1
            for j = 1:1:size(pic(i).cdata,1)
                hc(j) = sum(pic(i).cdata(j,:))/(512*65535);
            end
            hc = hc - min(hc(:));
            hc = hc/max(hc);
            plot(hc)
            hold on
            %pm.pkheight = input('pkheight? ');
            if ischar(pm.pkheight) ~= 1
                [pk,lc] = findpeaks(hc,'MinPeakHeight',0.8,'MinPeakDistance',pm.pkdist);
                if length(lc) ~= 2
                    pm.upbd = input('pm.upbd: ');
                    pm.lwbd = input('pm.lwbd: ');
                else
                    pm.upbd = lc(1);
                    pm.lwbd = lc(2);
                end
                plot(lc,pk,'rx')
            end
        end
    end
    % zproject = zproject/(pm.nFrames*65535);
    if ischar(pm.pkheight) == 1
        pm.upbd = input('upbd : ');
        pm.lwbd = input('lwbd : ');
    end
    disp(['cutoff boundary : ' , num2str(pm.upbd) ,' - ', num2str(pm.lwbd)])
    if abs(pm.upbd - pm.lwbd)<200
        plot(hc)
        pm.lwbd = input('lwbd : ');
        pm.upbd = input('upbd : ');
    end
    npic = std(kkk,0,3);
    npic = npic/max(npic(:));
    pm.umperpx = 400/abs(pm.upbd-pm.lwbd);
    pause(0.5)
    clf
    %pause
    clearvars -except pic pm.upbd pm.lwbd pm.tbf pm.nFrames zproject npic pm handles fldop listoftif frun fldname
    
    %prepare save path
    savename = handles.pathname(1:end-4);
    fldname = pm.fldname;
    fldname = [savename(1:max(findstr(savename,'\'))),fldname];
    mkdir(fldname)
    
    for i = 1:1:pm.nFrames
        i
        procpic = pic(i).cdata;
        
        hgamma = ...
            vision.GammaCorrector(2.2,'Correction','De-gamma');
        y = step(hgamma, procpic);
        
        ori = procpic;
        procpic = y;
        
        procpic(1:pm.upbd+pm.tlr,:) = 0;
        procpic(pm.lwbd-pm.tlr:end,:) = 0;
        
        tophatpic1 = imtophat(procpic, strel('disk', 20));
        
        logpic = imagefilt(tophatpic1,'log',10,2.5);
        I = logpic;
        thresh = multithresh(I,2);
        seg_I = imquantize(I,thresh);
        RGB = label2rgb(seg_I);
        BW = (seg_I == 3);
        [BW_out,properties] = filterRegions(BW,10);
        center = cat(1, properties.Centroid);
        
        figure(1)
        imshow(ori)
        hold on
        plot(center(:,1), center(:,2), 'r*');
        title(['Detected cells in frame',num2str(i)])
        
        pause(0.5)
        
        features{i} = [fliplr(center)];
        
    end
    
    %match cell of frame t and frame t+1
    run = 1;
    d.index = [];
    d.amp = [];
    z.index = [];
    mf = features;
    
    sidis = matchingframes(mf);
    
    [t_cell,sframe] = calculatevelocity(mf,sidis,pm);
    
    [result,ntf] = constrainsforcells(t_cell,sframe,pm);
    
    if isempty(t_cell) ~= 1
        %image projection
        for i = 1:size(result.cpos_new,2)
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
            clf
            imshow(ncpic{j})
            hold on
            plot(result.cpos_new{j}(:,1),result.cpos_new{j}(:,2),'r', 'LineWidth', 2)
            legend(['Traveling Path of Cell ', num2str(j) ,' | Starting frame: ' , num2str(sframe(j)),' , Last frame: ' , num2str(sframe(j)+ntf(j)-1)])
            pause(0.5)
            saveas(figure(2),[fldname,'\',pm.FileName,'\tj_cell',num2str(j),'.jpg'])
            if mean(cvel{j}) < 0.5
                if pm.removeslowspeed == 1
                    reply = input('keep this cell? (1=yes,2=no) : ');
                    if reply == 2
                        htd = [htd,j];
                    end
                else
                    htd = [htd,j];
                end
            end
        end
        
        %delete cell
        result.cpos_new(htd) = [];
        clear cvel cvelxy
        
        %plot cell trajectory in graph
        figure(3)
        clf
        for j = 1:size(result.cpos_new,2)
            clear dataplot
            dataplot = bsxfun(@minus,result.cpos_new{j},result.cpos_new{j}(1,:));
            cdist{j} = dataplot;
            plot(dataplot(:,1),dataplot(:,2),'color',pm.color, 'LineWidth', 1)
            hold on
        end
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
        
        pm.vidvel{frun} = cvel;
        
        save([handles.pathname(1:end-4),'.mat'],'fpname','pm','result')
    else
        disp('no moving cell in this file')
        figure(6)
        clf
        imshow(npic)
        pause
    end
    
end

if fldop == 1
    acv = get_acv(pm.vidvel);
    xlswrite([fldname,'\acv.xls'], acv, 'absolute cell velocity', 'A1');
    %create histogram,trajectory that accumulates every cell in every videos
    output = histsum(acv,pm.step,pm.lng,fldname,pm.color,1,1)
    trajsum(pm,savename)
end