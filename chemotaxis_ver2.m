%%%changelog for ver6
% / centroid x , y and area are used as feature instead of using only centroid x,y
% / no more histogram of each file (only histogram from accumulated average velocity)
% / collected result are more well organized
%% Problems
% X lose detected cell may be introduced leading to creation of new
% trajectory

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
pm.tlr = 10;
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
fldop = 0;
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
        procpic(1:pm.upbd+pm.tlr,:) = 0;
        procpic(pm.lwbd-pm.tlr:end,:) = 0;
        
        tophatpic1 = imtophat(procpic, strel('disk', 5));
        conpic = imadjust(tophatpic1,[0.45 1]);
        
        npic = conpic/double(max(conpic(:)))*65536;
        se = strel(ones(3,3));
        gdpic = imdilate(npic, se) - imerode(npic, se);
        procpic=imfill(gdpic);
        procpic = medfilt2(procpic,[3,3]);
        imshow(im2bw(procpic))
        CC = bwconncomp(im2bw(procpic), 8);
        obj = regionprops(CC,'Area','centroid','EquivDiameter','ConvexArea','MajorAxisLength','MinorAxisLength');
        center = cat(1, obj.Centroid);
        aa=cat(1, obj.Area);
        eqd = cat(1, obj.EquivDiameter);
        convexarea = cat(1, obj.ConvexArea);
        mjmn = [cat(1, obj.MajorAxisLength) cat(1, obj.MinorAxisLength)];
        
        %constrains
        dcen.index = find(aa<mean(aa)/2);
        dcen.value = aa(find(aa<mean(aa)/2));
        dcen.center = center(find(aa<mean(aa)/2),:);
        delcell = find(aa<mean(aa)/2);
        %         aa(delcell) = [];
        %         center(delcell,:) = [];
        
        %         imshow(procpic)
%         hold on
%         plot(center(:,1), center(:,2), 'b*');
%         hold on
%         plot(dcen.center(:,1), dcen.center(:,2), 'r*')
        pause(0.5)
        
        features{i} = [fliplr(center) aa];
        
    end
    
    %match cell of frame t and frame t+1
    run = 1;
    d.index = [];
    d.amp = [];
    z.index = [];
    mf = features;
    
    for i = 1 : pm.nFrames - 1
        clear dis
        %calculate between i objects in frame 1 to every j objects in frame 2
        if (isempty(mf{1,i}) == 0) && (isempty(mf{1,i+1}) == 0)
            [idis,dis] = knnsearch(mf{1,i},mf{1,i+1});
            
            %find new point for the objects in frame 1 that shared the same object
            %in frame 2
            
            q = 1;
            
            for j = 1 : size(dis,1)
                fsame = find(idis == j);
                
                if size(fsame,1) > 1
                    %objects in frame 2 matched with the same object in frame 1
                    iss{i,q} = j;
                    ss{i,q} = fsame;
                    q = q + 1;
                end
            end
            
            %select the shortest path in case there are more than one objects in
            %frame 2 belong to object in frame 1
            
            % number of deleted value
            nsv = 0;
            
            if exist('iss') == 1
                a = i;
                if size(ss,1) < i
                    disp('no repetitive cells here')
                else
                    for b = 1:size(ss(a,:),2)
                        if isempty(ss{a,b}) == 0
                            [iig,ig] = min(dis(ss{a,b}));
                            for c = 1 : length(dis(ss{a,b}))
                                %disp(['a: ',num2str(a),' b: ',num2str(b),' c: ',num2str(c)])
                                if c ~= ig
                                    %get repetitive index
                                    %frame
                                    repinx(run,1) = i;
                                    %index
                                    repinx(run,2) = ss{a,b}(c);
                                    %value
                                    repinx(run,3) = dis(ss{a,b}(c));
                                    dis(ss{a,b}(c)) = 0;
                                    idis(ss{a,b}(c)) = 0;
                                    run = run + 1;
                                    nsv = nsv + 1;
                                end
                            end
                        end
                    end
                end
            end
            
            %find maximum distances travel in each frame
            [maxamp(i),maxpoint(i)] = max(dis);
            dindex = [find(dis>10)];
            dindex = [ i*ones(1,length(dindex))' dindex];
            
            %find index and amplitude of distances
            d.index = [d.index ; dindex];
            d.amp = [d.amp ; dis(dis>10)];
            
            %find zero distances travel in each frame
            clear dindex
            dindex = [find(dis==0)];
            dindex = [ i*ones(1,length(dindex))' dindex];
            z.index = [z.index ; dindex];
            
            %calculate mean
            d2f(i) = sum(dis)/(length(dis)-nsv);
            sidis{i} = idis;
            sdis{i} = dis;
            vel{i} = dis/pm.tbf;
            ncell(i) = (length(dis)-nsv);
            cnsv{i} = nsv;
        end
    end
    
    inow = 1;
    t_cell{1} = [];
    sframe(1) = 1;
    d_cell{1} = [];
    figure(1)
    
    %calculate path of cells and label cells
    for inow = 1 : pm.nFrames - 1
        
        inow
        for i = 1:size(sidis{inow},1)
            if sidis{inow}(i) ~= 0
                if abs((mf{1,inow+1}(i,2))-mf{1,inow}(sidis{inow}(i),2)) >= 0 || abs(mf{1,inow+1}(i,1)-mf{1,inow}(sidis{inow}(i),1)) >= 0
                    clear vcell
                    ctf1 = [(mf{1,inow}(sidis{inow}(i),2)),(mf{1,inow}(sidis{inow}(i),1))];
                    ctf2 = [(mf{1,inow+1}(i,2)),(mf{1,inow+1}(i,1))];
                    vcell = (sqrt(sum((ctf2 - ctf1).^2))*pm.umperpx)/pm.tbf;
                    if inow == 1
                        if length(t_cell) < i
                            if (ctf1(1) - ctf2(1)) >= 0
                                t_cell{i} = [ctf1;ctf2];
                                sframe(i) = inow;
                            end
                        else
                            if (ctf1(1) - ctf2(1)) >= 0
                                t_cell{i} = [t_cell{i};[ctf1;ctf2]];
                            end
                        end
                    else
                        if isempty(t_cell) ~= 1
                            if (isempty(t_cell{1}) == 1)
                                if (ctf1(1) - ctf2(1)) >= 0
                                    t_cell{1} = [ctf1;ctf2];
                                    sframe(1) = inow;
                                end
                            else
                                for j = 1:size(t_cell,2)
                                    if t_cell{j}(end,:) == ctf1
                                        if (t_cell{j}(end,1)- ctf2(1)) >= 0
                                            t_cell{j} = [t_cell{j};ctf2];
                                            break;
                                        end
                                    else
                                        if j == size(t_cell,2)
                                            if (ctf1(1) - ctf2(1)) >= 0
                                                t_cell{j+1} = [ctf1;ctf2];
                                                sframe(j+1) = inow;
                                            end
                                        end
                                    end
                                end
                            end
                        else
                            if (ctf1(1) - ctf2(1)) >= 0
                                t_cell{1} = [ctf1;ctf2];
                                sframe(1) = inow;
                            end
                        end
                    end
                    
                    %check whether empty or not?
                    for k = size(t_cell,2):-1:1
                        if isempty(t_cell{k}) == 1
                            t_cell(k) = [];
                            sframe(k) = [];
                        end
                    end
                    
                else
                    %non moving cells
                    disp('hey')
                    %plot((mf{1,inow}(sidis{inow}(i),2)), (mf{1,inow}(sidis{inow}(i),1)), 'b*', 'MarkerSize', 12)
                end
            else
                disp('Tada!!!')
            end
        end
    end
         for i = 1:size(t_cell,2)
            ntf(i) = size(t_cell{i},1);
        end       
        %delete cell that appear less than x (pm.mft) frames
        t_cell = t_cell(find(ntf>5));
        sframe = sframe(find(ntf>5));
        ntf = [];
        for i = 1:size(t_cell,2)
            ntf(i) = size(t_cell{i},1);
        end
    
    %cell number (1) starting frame (2) - end frame (3)
    result.cpos_original = t_cell;
    result.cse_original = [[1:size(t_cell,2)];sframe;sframe + ntf - ones(size(ntf,1),1)];
    
    %find 2 detected objs that seem to be the same (disconnected cells)
    cccs = [];
    for cc = 1 : size(result.cse_original,2)
        result.cse_original(3,cc);
        for ccc = cc+1 : size(result.cse_original,2)
            if (result.cse_original(2,ccc) - result.cse_original(3,cc)) > 0 ...
            && (result.cse_original(2,ccc) - result.cse_original(3,cc)) < pm.ccf
                if abs(result.cpos_original{1,cc}(end,2) - result.cpos_original{1,ccc}(1,2)) < 1
                    if abs(result.cpos_original{1,cc}(end,1) - result.cpos_original{1,ccc}(1,1)) < pm.d_threshold
                        cccs = [cccs;cc,ccc];
                        break
                    end
                end
            end
        end
    end
    
    %connect seem to be connected cells
    if isempty(cccs) == 0
        for i = size(cccs,1):-1:1
            t_cell{cccs(i,1)} = [t_cell{cccs(i,1)};t_cell{cccs(i,2)}];
        end
        
        t_cell(cccs(:,2)) = [];
        sframe(cccs(:,2)) = [];
        
        ntf = [];
        for i = 1:size(t_cell,2)
            ntf(i) = size(t_cell{i},1);
        end
    end
    result.cpos_original = t_cell;
    result.cse_original = [[1:size(t_cell,2)];sframe;sframe + ntf - ones(size(ntf,1),1)];
    
    if isempty(t_cell) ~= 1
        
        for i = 1:size(t_cell,2)
            for j = 1 : size(t_cell{i},1)-1
                acvel{i}(j) = (sqrt(sum((t_cell{i}(j,:) - t_cell{i}(j+1,:)).^2)));
                if acvel{i}(j) > pm.d_threshold
                    acvel{i}(j) = [];
                    t_cell{i}(j+1:end,:) = [];
                    disp('cell moving faster than 25 is detected')
                    break
                else
                    cvel{i}(j) = (sqrt(sum((t_cell{i}(j,:) - t_cell{i}(j+1,:)).^2))*pm.umperpx)/pm.tbf;
                    cvelxy{i}(j,:) = (t_cell{i}(j,:) - t_cell{i}(j+1,:)).*(pm.umperpx/pm.tbf);
                end
            end
            ntf(i) = size(t_cell{i},1);
        end
        
        %delete cell that appear less than x (pm.mft) frames
        t_cell = t_cell(find(ntf>pm.mft));
        sframe = sframe(find(ntf>pm.mft));
        ntf = [];
        for i = 1:size(t_cell,2)
            ntf(i) = size(t_cell{i},1);
        end
        
        %cell number (1) starting frame (2) - end frame (3)
        result.cse_new = [[1:size(t_cell,2)];sframe;sframe + ntf - ones(size(ntf,1),1)];
        result.cpos_new = t_cell;
        
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
        
        %show each cell trajectory on picture
        figure(2)
        %mkdir(savename)
        mkdir([fldname,'\',pm.FileName])
        for j = 1:size(result.cse_new,2)
            clf
            imshow(ncpic{j})
            hold on
            plot(result.cpos_new{j}(:,1),result.cpos_new{j}(:,2),'r', 'LineWidth', 2)
            legend(['Traveling Path of Cell ', num2str(j) ,' | Starting frame: ' , num2str(sframe(j)),' , Last frame: ' , num2str(sframe(j)+ntf(j)-1)])
            pause(0.5)
            saveas(figure(2),[fldname,'\',pm.FileName,'\tj_cell',num2str(j),'.jpg'])
        end
        
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
    get_acv
    xlswrite([fldname,'\acv.xls'], acv, 'absolute cell velocity', 'A1');
    %create histogram,trajectory that accumulates every cell in every videos
    histsum(pm.vidvel,pm.step,pm.lng,fldname,pm.color,1,1)
    trajsum(pm,savename)
end