%used for patient data 13 14

clear
clc
clf

filetype{1} = '*.tif';
filetype{2} = '*.avi';

dt = datetime;
fldname = [num2str(dt.Day),'_',num2str(dt.Month),'_',num2str(dt.Year),'_H',num2str(dt.Hour),'M',num2str(dt.Minute)];
fldname = ['G:\drhou\collected results\',fldname];
mkdir(fldname)

%associated variables
pm.tbf = 0.5;
pm.contrast = [0.6 1];
pm.nFrames = 61;
pm.pkheight = 0.15;
pm.tlr = 5;
pm.d_threshold = 25;
pm.mft = 10;
pm.pkdist = 5;
pm.input_binsize = 9;
pm.convexarea = [25,40];
pm.umperpx = 0;
pm.step = 0.5;
pm.lng = 20;
pm.color = 'r';
pm.fldname = fldname;
fldop = 1;

[FileName,PathName] = uigetfile(filetype','Select the video file');
%
handles.pathname = [PathName,FileName];
%handles.pathname = 'C:\Users\nwbi\Google Drive\NTU\Progress\drhou\HOUHANWEI TTSH\TTSH clinical testing (April 13 onwards)\04132015 WBCs-0005 n 0006\0005\20x-1million mL\Time Series-9.tif'
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
            pm.pkheight = input('pkheight? ');
            if ischar(pm.pkheight) ~= 1
                [pk,lc] = findpeaks(hc,'MinPeakHeight',pm.pkheight,'MinPeakDistance',pm.pkdist);
                
                for j = 1:1:length(pk)-1
                    diffpeak(j) = abs(pk(j) - pk(j+1));
                end
                [~,idx] = min(diffpeak);
                pm.upbd = lc(idx);
                pm.lwbd = lc(idx+1);
                
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
        pm.lwbd = input('lwbd : ');
        pm.upbd = input('upbd : ');
    end
    npic = std(kkk,0,3);
    npic = npic/max(npic(:));
    pm.umperpx = 400/abs(pm.upbd-pm.lwbd);
    pause(0.5)
    clf
    %pause
    clearvars -except pic pm.upbd pm.lwbd pm.tbf pm.nFrames zproject npic pm handles fldop listoftif frun
    
    for i = 1:1:pm.nFrames
        i
        procpic = pic(i).cdata;
        procpic(1:pm.upbd+pm.tlr,:) = 0;
        procpic(pm.lwbd-pm.tlr:end,:) = 0;
        
        tophatpic1 = imtophat(procpic, strel('disk', 5));

        conpic = imadjust(tophatpic1,[0.45 1]);
        npic = conpic/double(max(conpic(:)))*65536;
        se = strel(ones(3,3));
        %         gdpic = imdilate(npic, se) - imerode(npic, se);
        %         procpic=imfill(gdpic);
        procpic = medfilt2(npic,[3,3]);
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
        aa(find(aa<mean(aa)/2)) = [];
        center(find(aa<mean(aa)/2),:) = [];
             
        imshow(procpic)
        hold on
        plot(center(:,1), center(:,2), 'b*');
        hold on
        plot(dcen.center(:,1), dcen.center(:,2), 'r*')
        pause(0.5)
        
        centers{i} = fliplr(center);
    end
    
    %match cell of frame t and frame t+1
    run = 1;
    d.index = [];
    d.amp = [];
    z.index = [];
    mf = centers;
    
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
                            if (vcell < pm.d_threshold) && ((ctf1(1) - ctf2(1)) >= 0)
                                t_cell{i} = [ctf1;ctf2];
                                sframe(i) = inow;
                            end
                        else
                            if (vcell < pm.d_threshold) && ((ctf1(1) - ctf2(1)) >= 0)
                                t_cell{i} = [t_cell{i};[ctf1;ctf2]];
                            end
                        end
                    else
                        if isempty(t_cell) ~= 1
                            if (isempty(t_cell{1}) == 1)
                                if (vcell < pm.d_threshold) && ((ctf1(1) - ctf2(1)) >= 0)
                                    t_cell{1} = [ctf1;ctf2];
                                    sframe(1) = inow;
                                end
                            else
                                for j = 1:size(t_cell,2)
                                    if t_cell{j}(end,:) == ctf1
                                        if (sqrt(sum((t_cell{j}(end,:) - ctf2).^2)) < pm.d_threshold) && ((t_cell{j}(end,1)- ctf2(1)) >= 0)
                                            t_cell{j} = [t_cell{j};ctf2];
                                            break;
                                        end
                                    else
                                        if j == size(t_cell,2)
                                            if (vcell < pm.d_threshold) && ((ctf1(1) - ctf2(1)) >= 0)
                                                t_cell{j+1} = [ctf1;ctf2];
                                                sframe(j+1) = inow;
                                            end
                                        end
                                    end
                                end
                            end
                        else
                            if (vcell < pm.d_threshold) && ((ctf1(1) - ctf2(1)) >= 0)
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
    
    %cell number (1) starting frame (2) - end frame (3)
    result.cpos_original = t_cell;
    result.cse_original = [[1:size(t_cell,2)];sframe;sframe + ntf - ones(size(ntf,1),1)];
    
    %delete cell that appear less than x (pm.mft) frames
    t_cell = t_cell(find(ntf>pm.mft));
    sframe = sframe(find(ntf>pm.mft));
    clear ntf
    
    if isempty(t_cell) ~= 1
        for i = 1:size(t_cell,2)
            ntf(i) = size(t_cell{i},1);
            for j = 1 : size(t_cell{i},1)-1
                cvel{i}(j) = (sqrt(sum((t_cell{i}(j,:) - t_cell{i}(j+1,:)).^2))*pm.umperpx)/pm.tbf;
                cvelxy{i}(j,:) = (t_cell{i}(j,:) - t_cell{i}(j+1,:)).*(pm.umperpx/pm.tbf);
            end
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
        clf
        for j = 1:size(result.cse_new,2)
            imshow(ncpic{j})
            hold on
            plot(result.cpos_new{j}(:,1),result.cpos_new{j}(:,2),'r', 'LineWidth', 2)
            legend(['Traveling Path of Cell ', num2str(j) ,' | Starting frame: ' , num2str(sframe(j)),' , Last frame: ' , num2str(sframe(j)+ntf(j)-1)])
            pause(0.5)
            saveas(figure(2),[savename,'_tj_cell',num2str(j),'.jpg'])
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
        
        saveas(figure(3),[savename,'_tj_graph.jpg'])
        
        figure(4)
        clf
        [histFreq, histXout] = hist(cell2mat(cvel),[0:pm.step:pm.lng]);
        
        %delete cells without pair
        bar(histXout, histFreq/sum(histFreq)*100,1,pm.color);
        xlabel('Rolling cel velocity (um/s)');
        ylabel('Frequency (percent)');
        %saving part%
        sv = pm.step/2;
        for i = 1 : length(histXout)
            if i == 1
                bandwidth{i} = ['0 - ',num2str(sv)];
            elseif i == length(histXout)
                bandwidth{i} = [num2str(sv),' - ','more than ',num2str(pm.d_threshold)];
            else
                bandwidth{i} = [num2str(sv),' - ',num2str(sv+pm.step)];
                sv = sv + pm.step;
            end
        end
        tb = {'Bin Center (X)','Frequency (Y)', 'Percentage'};
        tbb = [histXout',histFreq',(histFreq/sum(histFreq)*100)'];
        tb = [tb;num2cell(tbb)];
        tb = [['Binwidth' ; bandwidth'] tb];
        xlswrite([savename,'_histdata.xls'], tb, 'Histogram', 'A1');
        
        saveas(figure(4),[savename,'_histogram.jpg'])
        
        %x,y velocity plot
        figure(5)
        clf
        
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
    savename = handles.pathname(1:end-4);
    %create histogram,trajectory that accumulates every cell in every videos
    histsum(pm.vidvel,pm.step,pm.lng,savename,pm.color,1,1)
    trajsum(pm,savename)
end