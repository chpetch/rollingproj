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
        procpic = pic(i).cdata;
        procpic = imcomplement(procpic);

        I2 = imtophat(procpic, strel('disk', 5));
        
        procpic = medfilt2(procpic);
        procpic = procpic*(65535/double(max(procpic(:))));
        imshow(procpic)
        
        %
        %procpic = bwareaopen(procpic,20);
        %         CC = bwconncomp(im2bw(procpic), 8);
        obj = regionprops(im2bw(procpic),'centroid','EquivDiameter','ConvexArea','MajorAxisLength','MinorAxisLength');
        %         center = cat(1, obj.Centroid);
        %         eqd = cat(1, obj.EquivDiameter);
        convexarea = cat(1, obj.ConvexArea);
        
        %         mjmn = [cat(1, obj.MajorAxisLength) cat(1, obj.MinorAxisLength)];
        %         if isempty(mjmn) == 1
        %             disp('error spotted > initializing self destruction')
        %             pause
        %         end
        %         disp(['Before > number of obj detected : ',num2str(length(center))])
        %
        %         %determining circles
        %         clear rmjmn
        %         rmjmn = mjmn(:,1) ./ mjmn(:,2);
        %         ncc = find(rmjmn>2);
        %         if isempty(find(rmjmn>2)) ~= 1
        %             center(ncc,:) = [];
        %             convexarea(ncc) = [];
        %         end
        %         %determining area
        %          if isempty(convexarea) ~= 1
        % %             convexarea(find(convexarea<15)) = [];
        %              convexarea(find(convexarea>100)) = [];
        %              center(find(convexarea>100)) = [];
        %          end
        
        % disp(['After > number of obj detected : ',num2str(length(center))])
        
        %centers{i} = fliplr(center);
        %     imshow(procpic)
        %     plot(center(:,1), center(:,2), 'b*');
        %     pause(0.5)
        %     %     %procpic = imadjust(procpic);
        [center, radii, metric] = imfindcircles(procpic,[10,50]);
        %     %     %     imshow(pic2(i).cdata)
        hold on
        viscircles(center, radii,'EdgeColor','b');
        
        pause
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
                
                ang{i}(j) = mod(( atan2([round(mf{1,i}(idis(j),2))-round(mf{1,i+1}(j,2))], ...
                    [round(mf{1,i}(idis(j),1))-round(mf{1,i+1}(j,1))]) + 360 ),360);
                
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
        %     imshow(pic(inow).cdata)
        %     hold on
        inow
        for i = 1:size(sidis{inow},1)
            if sidis{inow}(i) ~= 0
                if abs((mf{1,inow+1}(i,2))-mf{1,inow}(sidis{inow}(i),2)) > 0 || abs(mf{1,inow+1}(i,1)-mf{1,inow}(sidis{inow}(i),1)) > 0
                    clear vcell
                    ctf1 = [(mf{1,inow}(sidis{inow}(i),2)),(mf{1,inow}(sidis{inow}(i),1))];
                    ctf2 = [(mf{1,inow+1}(i,2)),(mf{1,inow+1}(i,1))];
                    vcell = (sqrt(sum((ctf2 - ctf1).^2))*pm.umperpx)/pm.tbf;
                    %draw circle on cell in frame 1
                    %                 plot(ctf1(1), ctf1(2), 'ro', 'MarkerSize', 12)
                    %                 hold on
                    %                 %draw circle on cell in frame 2
                    %                 plot(ctf2(1), ctf2(2), 'go', 'MarkerSize', 13)
                    %                 hold on
                    %draw line , trajectory
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
                    %                 plot([ctf2(1),ctf1(1)],[ctf2(2),ctf1(2)], 'b', 'LineWidth', 1)
                    %                 hold on
                else
                    %draw *
                    %plot((mf{1,inow}(sidis{inow}(i),2)), (mf{1,inow}(sidis{inow}(i),1)), 'b*', 'MarkerSize', 12)
                end
            else
                disp('Tada!!!')
            end
            %         hold on
            %         legend(['Position of cell in current frame (',num2str(inow),')'],['Future position of cell (',num2str(inow+1),')'],'Displacement','Non-moving cell','Location','southoutside')
        end
        %   pause(0.1)
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
        %histFreq(1) = histFreq(1) - sum(cell2mat(handles.nsv));
        bar(histXout, histFreq/sum(histFreq)*100,1,pm.color);
        xlabel('Rolling cell velocity (um/s)');
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
        %         for i = 1:size(cvelxy,2)
        %             plot(cvelxy{i}(:,1),cvelxy{i}(:,2),'*','color',rand(1,3))
        %             hold on
        %         end
        %         saveas(figure(5),[savename,'_velscatterplot.jpg'])
        
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

savename = handles.pathname(1:end-4);
%create histogram,trajectory that accumulates every cell in every videos
histsum(pm.vidvel,pm.step,pm.lng,savename,pm.color)
trajsum(pm,savename)

