clear
clc
clf

filetype{1} = '*.tif';
filetype{2} = '*.avi';

%associated variables
pm.tbf = 0.5;
pm.contrast = [0.6 1];
pm.nFrames = 61;
pm.pkheight = 0.1;
pm.tlr = 5;
pm.d_threshold = 25;
pm.mft = 10;
pm.pkdist = 5;
pm.input_binsize = 9;

[FileName,PathName] = uigetfile(filetype','Select the video file');

handles.pathname = [PathName,FileName];

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
            %vc(j) = sum(pic(i).cdata(:,j))/(512*65535);
        end
        hc = hc - min(hc(:));
        hc = hc/max(hc);
        [pk,lc] = findpeaks(hc,'MinPeakHeight',pm.pkheight,'MinPeakDistance',pm.pkdist);
        for j = 1:1:length(pk)-1
            diffpeak(j) = abs(pk(j) - pk(j+1));
        end
        [~,idx] = min(diffpeak);
        pm.upbd = lc(idx);
        pm.lwbd = lc(idx+1);
        %         zproject = double(pic(1).cdata);
        %         zproject = zproject-min(zproject(:));
        %         zproject = zproject/max(zproject(:));
        plot(hc)
        hold on
        plot(lc,pk,'rx')
        %         newimage = double(pic(i).cdata);
        %         newimage = newimage-min(newimage(:));
        %         zproject = zproject + newimage;
        pause(0.5)
    end
end
% zproject = zproject/(pm.nFrames*65535);
disp(['cutoff boundary : ' , num2str(pm.upbd) ,' - ', num2str(pm.lwbd)])
npic = std(kkk,0,3);
npic = npic/max(npic(:));
% imshow(pic(1).cdata)
% hold on
% plot([1,512],[pm.upbd,pm.upbd],[1,512],[pm.lwbd,pm.lwbd],'r-')
clearvars -except pic pm.upbd pm.lwbd pm.tbf pm.nFrames zproject npic pm handles

% %background subtraction
% diff = pic(1).cdata - pic(pm.nFrames).cdata;
% for i = 1:size(diff,2)
%     sol(i) = length(find(diff(i,:)==0));
% end
% sol = sol*512/max(sol);
% bdy = find(sol>500);
% for i = 1:length(bdy)-1
%     sbdy(i) = bdy(i+1) - bdy(i);
% end
% [~,chc] = max(sbdy);
% pm.upbd = bdy(chc);
% pm.lwbd = bdy(chc+1);
% bgtemp_up = pic(pm.nFrames).cdata(1:pm.upbd+pm.tlr,:);
% bgtemp_low = pic(pm.nFrames).cdata(pm.lwbd-pm.tlr:end,:);
%test
% demo = pic(1).cdata;
% demo(1:pm.upbd+pm.tlr,:) = 0;
% demo(pm.lwbd-pm.tlr:end,:) = 0;
% imshow(demo)

for i = 1:1:pm.nFrames
    i
    procpic = pic(i).cdata;
    procpic(1:pm.upbd+pm.tlr,:) = 0;
    procpic(pm.lwbd-pm.tlr:end,:) = 0;
    %procpic = procpic - mean(procpic(:));
    
    %watershed segmentation
    I2 = imtophat(procpic, strel('disk', 10));
    %     I2 = imadjust(I2,[0.53 0.9]);
    %     level = graythresh(I2);
    %     BW = im2bw(I2,level);
    %     D = -bwdist(~BW);
    %     D(~BW) = -Inf;
    %     L = watershed(D);
    %     subplot(1,2,1)
    %     imshow(pic(i).cdata)
    %     subplot(1,2,2)
    %     imshow(label2rgb(L,'jet','w'))
    
    procpic = procpic*(65535/double(max(procpic(:))));
    %contrast adjust
    %imtool(procpic)
    %procpic = edge(procpic)
    procpic = imadjust(procpic,pm.contrast);
    %    pause
    %     procpic = medfilt2(procpic,[3,3]);
    %     %procpic = edge(procpic);
    %     % procpic = medfilt2(procpic,[2,2]);
    %     %     procpic = im2bw(procpic,threshold);
    %
    %procpic = bwareaopen(procpic,20);
    CC = bwconncomp(im2bw(procpic), 8);
    obj = regionprops(im2bw(procpic),'centroid','EquivDiameter','ConvexArea','MajorAxisLength','MinorAxisLength');
    center = cat(1, obj.Centroid);
    eqd = cat(1, obj.EquivDiameter);
    convexarea = cat(1, obj.ConvexArea);
    mjmn = [cat(1, obj.MajorAxisLength) cat(1, obj.MinorAxisLength)];
    if isempty(mjmn) == 1
        disp('your error')
    end
    disp(['Before > number of obj detected : ',num2str(length(center))])
    
    if length(center) == 0
        
    end
    
    %determining circles
    clear rmjmn
    rmjmn = mjmn(:,1) ./ mjmn(:,2);
    ncc = find(rmjmn>2);
    if isempty(find(rmjmn>2)) ~= 1
        center(ncc,:) = [];
    end
    disp(['After > number of obj detected : ',num2str(length(center))])
    %     if length(convexarea) > 1 && isempty(convexarea) == 0
    %         [idx ,C] = kmeans(convexarea,2,'EmptyAction','singleton');
    %         [~,idm]= min(C);
    %         center(find(idx==idm),:) = [];
    %         %hold on
    %         %plot(center(:,1), center(:,2), 'b*')
    %     else
    %         disp('no obj detected')
    %     end
    centers{i} = fliplr(center);
        imshow(procpic)
        hold on
        plot(center(:,1), center(:,2), 'b*');
                pause(0.5)
    %     %     %procpic = imadjust(procpic);
    %     %     [centers, radii, metric] = imfindcircles(procpic,[0.1 10]);
    %     %     %     imshow(pic2(i).cdata)
    %     %viscircles(centers, radii,'EdgeColor','b');
    %     pause(0.5)
    
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
                ctf1 = [(mf{1,inow}(sidis{inow}(i),2)),(mf{1,inow}(sidis{inow}(i),1))];
                ctf2 = [(mf{1,inow+1}(i,2)),(mf{1,inow+1}(i,1))];
                %draw circle on cell in frame 1
                %                 plot(ctf1(1), ctf1(2), 'ro', 'MarkerSize', 12)
                %                 hold on
                %                 %draw circle on cell in frame 2
                %                 plot(ctf2(1), ctf2(2), 'go', 'MarkerSize', 13)
                %                 hold on
                %draw line , trajectory
                if inow == 1
                    if length(t_cell) < i
                        if (sqrt(sum((ctf2 - ctf1).^2)) < pm.d_threshold) && ((ctf1(1) - ctf2(1)) > 0)
                            t_cell{i} = [ctf1;ctf2];
                            sframe(i) = inow;
                        end
                    else
                        if (sqrt(sum((ctf2 - ctf1).^2)) < pm.d_threshold) && ((ctf1(1) - ctf2(1)) > 0)
                            t_cell{i} = [t_cell{i};[ctf1;ctf2]];
                        end
                    end
                else
                    if isempty(t_cell{1}) == 1
                        if (sqrt(sum((ctf2 - ctf1).^2)) < pm.d_threshold) && ((ctf1(1) - ctf2(1)) > 0)
                            t_cell{1} = [ctf1;ctf2];
                            sframe(1) = inow;
                        end
                    else
                        for j = 1:size(t_cell,2)
                            if t_cell{j}(end,:) == ctf1
                                if (sqrt(sum((t_cell{j}(end,:) - ctf2).^2)) < pm.d_threshold) && ((t_cell{j}(end,1)- ctf2(1)) > 0)
                                    t_cell{j} = [t_cell{j};ctf2];
                                    break;
                                end
                            else
                                if j == size(t_cell,2)
                                    if (sqrt(sum((ctf2 - ctf1).^2)) < pm.d_threshold) && ((ctf1(1) - ctf2(1)) > 0)
                                        t_cell{j+1} = [ctf1;ctf2];
                                        sframe(j+1) = inow;
                                    end
                                end
                            end
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

for i = 1:size(t_cell,2)
    ntf(i) = size(t_cell{i},1);
    for j = 1 : size(t_cell{i},1)-1
        cvel{i}(j) = sqrt(sum((t_cell{i}(j,:) - t_cell{i}(j+1,:)).^2))/pm.tbf;
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
for j = 1:size(result.cpos_new,2)
    dataplot = bsxfun(@minus,result.cpos_new{j},result.cpos_new{j}(1,:));
    plot(dataplot(:,1),dataplot(:,2),'color',rand(1,3), 'LineWidth', 1)
    hold on
end

saveas(figure(3),[savename,'_tj_graph.jpg'])

[histFreq, histXout] = hist(cell2mat(cvel),[0:(pm.d_threshold/pm.input_binsize):pm.d_threshold]);
figure(4)
%delete cells without pair
%histFreq(1) = histFreq(1) - sum(cell2mat(handles.nsv));
bar(histXout, histFreq/sum(histFreq)*100);
xlabel('Velocity (um/s)');
ylabel('Occurance (percent)');

%saving part%
step = pm.d_threshold/pm.input_binsize;
sv = step/2;
for i = 1 : length(histXout)
    if i == 1
        bandwidth{i} = ['0 - ',num2str(sv)];
    elseif i == length(histXout)
        bandwidth{i} = [num2str(sv),' - ','more than ',num2str(pm.d_threshold)];
    else
        bandwidth{i} = [num2str(sv),' - ',num2str(sv+step)];
        sv = sv + step;
    end
end
tb = {'Bin Center (X)','Frequency (Y)', 'Percentage'};
tbb = [histXout',histFreq',(histFreq/sum(histFreq)*100)'];
tb = [tb;num2cell(tbb)];
tb = [['Binwidth' ; bandwidth'] tb];
xlswrite([savename,'_histdata.xls'], tb, 'Histogram', 'A1');

saveas(figure(4),[savename,'_histogram.jpg'])

fpname = handles.pathname;

save([handles.pathname(1:end-4),'.mat'],'fpname','pm','result')

