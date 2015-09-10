function sidis = matchingframes(mf)
% %output
% sidis - cell containing matching index
% %input
% mf - cell containing feature vectors

run = 1;
d.index = [];
d.amp = [];
z.index = [];
for i = 1 : size(mf,2) - 1
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
        
        %create output for this function
        sidis{i} = idis;
        
    end
end