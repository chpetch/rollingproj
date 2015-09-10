function [t_cell,sframe] = calculatevelocity(mf,sidis,pm)
% %calculate path of cells and label cells
% %input
%  sidis - cell containing matching index
%  mf - cell containing feature vectors
%  pm - structure containing parameters
% %output
%  t_cell - cell containing position of all cell in all frame
%  sframe - array containing cells' starting frame
 
t_cell{1} = [];
sframe(1) = 1;
figure(1)

for inow = 1 : size(mf,2) - 1
    
    inow
    for i = 1:size(sidis{inow},1)
        if sidis{inow}(i) ~= 0
            if abs((mf{1,inow+1}(i,2))-mf{1,inow}(sidis{inow}(i),2)) >= 0 || abs(mf{1,inow+1}(i,1)-mf{1,inow}(sidis{inow}(i),1)) >= 0
                ctf1 = [(mf{1,inow}(sidis{inow}(i),2)),(mf{1,inow}(sidis{inow}(i),1))];
                ctf2 = [(mf{1,inow+1}(i,2)),(mf{1,inow+1}(i,1))];
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
