function [result,ntf,sframe] = constrainsforcells(t_cell,result,sframe,pm)
%constrains that applied to the cell
% %input
% t_cell - cell containing position of all cell in all frame
% sframe - array containing cells' starting frame
% %output
% result - cpos_original/cpos_new - cells' centroid before / after applying constrain
%          cse_original/cse_new - cells' info label | starting frame | ending frame before / after applying constrain
%

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
        %check the ending frame of cell 1 and starting frame of cell 2 if
        %they are the same or less than pm.ccf will proceed to the next step
        if (result.cse_original(2,ccc) - result.cse_original(3,cc)) > 0 ...
                && (result.cse_original(2,ccc) - result.cse_original(3,cc)) < pm.ccf
            %check y distance between 2 cells if less than pm.d_ydis (1) proceed to the
            %next step
            if abs(result.cpos_original{1,cc}(end,2) - result.cpos_original{1,ccc}(1,2)) < pm.d_ydis
                %determine whether cell move forward or backward in x direction if
                %backward distance less than pm.d_threshold_bwd 
                %or forward distance less than pm.d_threshold
                %connect the cell
                if result.cpos_original{1,cc}(end,1) - result.cpos_original{1,ccc}(1,1) > -pm.d_threshold_bwd ...
                        && result.cpos_original{1,cc}(end,1) - result.cpos_original{1,ccc}(1,1) < pm.d_threshold                   cccs = [cccs;cc,ccc];
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
    
    %delete cells that move faster than pm.d_threshold
    for i = 1:size(t_cell,2)
        for j = 1 : size(t_cell{i},1)-1
            acvel{i}(j) = (sqrt(sum((t_cell{i}(j,:) - t_cell{i}(j+1,:)).^2)));
            if acvel{i}(j) > pm.d_threshold
                acvel{i}(j) = [];
                t_cell{i}(j+1:end,:) = [];
                disp(['cell moving faster than ',num2str(pm.d_threshold),' is detected'])
                break
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
end