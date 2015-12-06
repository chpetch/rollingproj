        for i = 1:size(result.cpos_new,2)
            for j = 1 : size(result.cpos_new{i},1)-1
                cvel{i}(j) = (sqrt(sum((result.cpos_new{i}(j,:) - result.cpos_new{i}(j+1,:)).^2))*pm.umperpx)/pm.tbf;
                cvelxy{i}(j,:) = (result.cpos_new{i}(j,:) - result.cpos_new{i}(j+1,:)).*(pm.umperpx/pm.tbf);
            end
        end