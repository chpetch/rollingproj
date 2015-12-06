function result_edited = delcell(result,delcellarray)
result_edited = result;
result_edited.cse_new(:,delcellarray) = [];
result_edited.cpos_new(delcellarray) = [];
end
