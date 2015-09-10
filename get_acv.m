function acv = get_acv(ipt)

acv = [];

for iii = 1:length(ipt)
    iii
    for jjj = 1:length(ipt{iii})
        jjj
        if mean(ipt{iii}{jjj}) == 0
            disp('hey')
        else
            acv = [acv;mean(ipt{iii}{jjj})];            
        end
    end
end