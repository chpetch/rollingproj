function acv = get_acv(ipt,ntf)

acv = [];

for iii = 1:length(ipt)
    iii
    for jjj = 1:length(ipt{iii})
        jjj
        if sum(ipt{iii}{jjj})./(ntf{iii}(jjj)-1) == 0
            disp('hey')
        else
            acv = [acv;sum(ipt{iii}{jjj})./(ntf{iii}(jjj)-1)];            
        end
    end
end