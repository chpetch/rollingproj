function seqcm = tracecellmotion(pic,cpos,startingframe,plf,savepath)
for jj = 1:size(cpos,1)
    
    cmpic = uint16(zeros(512,512));
    corner(1) = cpos(jj,1)-20;
    corner(2) = cpos(jj,1)+20;
    corner(3) = cpos(jj,2)-20;
    corner(4) = cpos(jj,2)+20;
    
    if corner(1) < 1
        corner(1) = 1;
    end
    if corner(3) < 1
        corner(3) = 1;
    end
    if corner(2) > 512
        corner(2) = 512;
    end
    if corner(4) > 512
        corner(4) = 512;
    end
    repim = pic(startingframe+jj-1).cdata(corner(3):corner(4),corner(1):corner(2));
    cmpic(corner(3):corner(4),corner(1):corner(2)) = repim;
    if plf == 1
        figure(3)
        imshow(cmpic)
        saveas(figure(3),[savepath,'_',num2str(jj),'.jpg'])        
        pause(0.5)
    end
    seqcm{jj} = cmpic;
end