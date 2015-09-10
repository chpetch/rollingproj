function proj=projstack(stack)
% PROJSTACK Methods for projecting 3-D image stacks into 2-D
%
% IN:   stack         x*y*z matrix, 3-D stack of images
%
% This Matlab function is part of the supplementary material for the 
% article: 
%
% "Bright field microscopy as an alternative to whole cell staining in % automated analysis of macrophage images" by J. Selinummi et al. 
%
% We kindly request you to acknowledge the authors properly 
% (citation or request for permission from the authors) when using this
% function.
%
% Website: http://sites.google.com/site/brightfieldorstaining
%
%
% Author: Jyrki Selinummi <jyrki.selinummi@tut.fi>

stack=double(stack);


% MAD is the only projection requiring some processing time, now
% implemented just using for loops
% If you don´t necessarily need this method, comment out for significantly
% faster processing
proj.mad=zeros(size(stack,1),size(stack,2));
proj.mad_mean=zeros(size(stack,1),size(stack,2));
for iter=1:size(stack,1)
    for iter2=1:size(stack,2)
        proj.mad(iter,iter2)=mad(stack(iter,iter2,:),1);
        proj.mad_mean(iter,iter2)=mad(stack(iter,iter2,:),0);
    end
end
proj.mad=scale_image(proj.mad);
proj.mad_mean=scale_image(proj.mad_mean);


% IQR
proj.iqr=scale_image(iqr(stack,3));
% average projection
proj.mean=scale_image(mean(stack,3));
% STD projection
proj.std=scale_image(std(stack,0,3));
% COV projection
proj.cov=scale_image(std(stack,0,3)./mean(stack,3));
% Scaled COV projection
proj.cov_scaled=imadjust(std(stack,0,3)./mean(stack,3),[0.01 0.99]);


function out=scale_image(in)
% Scale data between 0 and 1
in=double(in);
in=in-min(in(:));
out=in/max(in(:));

