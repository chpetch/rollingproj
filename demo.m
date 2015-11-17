%% Correcting Nonuniform Illumination
% This example shows how to correct nonuniform illumination in an image to make
% it easy to identify individual grains of rice in the image. You can then learn
% about the characteristics of the grains and easily compute statistics for all
% the grains in the image.
%
% Copyright 1993-2014 The MathWorks, Inc.

%% Step 1: Read Image

I = imcomplement(procpic);
imshow(I)

%% Step 2: Use Morphological Opening to Estimate the Background
% Notice that the background illumination is brighter in the center of the
% image than at the bottom. Use |imopen| to estimate the background
% illumination.

background = imopen(I,strel('disk',7));
imshow(background)
% Display the Background Approximation as a Surface
figure
surf(double(background(1:16:end,1:16:end))),zlim([0 max(background(:))]);
ax = gca;
ax.YDir = 'reverse';

%% Step 3: Subtract the Background Image from the Original Image

I2 = I - background;
imshow(I2)

%%
% Note that step 2 and step 3 together could be replaced by a single step
% using |imtophat| which first calculates the morphological opening and then
% subtracts it from the original image.
%
% I2 = imtophat(I,strel('disk',15));

%% Step 4: Increase the Image Contrast

I3 = imadjust(I2);
imshow(I3);

%% Step 5: Threshold the Image
% Create a new binary image by thresholding the adjusted image. Remove
% background noise with |bwareaopen|.

level = graythresh(I3);
bw = im2bw(I3,level);
bw = imfill(bw,'holes');
bw = bwareaopen(bw, 150);
imshow(bw)

%% Step 6: Identify Objects in the Image
% The function |bwconncomp| finds all the connected components (objects) in
% the binary image. The accuracy of your results depend on the size of the
% objects, the connectivity parameter (4,8,or arbitrary), and whether or
% not any objects are touching (in which case they may be labeled as one
% object). Some of the rice grains in |bw| are touching.

cc = bwconncomp(bw, 4)

%% Step 7: Examine One Object
% Each distinct object is labeled with the same integer value. Show the
% grain that is the 50th connected component.



%% Step 8: View All Objects
% One way to visualize connected components is to create a label matrix and
% then display it as a pseudo-color indexed image.
%
% Use  |labelmatrix| to create a label matrix from the output of
% |bwconncomp|. Note that |labelmatrix| stores the label matrix in the
% smallest numeric class necessary for the number of objects.

labeled = labelmatrix(cc);
whos labeled
