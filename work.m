a = imgaussfilt(thpic,2)
I = a;
hy = fspecial('laplacian',1);
hx = hy';
Iy = imfilter(double(I), hy, 'replicate');
Ix = imfilter(double(I), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);
gradmag = uint16(gradmag);
figure
imshow(gradmag,[]), title('Gradient magnitude (gradmag)')
