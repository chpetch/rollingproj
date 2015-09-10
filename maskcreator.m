mask = zeros(procpic);

mask = medfilt2(I2,[3,3]);

se = strel(ones(3,3));
gdpic = imdilate(mask, se) - imerode(mask, se);
imshow(gdpic)