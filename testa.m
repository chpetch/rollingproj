ori = pic(1).cdata;
se = strel('disk', 40);
opic = imopen(ori,se);
cpic = imclose(ori,se);
subplot(2,2,1)
imshow(opic)
subplot(2,2,2)
imshow(cpic)
subplot(2,2,3)
imshow(ori-opic)
subplot(2,2,4)
imshow(cpic - ori)

gdpic = imdilate(procpic, se) - imerode(procpic, se);
gdpic2 = imfill(gdpic,'holes');
gdpic3 = medfilt2(gdpic2,[3,3]);

cn = corner(cnpic);
imshow(vpic);
hold on
plot(cn(:,1), cn(:,2), 'r*');