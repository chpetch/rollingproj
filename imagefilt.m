function gradmag = imagefilt(I,ftype,a,b)

if strcmp(ftype,'log') == 1
    hy = fspecial(ftype,a,b);
else
    hy = fspecial(ftype);
end

hx = hy';
Iy = imfilter(double(I), hy, 'replicate');
Ix = imfilter(double(I), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);
gradmag = uint16(gradmag);

end