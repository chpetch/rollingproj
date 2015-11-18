function output = hessianridgefilter(input,pt)

[gx, gy] = gradient(double(input));
[gxx, ~] = gradient(gx);
[gxy, gyy] = gradient(gy);
output = ((gxx + gyy)+sqrt(gxx.^2 + 4.*(gxy.^2) - 2.*(gxx.*gyy) + gyy.^2))./2;

if pt == 1
    figure
    imshow(output)
end

end