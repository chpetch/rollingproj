function auc = linerotation(x,y,plf)
% define the x- and y-data for the original line we would like to rotate
angle = atan2(y(end)-y(1),x(end)-x(1));
% create a matrix of these points, which will be useful in future calculations
v = [x';y'];
x_center = x(1);
y_center = y(1);
% create a matrix which will be used later in calculations
center = repmat([x_center; y_center], 1, length(x));
% define pi-angle degree counter-clockwise rotation matrix
theta = pi-angle;       % pi/3 radians = 60 degrees
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
% do the rotation...
s = v - center;     % shift points in the plane so that the center of rotation is at the origin
so = R*s;           % apply the rotation about the origin
vo = so + center;   % shift again so the origin goes back to the desired center of rotation
% this can be done in one line as:
% vo = R*(v - center) + center
% pick out the vectors of rotated x- and y-data
x_rotated = vo(1,:);
y_rotated = vo(2,:);
% make a plot
if plf == 1
    figure(4)
    plot(abs(x), y, 'k-',abs(x_rotated), y_rotated, 'r-');
    ylim auto
    axis equal
end
%find area under curve
auc = trapz(abs(x_rotated),abs(y_rotated));
end