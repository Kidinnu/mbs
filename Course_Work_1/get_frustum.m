function [x, y, z, c] = get_frustum(radius1,radius2,length, color)

n = 32;
k = linspace(0,2*pi,n);

y = [radius1*cos(k);radius2*cos(k)];
z = [radius1*sin(k);radius2*sin(k)];
x = [zeros(1,n); length*ones(1,n)];

c = zeros(2,n,3);

c(:,:,1) = color(1); c(:,:,2) = color(2); c(:,:,3) = color(3);


end