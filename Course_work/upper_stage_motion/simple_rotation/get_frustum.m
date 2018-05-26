function [x, y, z, c] = get_frustum(R1,R2,L,cl)

n = 32;
k = linspace(0,2*pi,n);

y = [R1*cos(k);R2*cos(k)];
z = [R1*sin(k);R2*sin(k)];
x = [zeros(1,n); L*ones(1,n)];

c = zeros(2,n,3);
c(:,:,1) = cl(1); 
c(:,:,2) = cl(2); 
c(:,:,3) = cl(3);

end