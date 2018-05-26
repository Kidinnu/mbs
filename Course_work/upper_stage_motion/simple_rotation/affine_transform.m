function [x1,y1,z1]=affine_transform(x,y,z,r,A)
data = zeros(size(x,1),size(x,2),3);
data(:,:,1) = x;
data(:,:,2) = y;
data(:,:,3) = z;
data = permute(data,[3 2 1]);
r = reshape(r,3,1);
data(:,:,1) = A*data(:,:,1) + ...
              repmat(r,1,size(data(:,:,1),2));
data(:,:,2) = A*data(:,:,2) + ...
              repmat(r,1,size(data(:,:,1),2));
data = permute(data,[3 2 1]);
x1 = data(:,:,1);
y1 = data(:,:,2);
z1 = data(:,:,3);


