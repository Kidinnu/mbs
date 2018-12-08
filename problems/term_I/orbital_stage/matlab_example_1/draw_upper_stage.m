function draw_upper_stage(bc, rc, A)

% Nozzle 
L1 = 0.7; R11 = 0.8; R12 = 0.5;
[x,y,z,c] = get_frustum(R11,R12,L1,[1,0,0]);
[x,y,z  ] = affine_transform(x-bc(1), y-bc(2), z-bc(3), rc, A);
surf(x,y,z,c);

% Rear
L2 = 0.5; R21 = 0.8; R22 = 1.5;
[x,y,z,c] = get_frustum(R21,R22,L2,[0.5,0.5,0.5]);
x = x + L1;
[x,y,z  ] = affine_transform(x-bc(1), y-bc(2), z-bc(3), rc, A);
surf(x,y,z,c);

%fill3(x(1,:),y(1,:),z(1,:),c(1,:,:));

% Body
L3 = 6; R3 = 1.5;
[x,y,z,c] = get_frustum(R3,R3,L3,[0.0,0.5,0.0]);
x = x + L1 + L2;
[x,y,z  ] = affine_transform(x-bc(1), y-bc(2), z-bc(3), rc, A);
surf(x,y,z,c);
%fill3(x(2,:),y(2,:),z(2,:),c(2,:,:));


% Adapter
L4 = 0.5; R41 = 1.5; R42 = 0.7;
[x,y,z,c] = get_frustum(R41,R42,L4,[0.5,1.0,0.5]);
x = x + L1 + L2 + L3;
[x,y,z  ] = affine_transform(x-bc(1), y-bc(2), z-bc(3), rc, A);
surf(x,y,z,c);

% Axes
xc = A*[0, 7;
        0, 0;
        0, 0] + repmat(rc,1,2);
yc = A*[0, 0;
        0, 5;
        0, 0] + repmat(rc,1,2);
zc = A*[0, 0;
        0, 0;
        0, 5] + repmat(rc,1,2);  

line(xc(1,:),xc(2,:),xc(3,:),'Color',[1,0,0]);
line(yc(1,:),yc(2,:),yc(3,:),'Color',[0,1,0]);
line(zc(1,:),zc(2,:),zc(3,:),'Color',[0,0,1]);

end





