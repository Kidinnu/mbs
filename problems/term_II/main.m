
p.n  = 3;
p.na = [3;2;1];

p.S0 = [1 0 0];
p.S  = [-1  1  1
         0 -1  0
         0  0 -1];

p.N  = sum(p.na);
p.iq = [(cumsum(p.na)-p.na+1) p.na];
p.T  = inv(p.S);  

% Joint vectors
p.C  = zeros(3,p.n,p.n);
p.C(:,1,1)  = [-1; 0;0];
p.C(:,1,2)  = [+1; 0;0];
p.C(:,1,3)  = [+1; 0;0];
p.C(:,2,2)  = [-1; 0;0];
p.C(:,3,3)  = [-1; 0;0];

% Masses
p.mass = [1;1;1];
% Inertia tensors
p.I    = zeros(3,3,p.n);
p.I(:,:,1) = diag([1,1,1]);
p.I(:,:,2) = diag([1,1,1]);
p.I(:,:,3) = diag([1,1,1]);
% d_{ji} column vectors in j-CS
p.d  = zeros(3,p.n,p.n);
for i = 1:p.n
    for j = 1:p.n
        for a = 1:p.n
            p.d(:,j,i) = p.d(:,j,i) + p.T(a,i)*p.S(j,a)*p.C(:,j,a);
        end
    end
end
% b_{i0} column vectors in i-CS
p.b  = zeros(3,p.n);
for i = 1:p.n
    for j=1:p.n
        p.b(:,i) = p.b(:,i) + p.d(:,i,j)*p.mass(j);
    end
    p.b(:,i) = p.b(:,i)/sum(p.mass);
end
%
% Kinematics
%
p.A  = cell(p.n,1);
p.pa = cell(p.n,1);
p.pw = cell(p.n,1);

Ax = @(x) [1 0 0; 0 cos(x) -sin(x); 0 sin(x) cos(x)];
Ay = @(x) [cos(x) 0 sin(x); 0 1 0; - sin(x) 0 cos(x)];
Az = @(x) [cos(x) -sin(x) 0; sin(x) cos(x) 0; 0 0 1];

% 0-1 joint (spherical)
p.A{1} = @(q) Az(q(1))*Ax(q(2))*Az(q(3));
% 1-2 joint (universal)
p.A{2} = @(q) Az(q(1))*Ax(q(2));
% 2-3 joint (cylindrical)
p.A{3} = @(q) Ax(q(1));

% p-vectors
% 0-1 joint: 3 x 3
p.p{1} = @(q) [p.A{1}(q)*[0;0;1] p.A{1}([0,q(2),q(3)])*[1;0;0] [0;0;1]];
% 1-2 joint: 3 x 2
p.p{2} = @(q) [p.A{2}(q)*[0;0;1] [1;0;0]];
% 2-3 joint: 3 x 1
p.p{3} = @(q) [1;0;0];

% Relative angular velocity
p.Wr{1} = @(q,dq)  p.p{1}(q)*[dq(1);dq(2);dq(3)];
p.Wr{2} = @(q,dq)  p.p{2}(q)*[dq(1);dq(2)];
p.Wr{3} = @(q,dq)  p.p{3}(q)*dq;

% Partial derivatives of p vectors
p.pw{1} = @(q,dq) cross(p.Wr{1}(q,dq),p.p{1}(q));
p.pw{2} = @(q,dq) cross(p.Wr{2}(q,dq),p.p{2}(q)); 
p.pw{3} = @(q,dq) cross(p.Wr{3}(q,dq),p.p{3}(q)); 

% Diagonal blocks of K tensor in i-CS
p.Kii = zeros(3,3,p.n);
for i=1:p.n
    p.Kii(:,:,i) = p.I(:,:,i);
    for k=1:p.n
        p.Kii(:,:,i) = p.Kii(:,:,i) + p.mass(k)*(p.d(:,i,k)'*p.d(:,i,k)*eye(3) - p.d(:,i,k)*p.d(:,i,k)');
    end
end

squeeze( dqdt(0,[1;2;3;4;5;6;1;2;3;4;5;6],p) )


%%






