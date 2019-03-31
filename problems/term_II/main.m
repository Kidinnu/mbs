
p.n  = 3;

p.S0 = [1 0 0];

p.S  = [-1  1  1
         0 -1  0
         0  0 -1];

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
p.I    = cells(p.n,1);
p.I{1} = diag([1,1,1]);
p.I{2} = diag([1,1,1]);
p.I{3} = diag([1,1,1]);
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
p.A{1} = @(q11,q12,q13) Az(q11)*Ax(q12)*Az(q13);
% 1-2 joint (universal)
p.A{2} = @(q21,q22) Az(q21)*Ax(q22);
% 2-3 joint (cylindrical)
p.A{3} = @(q31) Ax(q31);

% p-vectors
% 0-1 joint: 3 x 3
p.p{1} = @(q11,q12,q13) [p.A{1}(q11,q12,q13)*[0;0;1] p.A{1}(0,q12,q13)*[1;0;0] [0;0;1]];
% 1-2 joint: 3 x 2
p.p{2} = @(q21,q22) [p.A{2}(q21,q22)*[0;0;1] [1;0;0]];
% 2-3 joint: 3 x 1
p.p{3} = @(q31) [1;0;0];

% Relative angular velocity
p.Wr{1} = @(q11,q12,q13,dq11,dq12,dq13)  p.p{1}(q11,q12,q13)*[dq11;dq12;dq13];
p.Wr{2} = @(q21,q22,dq21,dq22)  p.p{2}(q21,q22)*[dq21;dq22];
p.Wr{3} = @(q31,dq31)  p.p{3}(q31)*dq31;

% Partial derivatives of p vectors
p.pw{1} = @(q11,q12,q13,dq11,dq12,dq13) cross(p.Wr{1}(q11,q12,q13,dq11,dq12,dq13),p.p{1}(q11,q12,q13));
p.pw{2} = @(q21,q22,dq21,dq22) cross(p.Wr{2}(q21,q22,dq21,dq22),p.p{2}(q21,q22)); 
p.pw{3} = @(q31,dq31) cross(p.Wr{3}(q31,dq31),p.p{3}(q31)); 

% Diagonal blocks of K tensor in i-CS
p.Kii = cell(p.n,1);
for i=1:p.n
    p.Kii{i} = p.I{i};
    for k=1:p.n
        p.Kii{i} = p.Kii{i} + p.mass(k)*(p.d(:,i,k)'*p.d(:,i,k)*eye(3) - p.d(:,i,k)*p.d(:,i,k)');
    end
end


%%






