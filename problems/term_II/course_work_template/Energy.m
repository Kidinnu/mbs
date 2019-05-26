function [E, T, U] = Energy(q, p)
%% ========================================================================
% Angular velocity
w = zeros(3, p.n, size(q,1));

for k=1:size(q,1)
    phi  = @(i) q(k,p.iq(i,1) : p.iq(i,1)+p.iq(i,2)-1);
    dphi = @(i) q(k,p.iq(i,1)+p.N : p.iq(i,1)+p.iq(i,2)+p.N-1);
    
    A0i = zeros(3,3,p.n);
    for i=1:p.n
        A0i(:,:,i) = eye(3);
        for a=i:-1:1        
            if p.T(a,i) ~= 0
                A0i(:,:,i) = p.A{a}(phi(a))*A0i(:,:,i);
            end
        end
    end
    
    for i=1:p.n
        for j=1:p.n
            w(:,i,k) = w(:,i,k) - A0i(:,:,j)*(p.T(j,i)*p.Wr{j}(phi(j),dphi(j)));
        end
    end
end

%% ========================================================================
% Linear velocity 
v1 = cell(size(q,1),1);
v2 = cell(size(q,1),1);
v3 = cell(size(q,1),1);
v4 = cell(size(q,1),1);
v5 = cell(size(q,1),1);

for k=1:size(q,1)
    phi  = @(i) q(k,p.iq(i,1) : p.iq(i,1)+p.iq(i,2)-1);
    dphi = @(i) q(k,p.iq(i,1)+p.N : p.iq(i,1)+p.iq(i,2)+p.N-1);
    % transfomraion matrix (body i to 0)
    A0i = zeros(3,3,p.n);   
    for j=1:p.n
        A0i(:,:,j) = eye(3);
        for a=j:-1:1        
            if p.T(a,j) ~= 0
                A0i(:,:,j) = p.A{a}(phi(a))*A0i(:,:,j);
            end
        end
    end
    % d_{ij} vectors
    d011 = A0i(:,:,1)*p.d(:,1,1);
    d012 = A0i(:,:,1)*p.d(:,1,2);
    d022 = A0i(:,:,2)*p.d(:,2,2);
    d023 = A0i(:,:,2)*p.d(:,2,3);
    d033 = A0i(:,:,3)*p.d(:,3,3);
    d034 = A0i(:,:,3)*p.d(:,3,4);
    d035 = A0i(:,:,3)*p.d(:,3,5);
    d044 = A0i(:,:,4)*p.d(:,4,4);
    d055 = A0i(:,:,5)*p.d(:,5,5);
    % auxiliary angular velocity
    Wr = cell(5,1);
    for i=1:p.n
        Wr{i} = p.Wr{i}(phi(i),dphi(i));
    end
    w1 = A0i(:,:,1)*Wr{1};
    w2 = w1 + A0i(:,:,2)*Wr{2};
    w3 = w2 + A0i(:,:,3)*Wr{3};
    w4 = w3 + A0i(:,:,4)*Wr{4};
    w5 = w3 + A0i(:,:,5)*Wr{5};
    % linear velocity
    v1{k} = cross(w1,-d011);
    v2{k} = cross(w1,-d012) + cross(w2,-d022);
    v3{k} = cross(w1,-d012) + cross(w2,-d023) + cross(w3,-d033);
    v4{k} = cross(w1,-d012) + cross(w2,-d023) + cross(w3,-d034) + cross(w4,-d044);
    v5{k} = cross(w1,-d012) + cross(w2,-d023) + cross(w3,-d035) + cross(w5,-d055);
end

%% ========================================================================
% Kinetic energy
Tw = zeros(size(q,1),1);    % rotational component
for k=1:size(q,1)
    phi  = @(i) q(k,p.iq(i,1) : p.iq(i,1)+p.iq(i,2)-1);
    % transfomraion matrix (body i to 0)
    A0i = zeros(3,3,p.n);   
    for j=1:p.n
        A0i(:,:,j) = eye(3);
        for a=j:-1:1        
            if p.T(a,j) ~= 0
                A0i(:,:,j) = p.A{a}(phi(a))*A0i(:,:,j);
            end
        end
    end
    % inertia tensor in body-0 frame
    I01 = A0i(:,:,1)*p.I(:,:,1)*A0i(:,:,1)';
    I02 = A0i(:,:,2)*p.I(:,:,2)*A0i(:,:,2)';
    I03 = A0i(:,:,3)*p.I(:,:,3)*A0i(:,:,3)';
    I04 = A0i(:,:,4)*p.I(:,:,4)*A0i(:,:,4)';
    I05 = A0i(:,:,5)*p.I(:,:,5)*A0i(:,:,5)';
    % individual rotational energies
    Tw_1 = 0.5*w(:,1,k)'*I01*w(:,1,k); 
    Tw_2 = 0.5*w(:,2,k)'*I02*w(:,2,k); 
    Tw_3 = 0.5*w(:,3,k)'*I03*w(:,3,k); 
    Tw_4 = 0.5*w(:,4,k)'*I04*w(:,4,k); 
    Tw_5 = 0.5*w(:,5,k)'*I05*w(:,5,k); 
    % total rotationl energy
    Tw(k) = Tw_1 + Tw_2 + Tw_3 + Tw_4 + Tw_5;
end

Tl = zeros(size(q,1),1);    % linear component
for k=1:size(q,1)
    v_1 = norm(v1{k});
    v_2 = norm(v2{k});
    v_3 = norm(v3{k});
    v_4 = norm(v4{k});
    v_5 = norm(v5{k});
    Tl(k) = 0.5*p.mass(1)*v_1*v_1 + 0.5*p.mass(2)*v_2*v_2 + 0.5*p.mass(3)*v_3*v_3 + 0.5*p.mass(4)*v_4*v_4 + 0.5*p.mass(5)*v_5*v_5;
end

% total kinetic energy
T = Tw + Tl;

%% ========================================================================
% Potential energy
V = zeros(size(q,1),1);
for k=1:size(q,1)
    phi  = @(i) q(k,p.iq(i,1) : p.iq(i,1)+p.iq(i,2)-1);
    % transfomraion matrix (body i to 0)
    A0i = zeros(3,3,p.n);   
    for i=1:p.n
        A0i(:,:,i) = eye(3);
        for a=i:-1:1        
            if p.T(a,i) ~= 0
                A0i(:,:,i) = p.A{a}(phi(a))*A0i(:,:,i);
            end
        end
    end
    % d_{ij} vectors
    d011 = A0i(:,:,1)*p.d(:,1,1);
    d012 = A0i(:,:,1)*p.d(:,1,2);
    d022 = A0i(:,:,2)*p.d(:,2,2);
    d023 = A0i(:,:,2)*p.d(:,2,3);
    d033 = A0i(:,:,3)*p.d(:,3,3);
    d034 = A0i(:,:,3)*p.d(:,3,4);
    d035 = A0i(:,:,3)*p.d(:,3,5);
    d044 = A0i(:,:,4)*p.d(:,4,4);
    d055 = A0i(:,:,5)*p.d(:,5,5);
    % distance of body i from body 0
    h1 = - d011;
    h2 = - d012 - d022;
    h3 = - d012 - d023 - d033;
    h4 = - d012 - d023 - d034 - d044;
    h5 = - d012 - d023 - d035 - d055;
    % sum of potential energies of 5 bodies
    V(k) = 9.807*(p.mass(1)*h1(3) + p.mass(2)*h2(3) + p.mass(3)*h3(3) + p.mass(4)*h4(3) + p.mass(5)*h5(3));
end

% potential energy
U = V;

%% ========================================================================
% Total energy
E = T + U;

end