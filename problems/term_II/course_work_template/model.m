%% ========================================================================
% Описание модели
% =========================================================================
clc;
% Количество тел
p.n  = 5;

% Колчиество степеней свободы в каждом шарнире (от 1 до 3)
p.na = [3;3;1;1;1];

% Матрица инцидентности 
p.S0 = [1 0 0 0 0];
p.S  = [-1  1  0  0  0    
         0 -1  1  0  0
         0  0 -1  1  1
         0  0  0 -1  0
         0  0  0  0 -1];

% Шарнирные вектора
p.C  = zeros(3,p.n,p.n);
p.C(:,1,1)  = [0;0;+1];
p.C(:,1,2)  = [0;0;-1];
p.C(:,2,3)  = [0;0;-1];
p.C(:,2,2)  = [0;0;+1];
p.C(:,3,3)  = [0;0;+1];
p.C(:,3,4)  = [+0.5;0;0];
p.C(:,3,5)  = [-0.5;0;0];
p.C(:,4,4)  = [0;0;+0.5];
p.C(:,5,5)  = [0;0;+0.5];

% Masses
p.mass = [3;3;2;1.5;1.5];       

% Inertia tensors 
p.I    = zeros(3,3,p.n);
p.I(:,:,1) = p.mass(1)*diag([4/3,4/3,0.6*4/2]);     % Iz = 60% of (Ix,Iy)
p.I(:,:,2) = p.mass(2)*diag([4/3,4/3,0.6*4/3]);     % Iz = 60% of (Ix,Iy)
p.I(:,:,3) = p.mass(3)*diag([1/6,5/24,1/24]);
p.I(:,:,4) = p.mass(4)*diag([1/3,1/3,0.6*1/3]);     % Iz = 60% of (Ix,Iy)
p.I(:,:,5) = p.mass(5)*diag([1/3,1/3,0.6*1/3]);     % Iz = 60% of (Ix,Iy)

% Kinematics
p.A  = cell(p.n,1);
p.p  = cell(p.n,1);
p.pw = cell(p.n,1);

% Basic rotations
Ax = @(x) [1 0 0; 0 cos(x) -sin(x); 0 sin(x) cos(x)];
Ay = @(x) [cos(x) 0 sin(x); 0 1 0; - sin(x) 0 cos(x)];
Az = @(x) [cos(x) -sin(x) 0; sin(x) cos(x) 0; 0 0 1];

% 0-1 (1 to 0) joint (spherical)
p.A{1} = @(q) Az(q(1))*Ax(q(2))*Az(q(3));

% 1-2 (2 to 1) joint (spherical)
p.A{2} = @(q) Az(q(1))*Ax(q(2))*Az(q(3));

% 2-3 (3 to 2) joint (cylindrical)
p.A{3} = @(q) Az(q);

% 3-4 (4 to 3) joint (cylindrical)
p.A{4} = @(q) Ay(q);

% 3-5 (5 to 3) joint (cylindrical)
p.A{5} = @(q) Ay(q);

% p-vectors
% 0-1 spherical joint: 3 x 3 
%p.p{1} = @(q) [[sin(q(2))*sin(q(3));sin(q(2))*cos(q(3));cos(q(2))] [cos(q(3));-sin(q(3));0] [0;0;1]];
p.p{1} = @(q) [p.A{1}(q)'*[0;0;1] p.A{1}([0,q(2),q(3)])'*[1;0;0] [0;0;1]];
% 1-2 spherical joint: 3 x 3
%p.p{2} = @(q) [[sin(q(2))*sin(q(3));sin(q(2))*cos(q(3));cos(q(2))] [cos(q(3));-sin(q(3));0] [0;0;1]];
p.p{2} = @(q) [p.A{2}(q)'*[0;0;1] p.A{2}([0,q(2),q(3)])'*[1;0;0] [0;0;1]];
% 2-3 cylindrical joint: 3 x 1
p.p{3} = @(q) [0;0;1];

% 3-4 cylindrical joint: 3 x 1
p.p{4} = @(q) [0;1;0];

% 3-5 cylindrical joint: 3 x 1
p.p{5} = @(q) [0;1;0];

% Relative angular velocity
p.Wr{1} = @(q,dq)  p.p{1}(q)*[dq(1);dq(2);dq(3)];
p.Wr{2} = @(q,dq)  p.p{2}(q)*[dq(1);dq(2);dq(3)];
p.Wr{3} = @(q,dq)  p.p{3}(q)*dq;
p.Wr{4} = @(q,dq)  p.p{4}(q)*dq;
p.Wr{5} = @(q,dq)  p.p{5}(q)*dq;

% Relative angular acceleration
p.pw{1} = @(q,dq) [sin(q(3))*dq(2)*(cos(q(2))*dq(1)-dq(3))+cos(q(3))*sin(q(2))*dq(1)*dq(3);
                   cos(q(3))*dq(2)*(cos(q(2))*dq(1)-dq(3))-sin(q(3))*sin(q(2))*dq(1)*dq(3);
                   -sin(q(2))*dq(2)*dq(1)];
p.pw{2} = @(q,dq) [sin(q(3))*dq(2)*(cos(q(2))*dq(1)-dq(3))+cos(q(3))*sin(q(2))*dq(1)*dq(3);
                   cos(q(3))*dq(2)*(cos(q(2))*dq(1)-dq(3))-sin(q(3))*sin(q(2))*dq(1)*dq(3);
                   -sin(q(2))*dq(2)*dq(1)]; 

p.pw{3} = @(q,dq) [0;0;0]; 
p.pw{4} = @(q,dq) [0;0;0]; 
p.pw{5} = @(q,dq) [0;0;0]; 
% ========================================================================
% Preprocessing
% =========================================================================
preproc;
% ========================================================================
% Simulation
% =========================================================================
% Initial conditions
q0 = [pi/2;0.5;0.1;pi/2;0.2;0.1;0;0;0;0;0;0;0;0;0;0;0;0];
% Start integration process
[t, q] = ode113(@(t,q) ode_dqdt(t,q,p), [0 5], q0, odeset('RelTol',1e-7));
fprintf('OK\n');

%% ========================================================================
% Postprocessing
% =========================================================================
set(0,'defaultLineLineWidth',2);
set(0,'defaultFigureColor','w');
set(0,'defaultAxesFontSize',12);
set(0,'DefaultAxesXGrid','on');
set(0,'DefaultAxesYGrid','on');

figure;
subplot(511);
plot(t, q(:,1:3));legend('\phi_{11}','\phi_{12}','\phi_{13}');xlabel('t, s');
subplot(512);
plot(t, q(:,4:6));legend('\phi_{21}','\phi_{22}','\phi_{23}');xlabel('t, s');
subplot(513);
plot(t, q(:,7));legend('\phi_{33}');xlabel('t, s');
subplot(514);
plot(t, q(:,8));legend('\phi_{42}');xlabel('t, s');
subplot(515);
plot(t, q(:,9));legend('\phi_{52}');xlabel('t, s');
set(gca,'FontSize',12); grid on;

%% ========================================================================
% Validation : Checking conservation of the energy
% =========================================================================
figure;
[T, V] = Energy(q,p);
figure;
plot(t, T+V,'-',t, T,'--', t, V,':'); 
legend('Total energy', 'Kinetic energy', 'Potential energy'); 
xlabel('t, s'); ylabel('Energy, J');

%% ========================================================================
% Export data for animation 
% =========================================================================
csvwrite('results.csv',q)
csvwrite('vec_d.csv',p.d)
