% =========================================================================
% Model description
% =========================================================================
% Number of bodies
p.n  = 3;
% Numbers DOF for each joint
p.na = [3;2;1];
% Structure matrix 
p.S0 = [1 0 0];
p.S  = [-1  1  1
         0 -1  0
         0  0 -1];
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
% Kinematics
p.A  = cell(p.n,1);
p.pa = cell(p.n,1);
p.pw = cell(p.n,1);
% Basic rotations
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
p.pw{1} = @(q,dq) sum(tilde(p.Wr{1}(q,dq))*p.p{1}(q),2);
p.pw{2} = @(q,dq) sum(tilde(p.Wr{2}(q,dq))*p.p{2}(q),2); 
p.pw{3} = @(q,dq) [0;0;0]; 
% =========================================================================
% Preprocessing
% =========================================================================
preproc;
% =========================================================================
% Simulation
% =========================================================================
% Initial conditions
q0 = [1;1;1;1;1;1;0;0;0;0;0;0];
% Start integration process
[t, q] = ode113(@(t,q) dqdt(t,q,p), [0 1.3], q0);
%% =========================================================================
% Postprocessing
% =========================================================================
subplot(311);
plot(t, q(:,1:3));legend('\phi_{11}','\phi_{12}','\phi_{13}');xlabel('t, c');
subplot(312);
plot(t, q(:,4:5));legend('\phi_{21}','\phi_{22}');xlabel('t, c');
subplot(313);
plot(t, q(:,6));legend('\phi_{31}');xlabel('t, c');
% Checking conservation of the energy
%
%
%

% Export data for animation 
%
%
%







