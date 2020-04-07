function [dq, Rx, Ry]= dqdt_FairSep(t, q, P)

r = q(1:2);     %Position vector
phi = q(3);     %Angular position vector
v = q(4:5);     %Velocity vector
w = q(6);       %Angular velocity vector

%Shorthand
r0 = P.r0;  ra = P.ra;  rp = P.rp;
ax = P.ax;  m = P.m;    dh = P.dh;
J = P.J;    P0 = P.P0;  Pk = P.Pk;  
Stage = P.Stage;

%Transformation matrix (from Fairing-fixed CS to Hinge-fixed CS)
A = @(a) [cos(a), -sin(a); sin(a), cos(a)];

%Tilde matrix for cartesian components  
tilde2 = @(p) [-p(2), p(1)];

AP0 = -ra + A(0)*(rp - r0);     %Action vector of pusher-spring at t = 0
h0 = norm(AP0);

AP = -ra +A(phi)*(rp - r0);     %Action vector of pusher-spring
h = norm(AP);
ep = AP/h;                      %Direction vector of pusher-spring

P = P0 - (h - h0)*(P0 - Pk)/dh;     %Instantaneous reaction force due to pusher-spring

if h > h0 + dh
    P = 0;                          %Constraining action force due to pusher-spring
end

Iner = -[1; 0]*m*ax;                %Inertial force acting on fairing

Mp = tilde2(rp)*ep*P;               %Moment due to pusher-spring reaction force

W = [0, -1; 1, 0];                  %tilde matrix for angular components

Q = [eye(2), A(phi)*W*r0];          %Constraint equations

QT = -[eye(2); tilde2(A(phi)*r0)];  %Transpose of constraint matrix      !!!!! 

if Stage == 2
    M = [diag([m, m, J]), QT;
         zeros(2,3),        eye(2)];       %Coefficient matrix           !!!!!
end
if Stage == 1
    M = [diag([m, m, J]), QT;
         Q,   zeros(2)];    %Coefficient matrix
end


B = [ep*P + Iner; Mp; -A(phi)*W*W*r0*w*w];   %Force matrix    !!!!!! -A(phi)...

X = M\B;                            
%a = X(1:2), Acceleraation 
%d2phi = X(3), Angular acceleration
Rx = X(4);% Reaction force along X0 
Ry = X(5);% Reaction force along Y0

dq = [v; w; X(1:3)];                %Differential equivalence     

end

