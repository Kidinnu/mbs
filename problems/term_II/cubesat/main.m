%
% 2D model of the cubesat motion in the deployer
%   
% Vadim Yudintsev yudintsev@classmech.ru
% 
clc;
params = struct();

% Cubesat mass, kg
params.mass  = 1.8;
% Moment of inertia, kgm2
params.J     = params.mass*(0.340^2 + 0.1^2)/12;
% Center of mass position, m
params.yc    = 0.05    -0.02;
params.zc    = 0.340/2 -0.05;
% Cubesat length, m
params.L     = 0.340;
% Cubesat width, m
params.w     = 0.100;
% The gap between the guide rails and the cubesat, m
params.delta = 0.002;
% Pusher spring stroke, m
params.hp    = 0.339;
% Initial spring force, N
params.Fp0   = 10;
% Spring Stiffness
params.cp    = (params.Fp0-5)/params.hp;
% Final time, s
params.tk    = 0.8;

% Place for reaction forces
aux_results  = [];
% Initial conditions 
% Angle between the guide rails and the cubesat
fh = @(phi) params.w*cos(phi) + (params.L+params.w*sin(phi))*tan(phi) - params.w - params.delta;
phi0 = fzero(fh,0.001);

% Center of mass position
y0   = (params.yc-params.w)*cos(phi0) - params.zc*sin(phi0);
z0   = params.zc*cos(phi0) + (params.yc-params.w)*sin(phi0);

% Initial state vector
q0   = [y0;z0;phi0;0;0;0];

% Print initial conditions
fprintf('Initial conditions: \n');
fprintf('y0 = %5.1f mm, z0 = %5.1f mm, phi0 = %5.3f deg \n-----\n', y0*1000, z0*1000, phi0*180/pi);

% For t=0 the cubesat contacts with the deployer in two points A and B
%
% 1 - Contact is enabled
% 0 - Contact is disabled
% The model have to start with both equal to 1
params.isRb    = 1;
params.isRa    = 1;

% Set the accuracy, event function and output function
options = odeset('RelTol',1e-9, 'MaxStep', 0.02, 'Events', @(t,q) event(t,q,params), 'OutputFcn', @(t,q,flag) myOutputFcn(t,q,flag,params), 'Refine', 1);
% Run
[t1, q1, ~, ~, ie] = ode45(@(t,q) dqdt(t,q,params), [0,params.tk], q0, options);
for i=1:size(t1)
    [dq, Ra, Rb, zB, yA, zA, yB, Fp] = dqdt(t1(i), q1(i,:)', params);
    aux_results = [aux_results; Ra, Rb, zB, yA, zA, yB, Fp];
end
% The model supposes that 
if ie == 1 
    params.isRa = 0; 
    options = odeset('RelTol',1e-9, 'MaxStep', 0.02, 'Events', @(t,q) event(t,q,params), 'OutputFcn', @(t,q,flag) myOutputFcn(t,q,flag,params), 'Refine', 1);
    fprintf('-----\nt = %5.3f s. Contact point A is lost\n-----\n',t1(end));
else
    error('Rb = 0 before Ra = 0 \n');
end

[t2, q2, ~, ~, ie] = ode45(@(t,q) dqdt(t,q,params), [t1(end),params.tk], q1(end,:)', options);
for i=1:size(t2)
    [dq, Ra, Rb, zB, yA, zA, yB, Fp] = dqdt(t2(i), q2(i,:)', params);
    aux_results = [aux_results; Ra, Rb, zB, yA, zA, yB, Fp];
end
if ie == 2 
    params.isRb = 0; 
    options = odeset('RelTol',1e-9, 'MaxStep', 0.02, 'Events', @(t,q) event(t,q,params), 'OutputFcn', @(t,q,flag) myOutputFcn(t,q,flag,params), 'Refine', 1);
    fprintf('-----\nt = %5.3f s. Contact point B is lost\n-----\n',t2(end));
end
[t3, q3, ~, ~, ie] = ode45(@(t,q) dqdt(t,q,params), [t2(end),params.tk], q2(end,:)', options);
for i=1:size(t3)
    [dq, Ra, Rb, zB, yA, zA, yB, Fp] = dqdt(t3(i), q3(i,:)', params);
    aux_results = [aux_results; Ra, Rb, zB, yA, zA, yB, Fp];
end

% The results 
t = [t1;t2;t3];
q = [q1;q2;q3];

% Postprocessing 
% 
% Velocity
subplot(321);
plot(t,q(:,5));
xlabel('t, s');ylabel('Velocity, m/s');
% Angular velocity
subplot(322);
plot(t,[q(:,6)*180/pi,q(:,3)*180/pi]);
xlabel('t, s');ylabel('');
legend('Angular rate, deg/s','\phi, deg','Location','northwest');
% Reaction forces
subplot(323);
plot(t, aux_results(:,1:2));
xlabel('t, s');ylabel('Reaction forces, N');
legend('R_A','R_B');
subplot(324);
plot(t, aux_results(:,end));
xlabel('t, s');ylabel('Spring force, N');
legend('P');
subplot(325);
plot(t, aux_results(:,6)*1000);
xlabel('t, s');ylabel('y_B, mm');
% y_B should be equal to 0 while Rb > 0
subplot(326);
plot(t, aux_results(:,4)*1000);
xlabel('t, s');ylabel('y_A, mm'); 
% y_A should be equal to w + delta while Ra > 0
%% Animation

A = @(a) [cos(a) sin(a); -sin(a) cos(a)];
cubesat_points = [-params.zc, -params.yc; 
                  -params.zc, -params.yc+params.w;
                   params.L-params.zc, -params.yc+params.w;
                   params.L-params.zc, -params.yc;
                  -params.zc, -params.yc
                   ];
figure;
for i=1:size(t)
    cla;
    
    line([params.L params.L 0 0 params.L params.L],...
        [0.01 0 0 -params.w-params.delta -params.w-params.delta -params.w-params.delta-0.01],...
        'LineWidth',2,'Color','k');
    tpoints = (A(q(i,3))*cubesat_points')' + repmat([q(i,2),q(i,1)],5,1);
    hold on;
    plot(polyshape(tpoints(:,1),tpoints(:,2)));
    hold off;
    xlim([ 0.0, params.L*3  ]); 
    ylim([-params.w-0.1, +0.1]); 
    set(gca,'DataAspectRatio',[1 1 1]);
    
    getframe;
    pause(0.1);
end    




