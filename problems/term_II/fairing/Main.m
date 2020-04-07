P.r0 = [-4; 1];         %Location of CS-origin fixed to hinge wrt COM of the fairing, meters

P.rp = [-3.5; -0.8];    %Location of pusher-spring hinge (point P) on fairing wrt COM of fairing, meters

P.ra = [0; -1.8];       %Location of pusher-spring hinge (point A) wrt hinge-fixed CS, meters

P.ax = 9.81;            %Acceleration of launcher, m/s^2

P.P0 = 30000;           %Initial pusher-spring force, N

P.Pk = 10000;            %Final pusher-spring force, N

P.dh = 0.4;             %Pusher-spring stroke length, meters

P.m = 1000;             %Mass of head fairing, Kg

P.J = 10000;            %Moment of inertia of head fairirng, Kg*m^2

%Initial state vector, [x0; y0; phi0; dx0; dy0; dphi0]
q0 = [4; -1; 0; 0; 0; 0]; 

P.Stage = 1;
opt = odeset('Events', @PhiEvent);
[t1,q1] = ode45(@(t, q) dqdt_FairSep(t, q, P), [0 2], q0, opt);
P.Stage = 2;
[t2,q2] = ode45(@(t, q) dqdt_FairSep(t, q, P), [t1(end) 3], q1(end, :)');

clc;

t = [t1(1:end - 1); t2];
q = [q1(1:end - 1, :); q2];  


%%
contour_x = [0 7 10 0 0];
contour_y = [0 0 -2 -2 0];
contour_x  = contour_x + P.r0(1);
contour_y  = contour_y + P.r0(2);

for i=1:size(t)
    phi = q(i,3);
    ci = [cos(phi) -sin(phi); sin(phi) cos(phi)]*[contour_x;contour_y];    
    plot(ci(2,:)+q(i,2),ci(1,:)+q(i,1)); 
    hold on;
    plot([0 0 -2],[-2 0 0]);
    hold off;
    xlim([-2 11]);
    ylim([-2 11]);
    daspect([1 1 1]);
    getframe;
end



%%


%Position wrt hinge-fixed CS
subplot(2, 1, 1);
plot(t, q(:, 4));
xlabel('time, s');  ylabel('Vx, meters/sec');

%Angular position wrt hinge-fixed CS
subplot(2, 1, 2);
plot(t, q(:, 3)*180/pi);
xlabel('time, s');  ylabel('Phi, degrees');

%Cartesian vector of fairing (COM)
figure;
plot(q(:, 1), q(:, 2));
xlabel('x0, meters');    ylabel('y0, meters');








