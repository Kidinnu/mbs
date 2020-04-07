clc; clear all;
% Положение шарнирной точки отчносительно центра масс створки с ССК
p.c11 = [-4.0;  1.0; 0.0];
% Положение точки закрепления толкателя на створке относительно центра масс
% створки в ССК
p.rp1 = [-3.5; -0.8; 0.0];  
% Положение точки закрепления толкателя на РН относительно БСК
p.rp0 = [ 0.0;  0.2; 0.0];  
% Ускорение свободного падения
p.g   = 9.807;            
% Начальное усиле толкателей
p.P0  = 30000;            
% Конечное усилие толкателей
p.Pk  = 10000;            
% Ход толкателя
p.h   = 0.4;   
% Начальная длина толкателя
p.d0  = 0.5;
% Масса створки
p.m   = 1000;  
% Тензор инерции створки
p.J   = diag([10000,10000,10000]);            
%
p.phimax   = 60;            

% Начальные условия
q0 = [4.0; 1.0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]; 
% Первый этап
p.Stage = 1;
opt = odeset('Events', @(t,q) PhiEvent(t,q,p));
[t1,q1] = ode45(@(t, q) fairing3d(t, q, p), 0:0.02:2, q0, opt);
% Второй этап
p.Stage = 2;
[t2,q2] = ode45(@(t, q) fairing3d(t, q, p), [t1(end):0.02:3], q1(end, :)');
% Склейка результатов
t = [t1(1:end - 1); t2];
q = [q1(1:end - 1, :); q2];  

%% Анимация
contour_x = [0 0 7 10  0];
contour_y = [0 2 2 0   0];
contour_x  = contour_x - 4;
contour_y  = contour_y - 1;

for i=1:size(t)
    clf;
    phi = q(i,6);    
    ci = [cos(phi) -sin(phi); sin(phi) cos(phi)]*[contour_x;contour_y];    
    patch(ci(2,:)+q(i,2),ci(1,:)+q(i,1),[0.5 0.6 0.5]); 
    hold on;
    plot([-2 2 2],[0 0 -4]);
    hold off;
    xlim([-2 15]);
    ylim([-4 11]);
    daspect([1 1 1]);
    getframe;
end


%%

subplot(2, 1, 1);
plot(t, q(:, 4));
xlabel('time, s');  ylabel('Vx, meters/sec');


subplot(2, 1, 2);
plot(t, q(:, 3)*180/pi);
xlabel('time, s');  ylabel('Phi, degrees');

figure;
plot(q(:, 1), q(:, 2));
xlabel('x0, meters');    ylabel('y0, meters');








