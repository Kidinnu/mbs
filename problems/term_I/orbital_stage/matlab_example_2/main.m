clc;

% Масса
params.m = 3000.0;
% Матрица тензора инерции
params.J = [5000.0,     0,     0;
            0     , 10000,     0;
            0     , 0    , 10000];
% Матрциа обратная матрице тензора инерции
params.invJ = inv(params.J);

% Функция силы от времени
% Текущее значение силы интерполируется по таблице значений
params.sopForceTable = @(t) interp1([0 20 30 50],[1000 500 100 0],t);

% Точка приложения силы в ССК                     
params.pFsop = [3.0; 
                2.0*sin(45*pi/180.0); 
                2.0*cos(45*pi/180.0)];                    
% Единичный вектор направления силы в ССК
params.nFsop = [-cos(10*pi/180.0); 
                +sin(10*pi/180.0)*cos(40*pi/180.0);
                -sin(10*pi/180.0)*cos(40*pi/180.0)];                    
%
% Начальные условия
%

% Положение центра масс
r0 = [0;0;0];
% Скорость центра масс
v0 = [0;0;0];
% Начальная ориентация (углы Брайнта 1-2-3)
a0 = [0;0;0];
% Начальная угловая скорость
w0 = [0;0;0];

% Начальный вектор состояния
q0 = [r0;a0;v0;w0];
% Интервал интегрирования
tspan = 0:0.2:50;
% Относительная погрешность
option = odeset('RelTol',1e-6);

% Интегрирование
[t, q] = ode45(@(t,q) dqdt_orbital_stage(t,q,params), tspan, q0, option);


%% Обработка результатов
% Г рафики
% Координаты центра масс от времени
subplot(2,2,1);
plot(t,q(:,1:3));
xlabel('t, c'); ylabel('x, y, z, m');
legend('x','y','z');
grid on;
% Скорость центра масс от времени
subplot(2,2,2);
plot(t,q(:,7:9));
xlabel('t, c'); ylabel('Vx, Vy, Vz, m/c');
legend('Vx','Vy','Vz');
grid on;
% Углы 
subplot(2,2,3);
plot(t,q(:,4:6));
xlabel('t, c'); ylabel('a_1, a_2, a_3, градусы');
legend('a_1','a_2','a_3');
grid on;
% Угловая скорость 
subplot(2,2,4);
plot(t,q(:,10:12)*180/pi);
xlabel('t, c'); ylabel('\omega_1, \omega_2, \omega_3, градус/c');
legend('\omega_1','\omega_2','\omega_3');
grid on
% Вектор кинетического момента 
% три колонки Lx Ly Lz
% L = J*w
L = (params.J*q(:,10:12)')';
theta = acos(L(:,1)./sqrt(sum(L.^2,2)))*180/pi;

%plot(t,theta)


%%

figure;
axis([-10 10 -10 10 -10 10]*0.8);
hold on;

box;
axis vis3d;

for i=1:size(t,1)
    cla;
    r = q(i,1:3)';
    A = Axyz(q(i,4),q(i,5),q(i,6));
    
    draw_upper_stage([4; 0; 0], [0; 0; 0], A);
    
    shading flat;
    lighting gouraud;
    light('Position',[0 0 10]);
    
    getframe;
end


