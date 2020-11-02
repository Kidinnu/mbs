% Landing gear deployment model
% Copyright © 2020 Vadim Yudintsev. All rights reserved.

p.k     = 1.4;
p.d     = 0.3;
p.g     = 9.807;
p.Fp    = p.d*p.d*pi/4;
p.K     = sqrt(2*p.g*p.k/(p.k-1));
p.fe    = 0.0005*3;
p.pm    = 1e4;
p.R     = 29.27;
p.Tm    = 293;

p.m2m1  = 3.0; p.nx    = 2;
p.m1    = 100; p.l     = 6.8; p.eta = 5.0/6.0;

clc;

% Этап раскрытия опоры с закрытым дренажным клапаном
p.stage = 1; tI = 1.3;
[t1,q1] = ode113(@(t,q) dpdt(t,q,p),[0 tI],[0.0; 0.05; 10100*2.7],odeset('RelTol',1e-7));

% Этап раскрытия опоры с открытым дренажным клапаном
p.stage = 2;
[t,q] = ode113(@(t,q) dpdt(t,q,p),[t1(end) 3],q1(end,:)',odeset('RelTol',1e-7));

t = [t1;t(2:end)];
q = [q1;q(2:end,:)];

% Давление в демпфере
figure('Position',[100 100 600 350]);
axes('Position',[0.1 0.12 0.8 0.85]);
plot(t,q(:,3)/p.pm,'LineWidth',2);
xlabel('t, c');
ylabel('p_1, атм');
grid on;
set(gca,'FontSize',12);

% Угол поворота подкоса и угловая скорость
figure('Position',[100 100 600 350]);
axes('Position',[0.1 0.12 0.8 0.85]);
plot(t,q(:,1)*180/pi,'-','LineWidth',2);
ylabel('\phi, ...^o');
ylim([0,140]);
yyaxis right;
plot(t,q(:,2)*180/pi,'--','LineWidth',2);
ylabel('\omega ...^o/c');
ylim([0,140]);
yyaxis left;
xlabel('t, c');
grid on;
fprintf('de = %4.1f mm\n',sqrt(4*p.fe/pi)*1000);
fprintf('we = %4.1f deg/s\n',max(q(:,2))*180/pi);
set(gca,'FontSize',12);
legend('Угол поворода', 'Угловая скорость');

%
% Функция правых частей дифференциальных уравнений
%
function [dq, Fd] = dpdt(t,q,p)
    phi  = q(1);
    dphi = q(2);
    p1   = q(3);

    x    = sqrt(p.l^2+p.l^2*p.eta^2-2*p.l^2*p.eta*cos(phi));
    dx   = p.l^2*p.eta*sin(phi)*dphi/x;
    K    = sqrt(2*p.g*p.k/(p.k-1));

    if p.stage == 1
        dp1  = - p.k*p1*dx/(p.x0+x);
    else
        dp1  = p.k*p.fe*K*p.pm*sqrt(p.R*p.Tm)/(p.Fp*(p.x0+x))*phi_gas(p1/p.pm,p.k)...
                - p.k*p1*dx/(p.x0+x);
    end

    Fd   =-(p.pm-p1)*p.Fp*p.g;

    philim = 120;
    if phi > philim*pi/180
        Fd = Fd - (phi-philim*pi/180)*1e7 - dphi*4e6;
    end

    arm = p.l*p.eta*sin(phi)*sqrt(1/(1+p.eta^2-2*p.eta*cos(phi)));

    d2phi = 3*p.g*(1+p.m2m1)*p.m1*p.nx*sin(phi) +...
        6*Fd*p.eta*sin(phi)*sqrt(1/(1+p.eta^2-2*p.eta*cos(phi)));

    d2phi = d2phi/(2*p.m2m1*p.l*p.m1);

    dq  = [dphi; d2phi; dp1];
end

%
% Отановка интегрирования при развороте на 120 градусов
%
function [value,isterminal,direction]  = dpdt_ev(t, q)
    value = 120 - q(1)*180/pi;
    isterminal = 1;
    direction = -1;
end

%
% Функция расхода
%
function res = phi_gas(s, k)
    if s>= 0.528 && s<=1 
        res = sqrt(s^(2/k)-s^((k+1)/k));
    else
        res = 0;
    end

    if s < 0.528
        res = 0.2588;
    end
end

