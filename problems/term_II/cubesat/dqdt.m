function [dq, Ra, Rb, zB, yA, zA, yB Fp] = dqdt(t, q, params)
% ode function
%
w = params.w;
delta = params.delta;
L  = params.L;
zc = params.zc;
yc = params.yc;
m  = params.mass;
J  = params.J;

y   = q(1);
z   = q(2);
phi = q(3);

vy  = q(4);
vz  = q(5);
dphi= q(6);

% Сила толкателя зависит от перемещения спутника в внутри ТПК
% Будем считать, что это перемещение точки В от начального положения,
% которое равно нулю ( zB(0) = 0 ) -- это должно обеспечиваться начальными 
% условиями.
zB = z - (yc-w)*sin(phi)-zc*cos(phi);

% Проверяем, может быть толкатель уже закончил работу?
if zB <= params.hp  
    % нет, вычисляем силу пружины
    Fp  = params.Fp0 - params.cp*zB;
else
    % сила равна нулю, движемся "по инерции"
    Fp  = 0;
end

% Проекции активных сил (без реакций) на оси неподвижной системы (связанной с ТПК)
Fz = +Fp;
Fy =   0;

% Собираем матрицы А, В

%   А  *   X    =   B                     
%
% [m 0 0 k1 k2] [  y'']   [ b1 ]
% [0 m 0 k3 k4] [  z'']   [ b2 ]
% [0 0 J k5 k6]*[phi''] = [ b3 ]
% [0 0 0  0  0] [  Ra ]   [ b4 ]
% [0 0 0  0  0] [  Rb ]   [ b5 ]

BAp = (L - z + yc*sin(phi))/cos(phi) + zc;

k1 =+1;
k2 =-cos(phi);

k3 = 0;
k4 =-sin(phi);

k5 = zc*cos(phi) + (yc-w)*sin(phi);
k6 = (BAp-zc);

b1 = Fy;
b2 = Fz;
b3 = Fp*(zc*sin(phi)+(params.w*0.5-params.yc)*sin(phi));

% Q11*y'' + Q12*z''+ Q13 * phi'' = B1
% B
if params.isRb
    Q11 = 1;
    Q12 = 0;
    Q13 = sin(phi)*(yc-w)+zc*cos(phi);
    b4  = (zc*sin(phi) + (w-yc)*cos(phi))*dphi*dphi;
else
    Q11 = 0;
    Q12 = 0;
    Q13 = 0;
    b4  = 0;
end
    
% Q21*y'' + Q22*z''+ Q23 * phi'' = B2 
% A
if params.isRa
    Q21 =  0;
    Q22 = -tan(phi);
    Q23 = zc*cos(phi)+(yc-w)*sin(phi) + (L+yc*sin(phi)-z)*(sec(phi)^2);
    b5  = -(-16*vz-sec(phi)*dphi*(3*(w-5*yc)+4*w*cos(2*phi)+(w-yc)*cos(4*phi)+8*(zc*(cos(phi)^3)-2*L)*sin(phi)+16*sin(phi)*z))*(dphi*(sec(phi)^2))/8.0;
else    
    Q21 = 0;
    Q22 = 0;
    Q23 = 0;
    b5  = 0;
end

% 
A = [ m   0   0    k1   k2;
      0   m   0    k3   k4;
      0   0   J    k5   k6;
    Q11 Q12 Q13    1-params.isRb    0;
    Q21 Q22 Q23    0    1-params.isRa];
 
B = [b1;b2;b3;b4;b5];
%
% Решение системы (1) относительно ускорений и реакций
%
X = A\B;
% Записываем результат 
% ускорения центра масс (проекции)

% Acceleration of the center of mass
ay    = X(1);
az    = X(2);

% Angular acceleration
d2phi = X(3);

% За реакциями надо следить, чтобы "выключить" связь 
% Вернем их из функции последними двумя результатами
% см. первую строчку определения этой файл-функции
Ra    = X(5);
Rb    = X(4);

yB = y - (params.yc-params.w)*cos(phi) + params.zc*sin(phi);
yA = -params.w*cos(phi)-BAp*sin(phi);
zA = zB + params.w*sin(phi)+BAp*cos(phi);

% заполняем вектор результата
dq = [vy;vz;dphi;ay;az;d2phi];

end

