function dq = dqdt_orbital_stage(t,q,params)

r = q(1:3);
a = q(4:6);
v = q(7:9);
w = q(10:12); 

% Матрица преобразования координат из ССК в ИСК
A = Axyz(a(1),a(2),a(3));

F = params.sopForceTable(t);

%
% Вектор силы в ССК
% Значение * единичный вектор
%
Fsop_b = F*params.nFsop;
%
% Вектор момента этой силы относительно ССК
%
Msop_b = cross(params.pFsop, Fsop_b);

%
% Вектор силы в ИСК
%
Fsop_0 = A*Fsop_b;

%
% Ускорение центра масс
%
dv = Fsop_0/params.m;

%
% Угловое ускорение 
%
dw = params.invJ*(Msop_b - cross(w,params.J*w));

% Кинематические уравнения для углов Брайнта
% Производные углов
%
da = kinematicEq123(a, w);

% Собираем ответ 
% Значение правой части системы ДУ
dq = [v; da; dv; dw];

end

