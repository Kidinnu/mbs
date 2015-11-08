% Интеграл энергии
%
global m M g rc rc2 l1 l2 l3 yres;

% Моменты инерции
J1=[((m*4*l1^2)/12) + m*rc*rc/4,0,0;0,((m*4*l1^2)/12) + m*rc*rc/4,0;0,0,m*rc*rc/2];
J2=[((m*4*l2^2)/12) + m*rc*rc/4,0,0;0,((m*4*l2^2)/12) + m*rc*rc/4,0;0,0,m*rc*rc/2];
J3=[((m*l3^2)/3) + m*(rc*rc+rc2*rc2)/8,0,0;0,((m*l3^2)/3) + m*(rc*rc+rc2*rc2)/8,0;0,0,m*(rc*rc+rc2*rc2)/4];

N=length(t);

for j=1:N
    
  
  
dphi11=yres(j,1);
dphi12=yres(j,2);
dphi13=yres(j,3);
dphi21=yres(j,4);
dphi22=yres(j,5);
dphi31=yres(j,6);
dphi32=yres(j,7);
dphi33=yres(j,8);

phi11=yres(j,9);
phi12=yres(j,10);
phi13=yres(j,11);
phi21=yres(j,12);
phi22=yres(j,13);
phi31=yres(j,14);
phi32=yres(j,15);
phi33=yres(j,16);

%Вычисление матриц преобразования координат:

% матрица перехода из неподвижной системы в первую
A10=[cos(phi12)*cos(phi13), cos(phi11)*sin(phi13)+sin(phi11)*sin(phi12)*cos(phi13), sin(phi11)*sin(phi13)-cos(phi11)*sin(phi12)*cos(phi13);
    -cos(phi12)*sin(phi13), cos(phi11)*cos(phi13)-sin(phi11)*sin(phi12)*sin(phi13), sin(phi11)*cos(phi13)+cos(phi11)*sin(phi12)*sin(phi13);
    sin(phi12),             -sin(phi11)*cos(phi12),                                 cos(phi11)*cos(phi12)];

% матрица перехода из первой системы во вторую
A21=[cos(phi22), sin(phi21)*sin(phi22) , -cos(phi21)*sin(phi22);
        0      ,      cos(phi21)       ,        sin(phi21)     ;
     sin(phi22), -sin(phi21)*cos(phi22), cos(phi21)*cos(phi22)];

% матрица перехода из второй системы в третью
A32=[cos(phi32)*cos(phi33), cos(phi31)*sin(phi33)+sin(phi31)*sin(phi32)*cos(phi33), sin(phi31)*sin(phi33)-cos(phi31)*sin(phi32)*cos(phi33);
    -cos(phi32)*sin(phi33), cos(phi31)*cos(phi33)-sin(phi31)*sin(phi32)*sin(phi33), sin(phi31)*cos(phi33)+cos(phi31)*sin(phi32)*sin(phi33);
    sin(phi32),             -sin(phi31)*cos(phi32),                                 cos(phi31)*cos(phi32)];

% Дополнительные матрицы перехода
A31=A32*A21;
A13=A31';
A01=A10';
A12=A21';
A02=A01*A12;
A20=A02';
A23=A32';
A03=A01*A12*A23;
A30=A03';

% Вектора p
p11_1=A10*[1;0;0];
p12_1=[sin(phi13);cos(phi13);0];
p13_1=[0;0;1];
p21_2=[cos(phi22);0;sin(phi22)];
p22_2=[0;1;0];
p31_3=A32*[1;0;0];
p32_3=[sin(phi33);cos(phi33);0];
p33_3=[0;0;1];

% Относительные угловые скорости в своих системах
omega_ot1=p11_1*dphi11+p12_1*dphi12+p13_1*dphi13;
omega_ot2=p21_2*dphi21+p22_2*dphi22;
omega_ot3=p31_3*dphi31+p32_3*dphi32+p33_3*dphi33;

% Относительные угловые скорости в неподвижной системе
omega_ot10=A01*omega_ot1;
omega_ot20=A02*omega_ot2;
omega_ot30=A03*omega_ot3;

% векторы d
d11=[0;0;l1];
d12=[0;0;2*l1];
d13=d12;
d21=[0;0;0];
d22=[0;0;l2];
d23=[0;0;2*l2];
d31=[0;0;0];
d32=[0;0;0];
d33=[0;0;l3];

% компоненты d в неподвижной системе
d110=A01*d11;
d120=A01*d12;
d130=A01*d13;
d210=A02*d21;
d220=A02*d22;
d230=A02*d23;
d310=A03*d31;
d320=A03*d32;
d330=A03*d33;

% Компоненты абсолютной угловой скорости в неподвижной системе
omega10=omega_ot10;
omega20=omega_ot10+omega_ot20;
omega30=omega_ot10+omega_ot20+omega_ot30;

% Компоненты скоростей
V10=cross(omega10,-d110);
V1_2=V10'*V10;
V20=cross(omega10,-d120-d220)+cross(omega_ot20,-d220);
V2_2=V20'*V20;
V30=cross(omega10,-d120-d230-d330)+cross(omega_ot20,-d230-d330)+cross(omega_ot30,-d330);
V3_2=V30'*V30;

% Тензоры инерции в неподвижной системе
J10=A01*J1*A01';
J20=A02*J2*A02';
J30=A03*J3*A03';

% Кинетическая энергия W
W1=((omega10'*J10*omega10)/2)+m*V1_2/2;
W2=((omega20'*J20*omega20)/2)+m*V2_2/2;
W3=((omega30'*J30*omega30)/2)+m*V3_2/2;

W(j)=W1+W2+W3;

% Потенциальная энергия D. Нулевой уровень проходит через первый шарнир
h1=-d110;
h2=-d120-d220;
h3=-d120-d230-d330;

D1=m*g*h1(3);
D2=m*g*h2(3);
D3=m*g*h3(3);

D(j)=D1+D2+D3;

end;
%%
figure;
plot(t,[W;D;W+D]);

