function res=rside(t,y);

%глобальные переменные

global m M g l1 l2 l3 rc rc2;

%Объявляем псевдонимы для производных обобщенных координат:

dphi11=y(1);
dphi12=y(2);
dphi13=y(3);
dphi21=y(4);
dphi22=y(5);
dphi31=y(6);
dphi32=y(7);
dphi33=y(8);

%и самих обобщенных координат:

phi11=y(9);
phi12=y(10);
phi13=y(11);
phi21=y(12);
phi22=y(13);
phi31=y(14);
phi32=y(15);
phi33=y(16);

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

% Единичная матрица
E=[1,0,0;0,1,0;0,0,1];

% Вектора с
c11=[0;0;l1];
c12=[0;0;-l1];
c22=[0;0;l2];
c23=[0;0;-l2];
c33=[0;0;l3];

% Вектора b
b10=[0;0;5*l1/3];
b12=[0;0;-l1/3];
b13=b12;
b11=[0;0;2*l1/3];
b20=[0;0;l2];
b22=[0;0;0];
b23=[0;0;-l2];
b30=[0;0;l3/3];
b33=[0;0;-2*l3/3];

% Вектора d
d11=[0;0;l1];
d12=[0;0;2*l1];
d13=d12;
d21=[0;0;0];
d22=[0;0;l2];
d23=[0;0;2*l2];
d31=[0;0;0];
d32=[0;0;0];
d33=[0;0;l3];

% Вектора p
p11_1=A10*[1;0;0];
p12_1=[sin(phi13);cos(phi13);0];
p13_1=[0;0;1];
p21_2=[cos(phi22);0;sin(phi22)];
p22_2=[0;1;0];
p31_3=A32*[1;0;0];
p32_3=[sin(phi33);cos(phi33);0];
p33_3=[0;0;1];

% Вектора p в неподвижной системе
p110=A01*p11_1;
p120=A01*p12_1;
p130=A01*p13_1;
p210=A02*p21_2;
p220=A02*p22_2;
p310=A03*p31_3;
p320=A03*p32_3;
p330=A03*p33_3;

% Столбец нулей
z=[0;0;0];

% Матрица S
S=[-1,1,0;0,-1,1;0,0,-1];

% Матрица T
T=[-1,-1,-1;0,-1,-1;0,0,-1];

% Диагональные элементы матрицы K
J1=[((m*4*l1^2)/12) + m*rc*rc/4,0,0;0,((m*4*l1^2)/12) + m*rc*rc/4,0;0,0,m*rc*rc/2];
J2=[((m*4*l2^2)/12) + m*rc*rc/4,0,0;0,((m*4*l2^2)/12) + m*rc*rc/4,0;0,0,m*rc*rc/2];
J3=[((m*l3^2)/3) + m*(rc*rc+rc2*rc2)/8,0,0;0,((m*l3^2)/3) + m*(rc*rc+rc2*rc2)/8,0;0,0,m*(rc*rc+rc2*rc2)/4];

K11=J1+m*((d11'*d11)*E-d11*d11')+m*((d12'*d12)*E-d12*d12')+m*((d13'*d13)*E-d13*d13');
K22=J2+m*((d21'*d21)*E-d21*d21')+m*((d22'*d22)*E-d22*d22')+m*((d23'*d23)*E-d23*d23');
K33=J3+m*((d31'*d31)*E-d31*d31')+m*((d32'*d32)*E-d32*d32')+m*((d33'*d33)*E-d33*d33');

% не диагональные элементы K
K12=M*(b20'*(A21*d12)*E-b20*((A21*d12)'));
K13=M*(b30'*(A31*d13)*E-b30*((A31*d13)'));
K23=M*(b30'*(A32*d23)*E-b30*((A32*d23)'));
K21=M*(d12'*(A12*b20)*E-d12*((A12*b20)'));
K31=M*(d13'*(A13*b30)*E-d13*((A13*b30)'));
K32=M*(d23'*(A23*b30)*E-d23*((A23*b30)'));

% переводим каждый компонент матрицы K в неподвижную систему
K110=A01*K11*A10;
K120=A02*K12*A20;
K210=A01*K21*A10;
K130=A03*K13*A30;
K310=A01*K31*A10;
K220=A02*K22*A20;
K230=A03*K23*A30;
K320=A02*K32*A20;
K330=A03*K33*A30;

% матрица частей угловых скоростей(неполных) omeg 
omeg1 = [-cos(phi13)*sin(phi12);sin(phi12)*sin(phi13);cos(phi12)]*dphi11*dphi12+[-cos(phi12)*sin(phi13);-cos(phi12)*cos(phi13);0]*dphi11*dphi13+[cos(phi13);-sin(phi13);0]*dphi12*dphi13;
omeg2 = [     -sin(phi22)      ;          0          ;cos(phi22)]*dphi21*dphi22;
omeg3 = [-cos(phi33)*sin(phi32);sin(phi32)*sin(phi33);cos(phi32)]*dphi31*dphi32+[-cos(phi32)*sin(phi33);-cos(phi32)*cos(phi33);0]*dphi31*dphi33+[cos(phi33);-sin(phi33);0]*dphi32*dphi33;

% матрица omeg в неподвижной системе 
omeg10=A01*omeg1;
omeg20=A02*omeg2;
omeg30=A03*omeg3;
omeg=[omeg10;omeg20;omeg30];

% Относительные угловые скорости в своих системах
omega_ot1=p11_1*dphi11+p12_1*dphi12+p13_1*dphi13;
omega_ot2=p21_2*dphi21+p22_2*dphi22;
omega_ot3=p31_3*dphi31+p32_3*dphi32+p33_3*dphi33;

% Относительные угловые скорости в неподвижной системе
omega_ot10=A01*omega_ot1;
omega_ot20=A02*omega_ot2;
omega_ot30=A03*omega_ot3;

% Компоненты абсолютной угловой скорости в неподвижной системе
omega10=omega_ot10;
omega20=omega_ot10+omega_ot20;
omega30=omega_ot10+omega_ot20+omega_ot30;

% матрица omegaz (в неподвижной системе)
omegaz=[z;cross(omega10,omega20);cross(omega20,omega30)];

% матрица Fi (в неподвижной системе)
F1=[0;0;-m*g];
F2=[0;0;-m*g];
F3=[0;0;-m*g];


% компоненты b и d в неподвижной системе
b100=A01*b10;
b200=A02*b20;
b300=A03*b30;
d110=A01*d11;
d120=A01*d12;
d130=A01*d13;
d210=A02*d21;
d220=A02*d22;
d230=A02*d23;
d310=A03*d31;
d320=A03*d32;
d330=A03*d33;

% матрица Ms
Ms1=-cross(omega10,K110*omega10)-M*(cross(d120,cross(omega20,cross(omega20,b200))) + cross(d130,cross(omega30,cross(omega30,b300))))-cross(d110,F1)-cross(d120,F2)-cross(d130,F3);
Ms2=-cross(omega20,K220*omega20)-M*(cross(d230,cross(omega30,cross(omega30,b300)))+cross(b200,cross(omega10,cross(omega10,d120))))-cross(d220,F2)-cross(d230,F3);
Ms3=-cross(omega30,K330*omega30)-M*(cross(b300,cross(omega10,cross(omega10,d130))+cross(omega20,cross(omega20,d230))))-cross(d330,F3);

% блочная матрица P
P = cell(8,3);
for i=1:8
    for j=1:3
        P{i,j}=z;
    end
end;
P{1,1}=p110;
P{2,1}=p120;
P{3,1}=p130;
P{4,2}=p210;
P{5,2}=p220;
P{6,3}=p310;
P{7,3}=p320;
P{8,3}=p330;

% Блочная матрица K
K = cell(3,3);
K{1,1}=K110;
K{1,2}=K120;
K{1,3}=K130;
K{2,1}=K210;
K{2,2}=K220;
K{2,3}=K230;
K{3,1}=K310;
K{3,2}=K320;
K{3,3}=K330;

% Блочная матрица f
f = cell(3,1);
f{1,1}=omeg10+z;
f{2,1}=omeg20+cross(omega10,omega20);
f{3,1}=omeg30+cross(omega20,omega30);

% Блочная матрица MS
MS = cell(3,1);
MS{1,1}=Ms1;
MS{2,1}=Ms2;
MS{3,1}=Ms3;



% Расчет матрицы А
% (pT)
 PT = cell(8,3);
 for i=1:8
     for j=1:3
         PT{i,j} = [0;0;0];
         for k=1:3
             PT{i,j} = PT{i,j} + P{i,k}*T(k,j);
         end
     end
 end

% (pT)*K
PTK = cell(8,3);
for i=1:8
    for j=1:3
        PTK{i,j} = [0 0 0];
        for k=1:3
            PTK{i,j} = PTK{i,j} + PT{i,k}'*K{k,j};
        end
    end
end

% (pT)*K*(pT)'
PTKPT = zeros(8,8);
for i=1:8
    for j=1:8
        PTKPT(i,j) = 0; 
        for k=1:3
        PTKPT(i,j) = PTKPT(i,j)+PTK{i,k}*PT{j,k};
        end
    end
end

A = PTKPT;

% Расчет матрицы B

% T'f
Tf = cell(3,1);
for i=1:3
    Tf{i}=[0;0;0];
    for j=1:3
        Tf{i}=Tf{i}+T(j,i)*f{j,1};                                      
    end
end

% K(T'f)+MS
KTfM = cell(3,1);
for i=1:3
    KTfM{i}=z;
    for j=1:3
        KTfM{i}=KTfM{i}+K{i,j}*Tf{j};
    end
    KTfM{i}=KTfM{i}+MS{i,1};
end

% B
B = zeros(8,1);
for i=1:8
        B(i) = 0; 
        for k=1:3
        B(i) =B(i)-PT{i,k}'*KTfM{k};
        end
end

% 
ddphi=A\B;

res=[ddphi;y(1:8)];