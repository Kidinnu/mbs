%
% Интеграл энергии
%
global m M g rc l1 l2 l3 rezultat;
 
% Моменты инерции
J1=[(m*4*l1^2)/12,0,0;0,(m*4*l1^2)/12,0;0,0,(m*rc^2)/2];
J2=[(m*4*l2^2)/12,0,0;0,(m*4*l2^2)/12,0;0,0,(m*rc^2)/2];
J3=[(m*4*l3^2)/12,0,0;0,(m*4*l3^2)/12,0;0,0,(m*rc^2)/2];
 
N=length(t);
W=[]; D=[];
 
figure;
mov = moviein(N);
axis([-4 4 -8 2]);
hold on;
 
for j=1:N
    
cla;
    
dphi11=rezultat(j,1);
dphi12=rezultat(j,2);
dphi21=rezultat(j,3);
dphi31=rezultat(j,4);
dphi32=rezultat(j,5);

phi11=rezultat(j,6);
phi12=rezultat(j,7);
phi21=rezultat(j,8);
phi31=rezultat(j,9);
phi32=rezultat(j,10);


 
 
% матрица перехода из неподвижной системы в первую
A10=[cos(phi11)*cos(phi12), sin(phi12),    -sin(phi11)*cos(phi12);
    -cos(phi11)*sin(phi12), cos(phi12),    sin(phi11)*sin(phi12);
	sin(phi11),                 0,                    cos(phi11)];

% матрица перехода из первой системы во вторую
A21=[ cos(phi21),    0,      sin(phi21) ;
     0,1,0 ;
      -sin(phi21),0,cos(phi21)];
 
% матрица перехода из второй системы в третью
A32=[cos(phi32), sin(phi32)*sin(phi31), -sin(phi32)*cos(phi31);
              0,            cos(phi31),      sin(phi31);
     sin(phi32),-cos(phi32)*sin(phi31),cos(phi32)*cos(phi31)];
 
% матрица перехода из первой системы в третью
A31=A32*A21;
 
% матрица перехода из третьей системы в первую
A13=A31';
 
% матрица перехода из первой системы в неподвижную
A01=A10';
 
% матрица перехода из второй системы в первую
A12=A21';
 
% матрица перехода из второй системы в неподвижную
A02=A01*A12;
 
% матрица перехода из неподвижной системы во вторую
A20=A02';
 
% матрица перехода из третьей системы во вторую
A23=A32';
 
% матрица перехода из третьей системы в неподвижную
A03=A01*A12*A23;
 
% матрица перехода из неподвижной системы в третью
A30=A03';
 
% векторы p
p11_1=[sin(phi12);cos(phi12);0];
p12_1=[0;0;1];
p21_2=[1;0;0];
p31_3=[cos(phi32);0;sin(phi32)];
p32_3=[0;1;0];
% Относительные угловые скорости в своих системах
omega_ot1=p11_1*dphi11+p12_1*dphi12;
omega_ot2=p21_2*dphi21;
omega_ot3=p31_3*dphi31+p32_3*dphi32;
 
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
%V20=cross(omega10,-d120)+cross(omega_ot20,-d220);
V20=cross(omega10,-d120)+cross(omega20,-d220);
V2_2=V20'*V20;
%V30=cross(omega10,-d120)+cross(omega_ot20,-d230)-cross(omega_ot30,-d330);
V30=cross(omega10,-d120)+cross(omega20,-d230)+cross(omega30,-d330);
V3_2=V30'*V30;
 
% Тензоры инерции в неподвижной системе
J10=A01*J1*A10;
J20=A02*J2*A20;
J30=A03*J3*A30;
 
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
 
% 
%p0=[0;0;0];
%pA=p0+[2*l1*sin(phi11);-2*l1*cos(phi11);0];
%pB=pA+[2*l2*sin(phi21+phi11);-2*l2*cos(phi21+phi11);0];
%pC=pB+[2*l3*sin(phi11+phi21+phi31);-2*l3*cos(phi11+phi21+phi31);0];
 
p0=[0;0;0];
pA=-d120;
pB=-d120-d230;
pC=-d120-d230-d330*2;
 
rod1=[p0';pA'];
rod2=[pA';pB'];
rod3=[pB';pC'];
 
line(rod1(:,2),rod1(:,3),rod1(:,1),'LineWidth',4,'Color','Green'  );
line(rod2(:,2),rod2(:,3),rod2(:,1),'LineWidth',4,'Color','Black');
line(rod3(:,2),rod3(:,3),rod3(:,1),'LineWidth',4,'Color','Blue' );
 
mov(:,j) = getframe;
 
end;

figure;
%plot(t,[W;D;W+D]);
plot(t,W+D);
