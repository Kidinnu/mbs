clear;
clear global;
global m M g rc l1 l2 l3;
global rezultat;
 
% const
m=1;
M=3;
g=9.81;
rc=0.02;
l1=0.75;
l2=1;
l3=0.5;
 
% начальные условия
 
y0=[0; 
    0; 
0;
0;
    0; 
    0; 
    pi/3;    
    0;    
    pi/10;   
    pi/10];   
% начальное и конечное время
tspan=[0 5];
 
% функция интегрирования
[t,rezultat]=ode113(@ur,tspan,y0);
 
% строим phi11(t)
plot(t,rezultat(:,6),'LineWidth',2,'Color','Red');
subplot(2,3,1),plot(t,rezultat(:,6),'-k')