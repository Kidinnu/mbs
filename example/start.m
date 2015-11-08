%Объявление глобальных переменных и констант

clear;
global m M g l1 l2 l3 rc rc2;
global yres;

%Постоянные
m=1;
M=3*m;
g=9.81;
l1=1;
l2=1;
l3=1;
rc=0.01;
rc2=0.02;

%Начальные условия

y0=[0; %dphi11
    0; %dphi12
    0; %dphi13
    0; %dphi21
    0; %dphi22
    0; %dphi31
    0; %dphi32
    0; %dphi33
    0;    %phi11
    0;    %phi12
    0;    %phi13
    0; %phi21 - Отклонение тела 2 в плоскости y0z на угол
    0;    %phi22
    pi/6;    %phi31 
    pi/3;    %phi32
    0];   %phi33

%Начальное и конечное время интегрирования

tspan=[0 10];

%Вызываем функцию численного интегрирования и передаем ей ссылку на имя функции правой части диф. уравнений.

[t,yres]=ode113(@rside,tspan,y0);

%Строим необходимые графики, например угла phi11, записываем результаты
%численного интегрирования в текстовые файлы для последующего чтения при
%создании анимации

plot(t,yres(:,16));
csvwrite('phi11.txt', yres(:,9));
csvwrite('phi12.txt', yres(:,10));
csvwrite('phi13.txt', yres(:,11));
csvwrite('phi21.txt', yres(:,12));
csvwrite('phi22.txt', yres(:,13));
csvwrite('phi31.txt', yres(:,14));
csvwrite('phi32.txt', yres(:,15));
csvwrite('phi33.txt', yres(:,16));