%
%   Метод отдельных тел
%   Тело 1 не связано с "нулевым" телом  
% 
global Model;
Model=struct;
% Количество тел
n=4;
Model.n=n;
% Шарнирные векторы
Model.c=cell(n,n);
c=cell(n,n);
L=1;
for i=1:n
    for j=1:n
        c{i,j}=[0;0];        
    end
    c{i,i}=[-L;0]/n;
    if i~=n
        c{i,i+1}=[+L;0]/n;
    end
end
% Первый фиктивный шарнир имеет три степени свободы и расположен в
% центре масс тела 1
c{1,1}=[0;0];
Model.c=c;
% Количество дополнительных дифференциальных уравнений 
Model.n0=4;
% Число степеней свободы в шарнирах
Model.na=ones(n,1);
Model.na(1)=3;
% Матрицы масс
Model.mass=cell(n,1);
mass=cell(n,1);
for i=1:n
    mass{i}=[1 0 0; 0 1 0; 0 0 4/n/n/12]/n;
end
Model.mass=mass;
%
Model=init_MBSModel(Model);

%% Начальные условия
q0=zeros((n+2+2)*2,1);
q0(1)=L/2;
q0(4)=pi/2.0;
[t,q]=ode113(@dqdt,[0 20],q0);
%% Графики

plot(t,q(:,Model.iq{1}(3)));

%% Анимация

Lrod=c{1,2}(1,1)*2;

M = moviein(size(t,1));
axis([-n*Lrod n*Lrod -2*n*Lrod n*Lrod]*1);

hold on;

p=cell(n,1);

maxY=0;

for i=1:1:size(q,1) 
    cla; 
    a=q(i,3);
    rc1=[q(i,1);q(i,2);0];
    p0=[q(i,1);q(i,2);0]-Lrod/2*[cos(a);sin(a);0];        
    p{1}=p0+[cos(a);sin(a);0]*Lrod;
    rod=[p0';p{1}'];
    line(rod(:,1),rod(:,2),rod(:,3),'LineWidth',2);    
    for j=2:n
        a=a+q(i,j+2);
        p{j}=p{j-1}+[cos(a);sin(a);0]*Lrod;
        rod=[p{j-1}';p{j}'];
        line(rod(:,1),rod(:,2),rod(:,3),'LineWidth',2);        
        if (j>n-2)
            line(rod(:,1),rod(:,2),rod(:,3),'LineWidth',2);        
        end
        if (j==n && maxY<p{n}(2))
            maxY=p{n}(2);
        end
    end        
    a=getframe;
    %M(:,i) = getframe;    
    %imwrite(M(:,i),'img.png','PNG','Resolution',100);
    %mo= getframe;
    %mov=addframe(mov,mo);               
end
%%
% Energy

qabs=q*0;
v=zeros(size(q,1),2);
h=zeros(size(q,1),2);

T=0;
P=0;

for i=1:n
    v=zeros(size(q,1),2);
    h=zeros(size(q,1),2);
    for j=1:i
       qabs(:,i)=qabs(:,i)+q(:,j); 
       qabs(:,n+i)=qabs(:,n+i)+q(:,n+j);       
    end
    for j=1:i
       if j==i 
           b=-c{j,j};
       else
           b=-c{j,j}+c{j,j+1};
       end
       h=h+reshape(getA(qabs(:,j))*b,size(q,1),2);
       v=v-reshape(getA(qabs(:,j))*[0 1;-1 0]*b,size(q,1),2).*[qabs(:,n+j) qabs(:,n+j)];      
    end    
    T=T+mass{i}(1,1)*(v(:,1).^2+v(:,2).^2)*0.5+mass{i}(3,3)*qabs(:,n+i).^2*0.5;
    %P=P+Model.mass{i}(1,1)*9.81*h(:,2);    
    P=P+5*h(:,1);    
end

figure; 
P=-5*(q(:,1)-q(1,1));
E0=T+P;
plot(t,T+P-E0(1));
title('Energy error');




