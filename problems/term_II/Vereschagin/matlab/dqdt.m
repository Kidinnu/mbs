function dx = dqdt(t,x)
%DQDT Summary of this function goes here
%   Detailed explanation goes here

global Model;

n=Model.n;

q =x(3:3+n-1);
dq=x(n+7:n+7+n-1);

A0=cell(n,1); 
A0{1}=getA(x(Model.iq{1}(3)));

for i=2:n
    A0{i}=A0{i-1}*getA(x(Model.iq{i}));    
end

Q=cell(n,1);
F=[0;0;0];
m=0;
for i=1:n
    Q{i}=[0;0;0];    
    m=m+Model.mass{1}(1,1);
end
Q{1}=[-5;0;0];
ac=F/m;
%
dx=x;
dx(cell2mat(Model.iq)) = x(cell2mat(Model.idq));
iaux=Model.iauxq;
dx(iaux(1:2))=dx(iaux(3:4));
%
[Mk,Qk,iUk]=getMkQkiUk(x,Model,Q,A0);

w=[0;0;0];
for i=1:n
    Ck=getCk(i, x, Model, A0);
    Sk=getSk(i, x, Model, A0);
    Wp=getWprim(i, x, Model, A0);
    if i==1
       dx(Model.idq{1})=iUk{i}*Sk'*(Qk{i}-Mk{i}*(Ck*w+Wp));
       %dx(n+5:n+7)=iUk{i}*Sk'*(Qk{i}-Mk{i}*(Ck*w+Wp));
       w=inv(Mk{1})*Qk{1};
    else       
       dx(Model.idq{i})=iUk{i}*Sk'*(Qk{i}-Mk{i}*(Ck*w+Wp)); 
       % dx(n+6+i)   =iUk{i}*Sk'*(Qk{i}-Mk{i}*(Ck*w+Wp)); 
       w=Ck*w+Sk*dx(Model.idq{i})+Wp;
    end    
end
end

