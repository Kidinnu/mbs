function [Mk,Qk,iUk] = getMkQkiUk(x,Model,Q,A0)
%GETMK Summary of this function goes here
%   Detailed explanation goes here

n=Model.n;

Mk=cell(n,1);Qk=cell(n,1);iUk=cell(n,1);

Mk{n}=Model.mass{n};
Qk{n}=Q{n};

for k=n:-1:1
   Ck=getCk(k, x, Model, A0);
   Sk=getSk(k, x, Model, A0);   
   Wprim=getWprim(k, x, Model, A0);
   U=Sk'*Mk{k}*Sk;
   iUk{k}=inv(U);   
   if k>1 
        Mk{k-1}=Model.mass{k-1}+Ck'*Mk{k}*Ck-Ck'*Mk{k}*Sk*iUk{k}*Sk'*Mk{k}*Ck;
        Qk{k-1}=Q{k-1}-Ck'*(Mk{k}*(Sk*iUk{k}*Sk'*(Qk{k}-Mk{k}*Wprim)+Wprim)-Qk{k});
   end
end

