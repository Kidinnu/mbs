function [Mk,Qk,iUk] = getMkQkiUk(q,dq,masses,c,Q,A0)
%GETMK Summary of this function goes here
%   Detailed explanation goes here

n=size(q,1);

Mk=cell(n,1);Qk=cell(n,1);iUk=cell(n,1);

Mk{n}=masses{n};
Qk{n}=Q{n};

for k=n:-1:2
   Ck=getCk(k, q, dq, c, A0);
   Sk=getSk(k, q, dq, c, A0);   
   Wprim=getWprim(k, q, dq, c, A0);
   U=Sk'*Mk{k}*Sk;
   iUk{k}=inv(U);
    
   Mk{k-1}=masses{k-1}+Ck'*Mk{k}*Ck-Ck'*Mk{k}*Sk*iUk{k}*Sk'*Mk{k}*Ck;
   Qk{k-1}=Q{k-1}+Ck'*(Mk{k}*(Sk*iUk{k}*Sk'*(Qk{k}-Mk{k}*Wprim)+Wprim)-Qk{k});
end

