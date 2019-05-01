%
% 
%
function ck = getCk(k, x, Model, A0)
ck=eye(3,3);
if k==1
    ck=eye(3,3)*0;
else
    ck(1:2,3)=-[0 1;-1 0]*(A0{k-1}*Model.c{k-1,k}-A0{k}*Model.c{k,k});
end