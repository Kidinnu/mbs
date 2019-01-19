%
% 
%
function sk = getSk(k, x, Model, A0)    
    if k==1
        sk=eye(3);
    else
        sk=[0;0;1];
        sk(1:2)=[0 1;-1 0]*A0{k}*Model.c{k,k};
    end
