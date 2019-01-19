%
% 
%
function ck = getCk(q,c,A0)
ck=[0 1;-1 0]*A0(k)*c{k,k};