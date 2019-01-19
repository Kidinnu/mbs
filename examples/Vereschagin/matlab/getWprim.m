function  wprim = getWprim(k, x, Model, A0)

dphi = cell2mat(Model.idq);
dphi = x(dphi(3:size(dphi,1)));

if k==1       
   wprim=[0;0;0];
else
   wkp=sum(dphi(1:k-1)); 
   wprim=[(dphi(k)+wkp)^2*A0{k}*Model.c{k,k}-wkp*wkp*A0{k-1}*Model.c{k-1,k};0]; 
end


end

