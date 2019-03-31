function dq = dqdt(t,q,p)





pT = zeros(3,p.N,p.n);

for i=1:p.n
   for j=1:p.n
      pT(:,i,j) = p.p{i}()*p.T(i,j);       
   end
    
end


end

