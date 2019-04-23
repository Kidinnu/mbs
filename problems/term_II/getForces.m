function F = getForces(t,q,p)
% in 0 frame
F = zeros(3,p.n);

for i=1:p.n
   F(:,i) = [0;0;-p.m(i)*9.807]; 
end

end

