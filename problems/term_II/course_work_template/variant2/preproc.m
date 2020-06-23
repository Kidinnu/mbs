
p.N  = sum(p.na);
p.iq = [(cumsum(p.na)-p.na+1) p.na];
p.T  = inv(p.S);  

% d_{ji} column vectors in j-CS
p.d  = zeros(3,p.n,p.n);
for i = 1:p.n
    for j = 1:p.n
        for a = 1:p.n
            p.d(:,j,i) = p.d(:,j,i) + p.T(a,i)*p.S(j,a)*p.C(:,j,a);
        end
    end
end
% b_{i0} column vectors in i-CS
p.b  = zeros(3,p.n);
for i = 1:p.n
    for j=1:p.n
        p.b(:,i) = p.b(:,i) + p.d(:,i,j)*p.mass(j);
    end
    p.b(:,i) = p.b(:,i)/sum(p.mass);
end


% Diagonal blocks of K tensor in i-CS
p.Kii = zeros(3,3,p.n);
for i=1:p.n
    p.Kii(:,:,i) = p.I(:,:,i);
    for k=1:p.n
        p.Kii(:,:,i) = p.Kii(:,:,i) + p.mass(k)*(p.d(:,i,k)'*p.d(:,i,k)*eye(3) - p.d(:,i,k)*p.d(:,i,k)');
    end
end

%%






