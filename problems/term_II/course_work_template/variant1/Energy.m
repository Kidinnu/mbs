function [T, U] = Energy(q, p)
% Get Kinetic & potential energy 
%
%
r  = zeros(3, p.n, size(q,1));
w  = zeros(3, p.n, size(q,1));
V  = zeros(3, p.n, size(q,1));
T  = zeros(size(q,1),1);
U  = zeros(size(q,1),1);

for k=1:size(q,1)
    phi  = @(i) q(k,p.iq(i,1) : p.iq(i,1)+p.iq(i,2)-1);
    dphi = @(i) q(k,p.iq(i,1)+p.N : p.iq(i,1)+p.iq(i,2)+p.N-1);
    A0i = zeros(3,3,p.n);
    for i=1:p.n
        A0i(:,:,i) = eye(3);
        for a=i:-1:1        
            if p.T(a,i) ~= 0
                A0i(:,:,i) = p.A{a}(phi(a))*A0i(:,:,i);
            end
        end
    end
    for i=1:p.n
        for j=1:p.n
            w(:,i,k) = w(:,i,k) - A0i(:,:,j)*(p.T(j,i)*p.Wr{j}(phi(j),dphi(j)));
            r(:,i,k) = r(:,i,k) - A0i(:,:,j)*p.d(:,j,i);
            V(:,i,k) = V(:,i,k) - cross(w(:,j,k),A0i(:,:,j)*p.d(:,j,i));
        end
        T(k) = T(k) + 0.5*w(:,i,k)'*A0i(:,:,i)*p.I(:,:,i)*A0i(:,:,i)'*w(:,i,k);
        T(k) = T(k) + 0.5*V(:,i,k)'*V(:,i,k)*p.mass(i);
        U(k) = U(k) + p.mass(i)*9.807*r(3,i,k);
    end
end

end

