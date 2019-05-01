function A = getA0(k,angles)
%GETA0 Get transformation matrix from k to 0
%
    A=eye(2);
    for i=1:k
        A = A*getA(angles(i));    
    end
end

