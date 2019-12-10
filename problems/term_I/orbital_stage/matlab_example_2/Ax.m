function A = Ax(angle)

c = cos(angle);
s = sin(angle);

A = [1,  0,  0;
     0,  c, -s;
     0,  s,  c];


end

