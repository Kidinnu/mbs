function A = Ay(angle)

c = cos(angle);
s = sin(angle);

A = [c,  0,  s;
     0,  1,  0;
    -s,  0,  c];


end

