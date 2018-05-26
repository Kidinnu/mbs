function A = Az(angle)

c = cos(angle);
s = sin(angle);

A = [c, -s,  0;
     s,  c,  0;
     0,  0,  1];


end

