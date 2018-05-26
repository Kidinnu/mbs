function da = kinematicEq123(a, w)
     c3 = cos(a(3));
     s3 = sin(a(3));
     c2 = cos(a(2));
     s2 = sin(a(2));
     
     da = [ w(1)*c3/c2 - w(2)*s3/c2;
            w(1)*s3 + w(2)*c3;
           -w(1)*c3*s2/c2 + w(2)*s3*s2/c2 + w(3)];

end

