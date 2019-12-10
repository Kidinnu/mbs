figure;
axis([-7 7 -7 7 -7 7]);
hold on;

box;
axis vis3d;

for i=1:360
    cla;    
    xc = 4;
    r  = [0;0;0];
    % Rotate about y 
    draw_upper_stage([xc;0;0], r, Ay(i*pi/180.0));
    shading flat;
    lighting gouraud;
    light('Position',[0 0 10]);    
    getframe;
end