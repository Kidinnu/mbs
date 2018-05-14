
figure;
axis([-10 10 -10 10 -10 10]);
hold on;

box;
axis vis3d;

for i=1:360
    cla;
    draw_upper_stage([4;0;0], [0;0;0], Ay(i*pi/180.0));
    
    shading flat;
    lighting gouraud;
    light('Position',[0 0 10]);
    
    getframe;
end