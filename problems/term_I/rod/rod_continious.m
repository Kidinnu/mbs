% Модуль Юнга Н/м2
E  = 1000000;
% Площадь поперечного сечения м2
A  = 0.001;
% Погонная масса kg/м
mu = 10.0;
% Длина стержная м
L  = 1.0;
% Количество учитываемых собственных форм
k  = 1:10;
a  = (2*k-1)*pi*0.5/L;
% Массив частот
p  = a*sqrt(E*A/mu); % рад/с
% Начальная деформация
% К концу стержня приложена растягивающая сила F
F  = 100; % Н
y0 = @(xx) interp1(linspace(0,L,10),F*linspace(0,L,10)/(E*A),xx);
% Постоянные, определяемые по начальным условиям (для частоты i)
M  = @(i)  (2/L)*integral(@(xx) y0(xx).*sin(a(i).*xx),0, L);
M  = @(i)  (2*F/(L*E*A))*sin(a(i)*L)/(a(i)^2.0);
% Начальная скорость деформаций равна нулю, поэтому (для частоты i)
N  = @(i) 0;
% Слагаемые общего решения
yk = @(i,t,x) (M(i)*cos(p(i)*t)+N(i)*sin(p(i)*t))*sin((2*i-1)*pi*x*0.5/L);
% Слагаемые деформации
ek = @(i,t,x) (M(i)*cos(p(i)*t)+N(i)*sin(p(i)*t))*cos((2*i-1)*pi*x*0.5/L)*(2*i-1)*pi*0.5/L;
% Общее решение y(t,x)
y   = @(t,x) sum(arrayfun(@(i) yk(i,t,x),k));
% Деформация
eks = @(t,x) sum(arrayfun(@(i) ek(i,t,x),k));
% ------------------------------------------------------------------------
% Рисуем
% ------------------------------------------------------------------------
% Видео
v = VideoWriter('rod_c.avi');
open(v);

% Разделяем балку на 30 частей (для раскраски)
x = linspace(0,L,30);
figure('Position',[100 100 1920 1080]);
axis([0 1.3 -0.5 0.5]);
set(gca,'FontSize',20);
hold on;
% Цветовая палитра JET
cmap = colormap(jet(128));
% Максимальная деформация
maxd =  0.2;
mind = -0.2;
def2index = @(d) floor((d-mind)/((maxd-mind)/size(cmap,1)));
for t=0:0.002:2
    cla;
    yx = arrayfun(@(xx) y(t,xx),x);    
    pos = x + yx;    
    defx = arrayfun(@(xx) eks(t,xx),x);
    if 1 == 2
        plot(x, arrayfun(@(xx) y(t,xx),x));
    else
        for j = 1:size(x,2)
            col = cmap( def2index(defx(j)),:);
            prev_point = 0;
            if j ~= 1
                prev_point = pos(j-1);
            end
            patch([prev_point pos(j) pos(j) prev_point],...
                [0.1 0.1 -0.1 -0.1],col,'EdgeColor','none');            
        end
        plot(pos, arrayfun(@(xx) eks(t,xx),x),'k-','LineWidth',2);
        xlabel('x, м');ylabel('y, м');box('on');
    end
    %plot([0 L L 0],[0.103 0.103 -0.101 -0.101],'k:');
    text(0.1,0.45,sprintf('T=%5.3f c. Длина стержня L=%5.3f м',t,yx(end)+L),'FontSize',20);    
    text(0.1,0.35,sprintf('Частоты (Гц): '),'FontSize',20);    
    text(0.1,0.31,sprintf('%3.1f | ',(2*pi./p).^-1),'FontSize',20);    
    frame = getframe;    
    writeVideo(v,frame);
end
close(v);