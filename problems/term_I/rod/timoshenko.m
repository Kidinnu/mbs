

%
E = 10000;
A = 0.0001;
mu = 10.0;
L = 1.0;

k = 1:20;
a = (2*k-1)*pi*0.5/L;
p = a*sqrt(E*A/mu);

y0 = @(x) interp1(linspace(0,L,10),0.1*linspace(0,L,10)/L,x);

M  = @(i) (2/L)*integral(@(xx) y0(xx).*sin(a(i).*xx),0,L);
N  = @(i) 0;

yk = @(i,t,x) (M(i)*cos(p(i)*t)+N(i)*sin(p(i)*t))*sin((2*i-1)*pi*x*0.5/L);
dykdx = @(i,t,x) (M(i)*cos(p(i)*t)+N(i)*sin(p(i)*t))*cos((2*i-1)*pi*x*0.5/L)*((2*i-1)*pi*0.5/L);

y = @(t,x) sum(arrayfun(@(i) yk(i,t,x),k));
ep = @(t,x) sum(arrayfun(@(i) dykdx(i,t,x),k));

x = linspace(0,L,10);

figure('Position',[100 100 850 600]);
axis([0 1.2 -0.5 0.5]);
hold on;

cmap = colormap(jet(128));

maxd = ep(0,L/2)*1.5;
mind =-ep(0,L/2)*1.5;

def2index = @(d) floor((d-mind)/((maxd-mind)/128));

for t=0:0.1:100
    cla;
    if 1 == 1
        plot(x,arrayfun(@(xx) ep(t,xx),x));
    else
    def   = arrayfun(@(xx) y(t,xx),x);
    dydx  = arrayfun(@(xx) ep(t,xx),x);
    pos   = x + def;
    for j = 1:size(x,2)
        color = cmap(def2index(dydx(j)));
        if j==1 
            prev = 0; 
        else
            prev = pos(j-1); 
        end
        patch([prev pos(j) pos(j) prev], [0.1 0.1 -0.1 -0.1],color,'EdgeColor','White');
    end
    end
    getframe;
end
    
    

