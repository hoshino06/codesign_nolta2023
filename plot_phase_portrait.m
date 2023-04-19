% Draw phase portraits of the closed-loop system 

clear; clc;

%% Figure1: Initial system 
figure(1)
clf
width =  500;
height = 380;
set(gcf, "Position", [1000 500 width height])
fontsize(gca,15,"pixels")
box;
hold on
xlim([-1,1])
ylim([-1,1])
xlabel('x1')
ylabel('x2')

load("results.mat")
rho = rho_list(1);
u   = u_init;
sys = matlabFunction( f(rho) + g(rho)*u );
fun = @(t,x) sys(x(1), x(2));
hold on;
[X,Y] = meshgrid(-1:0.5:1, -1:2:1);
for i = 1:numel(X)
  [t,y] = ode45(fun,[0,10], [X(i); Y(i)]);
  plot(y(:,1), y(:,2))  
end
hold off;

saveas(gcf,'plot/portrait_init','epsc')


%% Figure2: Initial system 
figure(2)
clf
width =  500;
height = 380;
set(gcf, "Position", [1000 60 width height])
fontsize(gca,15,"pixels")
box;
hold on
xlim([-1,1])
ylim([-1,1])
xlabel('x1')
ylabel('x2')

rho = rho_list(end);
C_N   = (CN_list(end,:)).';
V_N   = Phi*C_N;

u   = u_init;
sys = matlabFunction( f(rho) + g(rho)*u );
fun = @(t,x) sys(x(1), x(2));
hold on;
[X,Y] = meshgrid(-1:0.5:1, -1:2:1);
for i = 1:numel(X)
  [t,y] = ode45(fun,[0,10], [X(i); Y(i)]);
  plot(y(:,1), y(:,2))  
end
hold off;

saveas(gcf,'plot/portrait_result','epsc')
