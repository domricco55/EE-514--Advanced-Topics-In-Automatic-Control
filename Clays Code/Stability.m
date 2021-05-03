%% EE 514 Module 3: Stability
close all; clear all; clc;
%% First example
N = 5;
x = linspace(-N,N,10); y = x;
[X,Y] = meshgrid(x,y);
U = -X-Y;
V = X-Y.^3;
quiver(X,Y,U,V)
%% Damped Pendulum
N = 100; eps = 0.001;
x1 = linspace(-pi+eps,pi-eps,N);
x2 = linspace(-2,2,N);
[X1,X2] = meshgrid(x1,x2);
f1 = X2;
f2 = -X2-sin(X1);
E = 1-cos(X1)+0.5*(X2.^2);
V = E+1-cos(X1)+0.5*(X1+X2).^2;
figure;
    contour3(X1,X2,E,'LineStyle','--')
    hold on;
    contour3(X1,X2,V,'LineStyle','-')
    xlabel('Angle'); ylabel('Speed')
    legend('System Energy','Ad-hoc Lyapunov',...
        'Location','SouthEastOutside')
%% Limit cycle
n = 5; N = 1000;
x1 = linspace(-n,n,N);
x2 = x1;
yp = sqrt(0.5*(10-x1.^4));
valid_set = find(abs(imag(yp))<0.1);
yp = real(yp(valid_set)); yn = -yp;    
% Pick a positive el.
el = 99.5;
[X1,X2] = meshgrid(x1,x2);
L = (X1.^4 +2*X2.^2 -10).^2;
[r,c] = find(L<el);
figure;
    scatter(x1(c),x2(r),'filled','MarkerFaceAlpha',.1,'MarkerEdgeColor','none');
    hold on;
    scatter(x1(valid_set),yp,10,'filled','MarkerFaceColor','k');
    scatter(x1(valid_set),yn,10,'filled','MarkerFaceColor','k');
    axis([-n,n,-n,n]); grid on;
    xlabel('$x_1$'); ylabel('$x_2$');
%% Limit cycle trajectories
x0 = [-.01;-.01];
f = @(t,X) [X(2)-(X(1).^7).*(X(1).^4 + 2*(X(2).^2)-10);
    -1*(X(1).^3)-3*(X(2).^5).*(X(1).^4 + 2*(X(2).^2)-10)];
[t,x] = ode45(f,[0,50000],x0);
scatter(x(:,1),x(:,2),10,'filled');
%% Variable Gradient example
n = 2; N = 30;
x1 = linspace(-n,n,N); x2=x1;
[X1,X2] = meshgrid(x1,x2);
f1 = -X1+2*(X1.^2).*X2;
f2 = -X2;
vg = figure;
    quiver(X1,X2,f1,f2,2.5); 
    hold on; grid on;
    plot(x1(x1<0),-1/2./abs(x1(x1<0)),'r');
    plot(x1(x1>0),1/2./x1(x1>0));
    axis([-n,n,-n,n])
%% Variable Gradient example (con't)    
x0 = [1.5;0.9];
[t,x] = ode45(@(t,X) [-X(1)+2*(X(1).^2).*X(2); -X(2)],[0,20],x0);
figure(vg)
    plot(x(:,1),x(:,2));