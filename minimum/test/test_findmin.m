T = 4 * pi;

%sys1
syms t x1 x2 omega x(t)
x(t) = x1 * cos(omega * t) + x2 * cos(2 * omega * t);

x(t,x1,x2,omega) = x;
centernode = 1/2 * (2 * pi / 1 * 3/4 + 2 * pi /1 * 3/4);

figure;
fplot(@(t)x(t,1,1,1),[0,T],LineStyle="-");
hold on;
fplot(@(t)x(t,1,0,1),[0,T],LineStyle=":");
hold on;
fplot(@(t)x(t,0,1,1),[0,T],LineStyle="-.");
hold off;
grid on;
lgd = legend('cos( t ) + cos( 2t )', 'cos( t )', 'cos( 2t )');
title('cos( t ) + cos( 2t )');

%sys2
syms t x1 x2 x3 omega x(t)
x(t) = x1 * cos(omega * t) + x2 * cos(2 * omega * t) + x3 * cos(3 * omega * t);

x(t,x1,x2,x3, omega) = x;

figure;
fplot(@(t)x(t,1,1,1,1),[0,T],LineStyle="-");
hold on;
fplot(@(t)x(t,1,0,0,1),[0,T],LineStyle=":");
hold on;
fplot(@(t)x(t,0,1,0,1),[0,T],LineStyle="-.");
hold on;
fplot(@(t)x(t,0,0,1,1),[0,T],LineStyle="--");
hold off;
grid on;
lgd = legend('cos( t ) + cos( 2t ) + cos( 3t )', 'cos( t )', 'cos( 2t )','cos( 3t )');
title('cos( t ) + cos( 2t ) + cos( 3t )');

%sys3
syms t x1 x2 x3 x4 omega x(t)
x(t) = x1 * cos(omega * t) + x2 * cos(2 * omega * t) + x3 * cos(3 * omega * t) + x4 * cos( 4 * omega * t);

x(t,x1,x2,x3,x4, omega) = x;

figure;
fplot(@(t)x(t,1,1,1,1,1),[0,T],LineWidth=2,Color='b');
hold on;
fplot(@(t)x(t,1,0,0,0,1),[0,T],LineStyle=":");
hold on;
fplot(@(t)x(t,0,1,0,0,1),[0,T],LineStyle="-.");
hold on;
fplot(@(t)x(t,0,0,1,0,1),[0,T],LineStyle="--");
hold on;
fplot(@(t)x(t,0,0,0,1,1),[0,T],LineStyle="-" );
hold off;
grid on;
lgd = legend('cos( t ) + cos( 2t ) + cos( 3t ) + cos( 4t )', 'cos( t )', 'cos( 2t )','cos( 3t )','cos( 4t )');
title('cos( t ) + cos( 2t ) + cos( 3t ) + cos( 4t )');



