%Initialize the parameters
T = 10 * pi ;

ampl1 = 1;
ampl2 = 0.5;
ampl3 = 0.2;%x1,x2,x3

damping = 0.8;

natural_freq = 1;%omega
freq = (0:0.25:1) * natural_freq;%Omega

%construct the system
syms x(t) zeta omega x1 x2 x3 Omega
x(t) = x1 * cos(Omega * t) + x2 * cos(2 * Omega * t) + x3 * cos(3 * Omega * t);
p = diff(x,t) * (diff(x,t,t) + 2 * zeta * omega * diff(x) + omega^2 * x);
f = diff(x,t,t) + 2 * zeta * omega * diff(x) + omega^2 * x;

x(t,x1,x2,x3,Omega) = x;
p(t,zeta,omega,x1,x2,x3,Omega) = p;
f(t,zeta,omega,x1,x2,x3,Omega) = f;

%plot
figure;
fplot(@(t)p(t,damping,natural_freq,ampl1,ampl2,ampl3,freq(1)),[0 T])
hold on;
for k = 2:numel(freq)
    fplot(@(t)p(t,damping,natural_freq,ampl1,ampl2,ampl3,freq(k)),[0 T])
end
hold off;
grid on;
xticks(T*linspace(0,1,5));
xticklabels({'0','\pi','2\pi','3\pi','4\pi'});
xlabel('t / \omega');
ylabel('f');
lgd = legend('0','0.25\omega','0.5\omega','0.75\omega','\omega');
title(lgd,'\Omega');
title('p-t: x_1 = 1,x_2 = 0.5,x_3 = 0.2, \zeta = 0.8');

figure;
fplot(@(t)x(t,ampl1,ampl2,ampl3,freq(1)),@(t)f(t,damping,natural_freq,ampl1,ampl2,ampl3,freq(1)),[0 T])
hold on;
for k = 2:numel(freq)
    fplot(@(t)x(t,ampl1,ampl2,ampl3,freq(k)),@(t)f(t,damping,natural_freq,ampl1,ampl2,ampl3,freq(k)),[0 T])
end
hold off;
grid on;
lgd = legend('0','0.25\omega','0.5\omega','0.75\omega','\omega');
%lgd = legend('0','0.1\omega','0.2\omega','0.3\omega','0.4\omega','0.5\omega','0.6\omega','0.7\omega','0.8\omega','0.9\omega','\omega');
title(lgd,'\Omega');
title('Work loop: x_1 = 1,x_2 = 0.5,x_3 = 0.2, \zeta = 0.8');

