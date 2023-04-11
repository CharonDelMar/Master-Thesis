%Initialize the Parameters
T = 10 * pi ;

ampl = 1;
damping = 0.8;
natural_freq = 1;
freq = (0:0.25:1) * natural_freq;

%Construct the system
%syms x(t) x0 Omega
%x(t) = x0 * cos(Omega * t);

syms x(t) zeta omega x0 Omega
x(t) = x0 * cos(Omega * t);
p = diff(x,t) * (diff(x,t,t) + 2 * zeta * omega * diff(x) + omega^2 * x);
f = diff(x,t,t) + 2 * zeta * omega * diff(x) + omega^2 * x;

p(t,zeta,omega,x0,Omega) = p;
f(t,zeta,omega,x0,Omega) = f;

figure;
fplot(@(t)p(t,damping,natural_freq,ampl,freq(1)),[0 T])
hold on;
for k = 2:numel(freq)
    fplot(@(t)p(t,damping,natural_freq,ampl,freq(k)),[0 T])
end
hold off;
grid on;
xticks(T*linspace(0,1,5));
xticklabels({'0','\pi','2\pi','3\pi','4\pi'});
xlabel('t / \omega');
ylabel('f');
lgd = legend('0','0.25\omega','0.5\omega','0.75\omega','\omega');
%lgd = legend('0','0.1\omega','0.2\omega','0.3\omega','0.4\omega','0.5\omega','0.6\omega','0.7\omega','0.8\omega','0.9\omega','\omega');
title(lgd,'\Omega');
title('p-t');

x(t,x0,Omega) = x;
figure;
fplot(@(t)x(t,ampl,freq(1)),@(t)f(t,damping,natural_freq,ampl,freq(1)),[0 T])
hold on;
for k = 2:numel(freq)
    fplot(@(t)x(t,ampl,freq(k)),@(t)f(t,damping,natural_freq,ampl,freq(k)),[0 T])
end
hold off;
grid on;
lgd = legend('0','0.25\omega','0.5\omega','0.75\omega','\omega');
%lgd = legend('0','0.1\omega','0.2\omega','0.3\omega','0.4\omega','0.5\omega','0.6\omega','0.7\omega','0.8\omega','0.9\omega','\omega');
title(lgd,'\Omega');
title('Work loop');
