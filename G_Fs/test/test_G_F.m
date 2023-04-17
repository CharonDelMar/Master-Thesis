%Initialize the Parameters
T = 10 * pi ;

ampl = 1;
damping = 0.8;
natural_freq = 1;
freq = (0:0.25:1) * natural_freq;

%Construct the system
%%system G(x(t))
syms x(t) zeta omega x0 Omega
x(t) = x0 * cos(Omega * t);
G = diff(x,t,t) + 2 * zeta * omega * diff(x);

%%system -Fs(x(t))
syms x(t) zeta omega Omega
x(t) = x0 * cos(Omega * t);
Fs = - omega^2 * x;

x(t,x0,Omega) = x;
Dx(t,x0,Omega) = diff(x);
G(t,zeta,omega,x0,Omega) = G;
Fs(t,zeta,omega,x0,Omega) = Fs;

%plot
figure;
fplot(@(t)x(t,ampl,freq(1)),@(t)G(t,damping,natural_freq,ampl,freq(1)),[0 T])
hold on;
for k = 2:numel(freq)
    fplot(@(t)x(t,ampl,freq(k)),@(t)G(t,damping,natural_freq,ampl,freq(k)),[0 T])
end
hold on;

fplot(@(t)x(t,ampl,freq(1)),@(t)Fs(t,damping,natural_freq,ampl,freq(1)),[0 T])
hold on;
for k = 2:numel(freq)
    fplot(@(t)x(t,ampl,freq(k)),@(t)Fs(t,damping,natural_freq,ampl,freq(k)),[0 T])
end

fplot(@(t)x(t,ampl,freq(1)),@(t)Dx(t,ampl,freq(1)),[0 T],LineStyleMode="manual",LineWidth=2)
hold on;
for k = 2:numel(freq)
    fplot(@(t)x(t,ampl,freq(k)),@(t)Dx(t,ampl,freq(1)),[0 T],LineStyleMode="manual",LineWidth=2)
end
hold off;

grid on;
lgd = legend('0','0.25\omega','0.5\omega','0.75\omega','\omega');
title(lgd,'\Omega');
title('G-Fs:Work loop');
