%Initialize the parameters -- change to two terms diff(x,t,t) =
%-w_{0}max(x) max --> t=0 what I need is to find the min
T = 10 * pi ;

ampl1 = 1;
ampl2 = 0.5;
ampl3 = 0.2;%x1,x2,x3

damping = 0.8;

natural_freq = 1;%omega
freq = (0:0.25:1) * natural_freq;%Omega

%construct the system
%%system G(x(t))
syms x(t) zeta omega x1 x2 x3 Omega
%x(t) = x1 * sin(Omega * t) + x2 * sin(2 * Omega * t) + x3 * sin(3 * Omega * t);
x(t) = x1 * cos(Omega * t) + x2 * cos(2 * Omega * t) + x3 * cos(3 * Omega * t);
G = diff(x,t,t) + 2 * zeta * omega * diff(x);

%%system -Fs(x(t))
syms x(t) zeta omega Omega
%x(t) = x1 * sin(Omega * t) + x2 * sin(2 * Omega * t) + x3 * sin(3 * Omega * t);
x(t) = x1 * cos(Omega * t) + x2 * cos(2 * Omega * t) + x3 * cos(3 * Omega * t);
Fs = - omega^2 * x;

x(t,x1,x2,x3,Omega) = x;
G(t,zeta,omega,x1,x2,x3,Omega) = G;
Fs(t,zeta,omega,x1,x2,x3,Omega) = Fs;

%plot
figure;
fplot(@(t)x(t,ampl1,ampl2,ampl3,freq(1)),@(t)G(t,damping,natural_freq,ampl1,ampl2,ampl3,freq(1)),[0 T],LineStyle="-.",LineWidth=2)
hold on;
for k = 2:numel(freq)
    fplot(@(t)x(t,ampl1,ampl2,ampl3,freq(k)),@(t)G(t,damping,natural_freq,ampl1,ampl2,ampl3,freq(k)),[0 T],LineStyle="-.",LineWidth=2)
end
hold on;

fplot(@(t)x(t,ampl1,ampl2,ampl3,freq(1)),@(t)Fs(t,damping,natural_freq,ampl1,ampl2,ampl3,freq(1)),[0 T])
hold on;
for k = 2:numel(freq)
    fplot(@(t)x(t,ampl1,ampl2,ampl3,freq(k)),@(t)Fs(t,damping,natural_freq,ampl1,ampl2,ampl3,freq(k)),[0 T])
end
hold off;

grid on;

lgd = legend('0','0.25\omega','0.5\omega','0.75\omega','\omega');
title(lgd,'\Omega');

%title('G-Fs_{three terms}:Work loop_{sin}');
title('G-Fs_{three terms}:Work loop_{cos}');





