%Initialize the parameters
T = 10 * pi ;

figure;
%for sigma = -1:0.1:1
for sigma = -0.089:0.005:0.081
    %sigma = 0.5;
    ampl1 = 1 + 9/16 * sigma;
    ampl2 = 7/12 * sigma;
    ampl3 = -1/16 * sigma;
    ampl4 = -1/12 * sigma;%x1,x2,x3,x4
    
    damping = 0.4;
    
    natural_freq = 1;%omega
    freq = (0:0.25:1) * natural_freq;%Omega
    
    %construct the system
    syms x(t) zeta omega x1 x2 x3 x4 Omega
    x(t) = x1 * cos(Omega * t) + x2 * cos(2 * Omega * t) + x3 * cos(3 * Omega * t) + x4 * cos(4 * Omega * t);
    p = diff(x,t) * (diff(x,t,t) + 2 * zeta * omega * diff(x) + omega^2 * x);
    f = diff(x,t,t) + 2 * zeta * omega * diff(x) + omega^2 * x;
    F = -omega^2 * x;
    G = diff(x,t,t) + 2*zeta * omega * diff(x);
    
    x(t,x1,x2,x3,x4,Omega) = x;
    F(t,zeta,omega,x1,x2,x3,x4,Omega) = F;
    G(t,zeta,omega,x1,x2,x3,x4,Omega) = G;
    
    %figure;
    fplot(@(t)x(t,ampl1,ampl2,ampl3,ampl4,freq(5)),@(t)F(t,damping,natural_freq,ampl1,ampl2,ampl3,ampl4,freq(5)),[0 T])
    hold on;

    fplot(@(t)x(t,ampl1,ampl2,ampl3,ampl4,freq(5)),@(t)G(t,damping,natural_freq,ampl1,ampl2,ampl3,ampl4,freq(5)),[0 T])
    hold on;

    %hold off;
    %grid on;
    %title('Work loop: x_1 = 1 + 9/16 \times \sigma,x_2 = 7/12 \times \sigma,x_3 = -1/16 \times \sigma, x_4 = -1/12 \times \sigma,\zeta = 0.4');

end

hold off;
grid on;
title('Work loop: x_1 = 1 + 9/16 \times \sigma,x_2 = 7/12 \times \sigma,x_3 = -1/16 \times \sigma, x_4 = -1/12 \times \sigma,\zeta = 0.4');
%(-0.089,0.081)
