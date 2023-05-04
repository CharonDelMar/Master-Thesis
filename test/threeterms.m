%Initialize the parameters
T = 10 * pi ;
figure;
for sigma = -0.1:0.005:0.175
   %sigma = 0.5;
   ampl1 = 1 + 1/2 * sigma;
   ampl2 = 5/8 * sigma;
   ampl4 = -1/8 * sigma;
   
   damping = 0.4;
   
   natural_freq = 1;%omega
   freq = (0:0.25:1) * natural_freq;%Omega
   
   %construct the system
   syms x(t) zeta omega x1 x2 x4 Omega
   x(t) = x1 * cos(Omega * t) + x2 * cos(2 * Omega * t) + x4 * cos(4 * Omega * t);
   p = diff(x,t) * (diff(x,t,t) + 2 * zeta * omega * diff(x) + omega^2 * x);
   f = diff(x,t,t) + 2 * zeta * omega * diff(x) + omega^2 * x;
   F = - omega^2 * x;
   G = diff(x,t,t) + 2 * zeta * omega * diff(x);
   
   x(t,x1,x2,x4,Omega) = x;
   F(t,zeta,omega,x1,x2,x4,Omega) = F;
   G(t,zeta,omega,x1,x2,x4,Omega) = G;
  
   fplot(@(t)x(t,ampl1,ampl2,ampl4,freq(5)),@(t)F(t,damping,natural_freq,ampl1,ampl2,ampl4, freq(5)),[0 T])
   hold on;
   %for k = 2:numel(freq)
   %    fplot(@(t)x(t,ampl1,ampl2,ampl4,freq(k)),@(t)F(t,damping,natural_freq,ampl1,ampl2,ampl4,freq(k)),[0 T])
   %end
   %hold on;
   fplot(@(t)x(t,ampl1,ampl2,ampl4,freq(5)),@(t)G(t,damping,natural_freq,ampl1,ampl2,ampl4,freq(5)),[0 T])
   hold on;
   %for k = 2:numel(freq)
   %    fplot(@(t)x(t,ampl1,ampl2,ampl4,freq(k)),@(t)G(t,damping,natural_freq,ampl1,ampl2,ampl4,freq(k)),[0 T])
   %end
   hold on;
 end

 hold off;
 grid on;
 title('Work loop: x_1 = 1 + 1/2 * \sigma,x_2 = 5/8 * \sigma,x_4 = -1/8 * \sigma,\zeta = 0.4');
