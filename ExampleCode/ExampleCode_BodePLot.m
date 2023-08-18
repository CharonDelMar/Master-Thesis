clear
clc
close all

% %the solutions is completely different from the samples, should be careful
% %to use this version
% 
% %step 1: Generate a list of desired frequencies
% w = logspace(-1,2,300);%starts at 10^-1 = 0.1, ends at 10^2 = 100 rad/s
% 
% %step 2: Compute amplification and phase shift at each frequency
% alpha = (12-3.*w.^2)./(16-(31/4)*w.^2 + w.^4);%re([G(jw)])
% beta = ((-3/2).*w)./(16-(31/4)*w.^2 + w./4); %im([G(jw)])
% 
% magGjw = (alpha.^2 + beta.^2).^(1/2);%-->linear amplification factor
% angleGjw = atan2(beta, alpha); %phase(rad)
% 
% %Step 3: Convert amplification to dB
% magGjw_dB = 20*log10(magGjw);
% angleGjw_deg = rad2deg(angleGjw);
% 
% %Step 4: Plot on a log10 x-axis
% figure
% subplot(2,1,1)
% semilogx(w,magGjw_dB)
% grid on
% ylabel('Magnitude (dB)')
% title('Manually creating a bode plot')
% 
% subplot(2,1,2)
% semilogx(w,angleGjw_deg)
% grid on
% ylabel('Phase (deg)')
% xlabel('\omega(rad/s)')

%%Use  Matlabe's built in function'bode' (This version works for now, use this one)
%Create a transfer function object of the system
%                    3
%   G(s) = -----------------------
%             s^2 +(1/2)*s + 4
num = [3];
den = [1 (1/2) 4];

G = tf(num,den);
%Use bode' command to generate the plot
figure
bode(G)
grid on