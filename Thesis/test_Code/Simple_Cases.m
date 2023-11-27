clc;
clear all;

% I. Check P = f * dx/dt >0 with specific values

% Define function x = a_1 * cos(omega * t) + a_2 * cos(2 * omega * t) * ...
% f(dx,t) = ddx + 2 * zeta * omega_0 * dx + omega_0^2 * x
% p = f * dx
% omega = lamda * omega_0  * zeta

%Coefficient Series A = [a_1, a_2, ...]
%Trigonometric function Matrix Cn = [cos(omega * t), cos(2 * omega * t), ...]

%Initialization
m = 3; % m terms Case
Sup = 10; % Upper Bound of a_k eg.10
A = zeros(1,m); % eg. Define A = [a_1, a_2] (<=> Two terms)
Cn = zeros(m,1);% eg. Define cos(omega * t), cos(2 * omega * t) (<=> Two terms)
x_t = A * Cn; % eg. x = a_1 * cos(omega * t) + a_2 * cos(2 * omega * t) 
zeta = rand(1);
omega_0  = 10 * rand(1);
omega_num = 10; % eg. Generate 10 omega
omega = zeros(1,omega_num); % Define an array to store different omega
Period = zeros(1,omega_num);
time_step = 80;
step_size = zeros(1,omega_num);
x = zeros(1,time_step); %Initialize X array
f = {omega_num};
p = {omega_num};
X = {omega_num};
delta_t = 0; % delta_t = Period / 2 * time_step

for iter = 1 : m
    A(iter) = Sup * rand(1); %Initialize A matrix
end

% Compute omega and period T (T is also the time interval)
for iter = 1:omega_num
    lamda = iter;
    omega(iter) = lamda * omega_0  * zeta;
    Period(iter) = 2 * pi / omega(iter);
    step_size(iter) = 2 * Period(iter)/time_step;
end

% Compute X
for omega_index = 1:omega_num
    for step_num = 1:time_step
        for k = 1:m
            Cn(k) = cos(k * omega(omega_index) * step_size(omega_index) * step_num);
        end
        x_t = A * Cn;
        x(step_num) = x_t;
    end
    X{omega_index} = x;
end

%Derivatives DDX, DX
ddx = zeros(1,time_step);
dx = zeros(1,time_step);
Sn = zeros(m,1);
CCn = zeros(m,1);
dx_t = A * Sn;
ddx_t = A * CCn;
for omega_index = 1:omega_num
    for step_num = 1:time_step
        for k = 1:m
            Sn(k) = - k * omega(omega_index) * sin(k * omega(omega_index) * step_size(omega_index) * step_num);
            CCn(k) = - k^2 * omega(omega_index)^2 * cos(k * omega(omega_index) * step_size(omega_index) * step_num);
        end
        ddx_t = A * CCn;
        dx_t = A * Sn;
        ddx(step_num) = ddx_t;
        dx(step_num) = dx_t;
    end
    DDX{omega_index} = dx;
    DX{omega_index} = dx;
end

%Compute f and p
for iter = 1:omega_num
    f{iter} = DDX{iter} + 2 * zeta * omega_0 * DX{iter} + omega_0^2 * X{iter} ;
    p{iter} = f{iter} .* DX{iter};
end



% Plot
%X
figure;
for iter = 1:omega_num
    Time_series = linspace(0, 2 * Period(iter), time_step);
    plot(Time_series, X{iter});
    title(sprintf('X Curve with different omega, omega_0 = %f, zeta = %f', omega_0, zeta));
    xlabel('Time');
    ylabel('x');
    hold on;
end
hold off;


for iter = 1:omega_num
    figure;
    Time_series = linspace(0, 2 * Period(iter), time_step);
    plot(Time_series, X{iter});
    title(sprintf('Details: X Curve with omega = %f', omega(iter)));
    xlabel('Time');
    set(gca,'XAxisLocation','origin');
    xticks(linspace(0 , 2 * Period(iter), 9))
    xticklabels({'0','T/4','T\2','3/4T', 'T','5\4T','3\2T','7\4T','2T'});
    ylabel('x');
    hold on;
    %grid on;
end
hold off;

%f
figure;
for iter = 1:omega_num
    Time_series = linspace(0, 2 * Period(iter), time_step);
    plot(Time_series, f{iter});
    title(sprintf('f Curve with different omega, omega_0 = %f, zeta = %f', omega_0, zeta));
    xlabel('Time');
    ylabel('f');
    hold on;
end
hold off;


for iter = 1:omega_num
    figure;
    Time_series = linspace(0, 2 * Period(iter), time_step);
    plot(Time_series, f{iter});
    title(sprintf('Details: f Curve with omega = %f', omega(iter)));
    xlabel('Time');
    set(gca,'XAxisLocation','origin');
    xticks(linspace(0 , 2 * Period(iter), 9))
    xticklabels({'0','T/4','T\2','3/4T', 'T','5\4T','3\2T','7\4T','2T'});
    ylabel('f');
    hold on;
    %grid on;
end
hold off;

%p
figure;
for iter = 1:omega_num
    Time_series = linspace(0, 2 * Period(iter), time_step);
    plot(Time_series, p{iter});
    title(sprintf('p Curve with different omega, omega_0 = %f, zeta = %f', omega_0, zeta));
    xlabel('Time');
    ylabel('p');
    hold on;
end
hold off;

for iter = 1:omega_num
    figure;
    Time_series = linspace(0, 2 * Period(iter), time_step);
    plot(Time_series, p{iter});
    title(sprintf('Details: p Curve with omega = %f', omega(iter)));
    xlabel('Time');
    set(gca,'XAxisLocation','origin');
    xticks(linspace(0 , 2 * Period(iter), 9))
    xticklabels({'0','T/4','T\2','3/4T', 'T','5\4T','3\2T','7\4T','2T'});
    ylabel('p');
    hold on;
    %grid on;
end
hold off;
