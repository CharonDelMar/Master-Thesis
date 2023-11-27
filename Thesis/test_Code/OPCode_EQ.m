% evaluate function power() at sol.a

% waveforms: plot x(t), f(t), p(t)
% p(t) > 0 (within machine precision)


% switch out from random parameter input (zeta = 0.4, omega_0 = 1, omega =
% 0.9 (lambda = 0.9)

clc;
clear all;

%Objective Function is p = f * Dx

%f = DDX + 2 * zeta * omega_0 * DX + omega_0^2 * X;
%p = f .* DX;

%Initialization
m = 2; % m terms Case
zeta = 0.4;
%for zeta = 0.1 : 0.1 : 1.0
omega_0 = 1;
omega = 0.8;

Period = 2 * pi / omega;
time_step = 1000;
t = linspace(0,2 * Period,time_step);

Cnt = cosharm(m,omega,time_step,t);
Snt = sinharm(m,omega,time_step,t);
CCnt = ccosharm(m,omega,time_step,t);

%Create Variables
a = optimvar('a',m); % m-by-1 variable

%Convert the function into optimization expression
[min_p,max_p,p,f,x,minx,maxx,F_minx,G_minx,F_maxx,G_maxx] = fcn2optimexpr(@power,zeta,omega_0,a,Cnt,Snt,CCnt);
fun = @(p)min(p);
funexpr = fcn2optimexpr(fun,p,'OutputSize',[1,1]);
funmaxp = @(p)max(p);
funexprmaxp = fcn2optimexpr(funmaxp,p,'OutputSize',[1,1]);
funxmin = @(x)min(x);
funexprxmin = fcn2optimexpr(funxmin,x,'OutputSize',[1,1]);
funxmax = @(x)max(x);
funexprxmax = fcn2optimexpr(funxmax,x,'OutputSize',[1,1]);

%Define Objective Functions
prob = optimproblem;
prob.Objective =  - min_p;

%Define constraints
%prob.Constraints.positive = min_p >= 0;
%prob.Constraints.equalityminx = F_minx == G_minx;
%prob.Constraints.equalitymaxx = F_maxx == G_maxx;
prob.Constraints.normalisation = maxx-minx == 2;

%Tolerance
%OptimalityTolerance = 1e-9;
options = optimoptions('fmincon','ConstraintTolerance', 1e-13, 'StepTolerance', 1e-13)

%Solve the problem
x0.a = [1 -0.5];
%show(prob);
[sol,fval,exitflag,output] = solve(prob,x0);

%%Plot

%Data
%x dx ddx
x_data = sol.a(1) * Cnt(1,:) + sol.a(2) * Cnt(2,:);
dx_data = sol.a(1) * Snt(1,:) + sol.a(2) * Snt(2,:);
ddx_data = sol.a(1) * CCnt(1,:) + sol.a(2) * CCnt(2,:);
%f p
f_data = ddx_data + 2 * zeta * omega_0 * dx_data + omega_0^2 * x_data;
p_data = f_data .* dx_data;

%Plot the Waveforms
figure;

subplot(3,1,1);
plot(t,x_data);
title('WaveForm: x(t)');
xlabel('Time','FontName','BIZ UDGothic','FontSize',8);
set(gca,'XAxisLocation','origin');
xticks(linspace(0 , 2 * Period, 9));
xticklabels({'0','T/4','T\2','3/4T', 'T','5\4T','3\2T','7\4T','2T'});
ylabel('Displacement: x','FontName','BIZ UDGothic','FontSize',8,'FontWeight','bold');
grid on;

subplot(3,1,2);
plot(t,f_data);
title('WaveForm: f(t)');
xlabel('Time','FontName','BIZ UDGothic','FontSize',8);
set(gca,'XAxisLocation','origin');
xticks(linspace(0 , 2 * Period, 9));
xticklabels({'0','T/4','T\2','3/4T', 'T','5\4T','3\2T','7\4T','2T'});
ylabel('force: f','FontName','BIZ UDGothic','FontSize',8,'FontWeight','bold');
grid on;

subplot(3,1,3);
plot(t,p_data);
title('WaveForm: p(t)');
xlabel('Time','FontName','BIZ UDGothic','FontSize',8);
set(gca,'XAxisLocation','origin');
xticks(linspace(0 , 2 * Period, 9));
xticklabels({'0','T/4','T\2','3/4T', 'T','5\4T','3\2T','7\4T','2T'});
ylabel('Power: p','FontName','BIZ UDGothic','FontSize',8,'FontWeight','bold');
hold on;
%Mark the minumium
[pmin_data, min_point] = min(p_data);
tmin_data = t(min_point);
hold on
style = 'm.'; 
markersize = 15; %// change as needed
plot(t(min_point), pmin_data, style, 'markersize', markersize);
grid on;
% Text with coordinates of minimum
offset = -.05; %// vertical offset as a fraction of y-axis span. Change as needed.
text(t(min_point),pmin_data + diff(ylim)*offset,['(' num2str(t(min_point)) ',' num2str(pmin_data) ')'])
% Enlarge y axis so that text is properly seen, if offset is negative
ylim(ylim+[diff(ylim)*offset*(offset<0) 0])

sgtitle(["2-Terms Case Test Waveform ", "Parameters: \xi = " + zeta, ", \omega_0 = 1, \omega =0.9;","Constraint: p\geq 0"], ...
    'HorizontalAlignmen','left', ...
    'FontName','Times Roman', ...
    'FontWeight','bold', ...
    'FontSize',12);
hold off;

%Plot WorkLoop: f-x
%Data Remix
F = - omega_0^2 * x_data;
G = ddx_data + 2 * zeta * omega_0 * dx_data;

figure;
F_Curve = plot(x_data,F,'LineWidth',1);
hold on;
G_Curve = plot(x_data,G,'LineWidth',1);

title('Wookloop: (F,G)-x');
xlabel('x(t)','FontName','BIZ UDGothic','FontSize',12,'FontWeight','bold');
set(gca,'XAxisLocation','origin');
ylabel('(F,G)','FontName','BIZ UDGothic','FontSize',12,'FontWeight','bold');
%Mark the xmin and xmax
hold on;
xmin_data = min(x_data);
xmax_data = max(x_data);
hold on
style = 'm.'; 
markersize = 15; %// change as needed
plot(xmin_data, 0, style, 'markersize', markersize);
hold on;
plot(xmax_data, 0, style, 'markersize', markersize);
hold on;

Fxmin_data = - omega_0^2 * xmin_data;
Fxmax_data = - omega_0^2 * xmax_data;
line([xmin_data xmin_data],[0 Fxmin_data],'linestyle','--','Color','b','LineWidth',0.5);
hold on;
line([xmax_data xmax_data],[0 Fxmax_data],'linestyle','--','Color','b','LineWidth',0.5);
hold off;
legend([F_Curve,G_Curve],{'F','G'});
grid on;

%end

function [min_p,max_p,p,f,x,minx,maxx,F_minx,G_minx,F_maxx,G_maxx] = power(zeta,omega_0,a,Cnt,Snt,CCnt)
x = a(1) * Cnt(1,:) + a(2) * Cnt(2,:);
dx = a(1) * Snt(1,:) + a(2) * Snt(2,:);
ddx = a(1) * CCnt(1,:) + a(2) * CCnt(2,:);
f = ddx + 2 * zeta * omega_0 * dx + omega_0^2 * x;
p = f .* dx;
min_p = min(p);
max_p = max(p);

%added Equality Condition
[minx, minx_index] = min(x);
F_minx = -omega_0 ^2 * minx;
G_minx = ddx(minx_index) + 2 * zeta * omega_0 * dx(minx_index);
[maxx, maxx_index] = max(x);
F_maxx = -omega_0 ^2 * maxx;
G_maxx = ddx(maxx_index) + 2 * zeta * omega_0 * dx(maxx_index);
end

function [Cnt] = cosharm(m,omega,time_step,t)
    Cnt = zeros(m,time_step);
    
    for index_term = 1:m
        for index_timepoint = 1:time_step
            time = t(index_timepoint);
            Cnt(index_term,index_timepoint) = cos(index_term * omega * time);
        end
    end
end

function [Snt] = sinharm(m,omega,time_step,t)
    Snt = zeros(m,time_step);
    
    for index_term = 1:m
        for index_timepoint = 1:time_step
            time = t(index_timepoint);
            Snt(index_term,index_timepoint) = - index_term * omega * sin(index_term * omega * time);
        end
    end
end

function [CCnt] = ccosharm(m,omega,time_step,t)
    CCnt = zeros(m,time_step);
    
    for index_term = 1:m
        for index_timepoint = 1:time_step
            time = t(index_timepoint);
            CCnt(index_term,index_timepoint) = - index_term^2 * omega^2 * cos(index_term * omega * time);
        end
    end
end




