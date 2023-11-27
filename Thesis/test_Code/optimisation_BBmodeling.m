clc;
clear all;

%Objective Function is p = f * Dx

%f = DDX + 2 * zeta * omega_0 * DX + omega_0^2 * X;
%p = f .* DX;

%Initialization
m = 2; % m terms Case
zeta = rand(1);
omega_0  = 10 * rand(1);
%Lambda_index = 10; % How many lambda we are interested in
%Lambda = linspace(1,10,Lambda_index);
Lambda = 0.2;
omega = Lambda * zeta * omega_0;
Period = 2 * pi / omega;
time_step = 80;
t = linspace(0,2 * Period,time_step);

Cnt = cosharm(m,omega,time_step,t);
Snt = sinharm(m,omega,time_step,t);
CCnt = ccosharm(m,omega,time_step,t);

%Create Variables
a = optimvar('a',m,'LowerBound',0,'UpperBound',100); % m-by-1 variable

%Convert the function into optimization expression
[min_p,p,f,x] = fcn2optimexpr(@power,zeta,omega_0,a,Cnt,Snt,CCnt);
fun = @(p)min(p);
funexpr = fcn2optimexpr(fun,p,'OutputSize',[1,1]);

%Define Objective Functions
prob = optimproblem;
prob.Objective = - min_p;

%Define constraints
%prob.Constraints.postive = min_p >= 0;

%Solve the problem
x0.a = [0 0];

show(prob);

[sol,fval,exitflag,output] = solve(prob,x0);

function [min_p,p,f,x] = power(zeta,omega_0,a,Cnt,Snt,CCnt)
x = a(1) * Cnt(1,:) + a(2) * Cnt(2,:);
dx = a(1) * Snt(1,:) + a(2) * Snt(2,:);
ddx = a(1) * CCnt(1,:) + a(2) * CCnt(2,:);
f = ddx + 2 * zeta * omega_0 * dx + omega_0^2 * x;
p = f .* dx;
min_p = min(p);
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

%
% evaluate function power() at sol.a
%
% waveforms: plot x(t), f(t), p(t)
% p(t) > 0 (within machine precision)


% switch out from random parameter input (zeta = 0.4, omega_0 = 1, omega =
% 0.9 (lambda = 0.9)

% IF optimisation doesn't given the true answer for prob above
% > Check that solver isn't defaulting to trivial solution (all a = 0), maybe
% include constrain to have a ne 0
% Change initial conditions (if solve is close to true soln it might even
% ignore the trivial one)

% > turning down the tolerance (1e-9)
% > next step: think about... add equality conditions....? (matrix conditions)

% > particle swarm optimisation (PSO) (toolboxes etc) - Dont do right
% away!!! next next step (or beyond)

% extra pretty plots - work loop (f against x)