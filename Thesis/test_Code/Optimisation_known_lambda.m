clc;
clear all;

% II. Optimisation Problem

%Initialization
m = 2; % m terms Case
inf_A = 0; % Lower Bound of a_k eg.10
zeta = rand(1);
omega_0  = 10 * rand(1);
Lambda_index = 10; % How many lambda we are interested in
Lambda = linspace(1,10,Lambda_index);
omega = Lambda * zeta * omega_0;

Period = zeros(1,Lambda_index);
for iter = 1:Lambda_index
    Period(iter) = 2 * pi / omega(iter);
end

time_step = 80;


CCnt = {m};
Cnt = {m}; 
Snt = {m};

%Construct the Problem

%Set up Variables
a = optimvar('a',m,'LowerBound',0,'Type','integer');
%t = optimvar('t','LowerBound',0,'UpperBound',2 * Period(2));%e.g. iter = 2

%Create an optimization problem and set the objective function
%prob = optimproblem;
%obj = p;%sum(a); %need to be reconsidered carefully
%prob.Objective = obj;

%Constraints
t = linspace(0,2 * Period(2),time_step);
for k = 1:m
    CCnt{k} = - k^2 * omega(2)^2 * cos(k * omega(2) * t);
    Cnt{k} = cos(k * omega(2) * t);
    Snt{k} = - k * omega(2) * sin(k * omega(2) * t);
end

X = 0 * t;
DX = 0 * t;
DDX = 0 * t;
for k = 1:m
    X = X + (a(k) * Cnt{k});
    DX = DX + (a(k) * Snt{k});
    DDX = DDX + (a(k) * CCnt{k});
end

f = DDX + 2 * zeta * omega_0 * DX + omega_0^2 * X;
p = f .* DX; %>= 0;
%prob.Constraints.power = p;
prob = optimproblem;
obj = p;

%sum(a); %need to be reconsidered carefully
prob.Objective = obj;

problem = prob2struct(prob);

show(prob)

%%Solve the Problem
%sol = solve(prob);
%
%%Check feasibility
%