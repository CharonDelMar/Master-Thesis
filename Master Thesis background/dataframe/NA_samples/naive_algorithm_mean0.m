clear all;
clc;

%Objective Function is p = f * Dx

%f = DDX + 2 * zeta * omega_0 * DX + omega_0^2 * X;
%p = f .* DX;

%CODE for simple case
sz = [400 3];
varTypes = ["double","double","cell"];
varNames = ["zeta","omega_ratio","solution"];
dat = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
%Initialization
m = 3; % m terms Case

%Gaussian_Sampling_initialization
num_samples = 10; % Number of sampling points
sigma = 0.5; % Standard deviation
range = [-1, 1]; % Allowed range for sampling values

samples = discrete_gaussian_sampling(m, num_samples, sigma, range);
writematrix(samples,'C:\Users\m1352\Desktop\Master Thesis background\dataframe\samples05.csv')

omega_0 = 1;
table_index = 1;


for zeta= 0.01:0.01:1
    for omega = 0.01:0.01:2
        
        Period = 2 * pi / omega;
        time_step = 1000;
        t = linspace(0,Period,time_step);
        
        Cnt = cosharm(m,omega,time_step,t);
        Snt = sinharm(m,omega,time_step,t);
        CCnt = ccosharm(m,omega,time_step,t);
        
        digits(9);
        
        %Create Variables
        a = optimvar('a',m); % m-by-1 variable
        
        %Convert the function into optimization expression
        [min_p,max_p,p,f,x,minx,maxx,F,G] = fcn2optimexpr(@power,zeta,omega_0,a,Cnt,Snt,CCnt,m,time_step);
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
        prob.Objective =  - min_p / max_p; % question
        
        %Define constraints
        prob.Constraints.normalisation = maxx-minx == 2;
        prob.Constraints.a1 = a(1) * a(1) >= 0.25;
        
        %Tolerance
        %OptimalityTolerance = 1e-9;
        options = optimoptions(...
               @fmincon, ...
               'ConstraintTolerance', 1e-9, ...
               'OptimalityTolerance', 1e-9, ...
               'StepTolerance', 1e-10, ...
               'MaxFunctionEvaluations', 500, ...
               'Algorithm', 'sqp', ...
               'Display','iter'); %,'PlotFcn','optimplotfval');
        
        %Solve the problem
        
        for iter = 1:1:10
            x0.a = samples(iter,:);
            %show(prob);
            [sol,fval,exitflag,output] = solve(prob,x0,'Options', options);
            p_val = evaluate(min_p,sol);
            
            %Export Data to EXCEL
            if  abs(p_val) <= 1e-9     
                dat(table_index,:) = {zeta,omega/omega_0,sol.a};   
                table_index = table_index + 1;
                break;
            end
        end
        
    end
end
%mkdir('C:\Users\m1352\Desktop\Master Thesis background\dataframe\dataset\',num2str(14))
writetable(dat,'C:\Users\m1352\Desktop\Master Thesis background\dataframe\NA_sample_sqp05.csv')

for zeta= 0.01:0.01:1
    for omega = 0.01
        
        Period = 2 * pi / omega;
        time_step = 1000;
        t = linspace(0,Period,time_step);
        
        Cnt = cosharm(m,omega,time_step,t);
        Snt = sinharm(m,omega,time_step,t);
        CCnt = ccosharm(m,omega,time_step,t);
        
        digits(9);
        
        %Create Variables
        a = optimvar('a',m); % m-by-1 variable
        
        %Convert the function into optimization expression
        [min_p,max_p,p,f,x,minx,maxx,F,G] = fcn2optimexpr(@power,zeta,omega_0,a,Cnt,Snt,CCnt,m,time_step);
         funminp = @(p)min(p);
         funexpr = fcn2optimexpr(funminp,p,'OutputSize',[1,1]);
         funmaxp = @(p)max(p);
         funexprmaxp = fcn2optimexpr(funmaxp,p,'OutputSize',[1,1]);
         funminx = @(x)min(x);
         funexprminx = fcn2optimexpr(funminx,x,'OutputSize',[1,1]);
         funmaxx = @(x)max(x);
         funexprmaxx = fcn2optimexpr(funmaxx,x,'OutputSize',[1,1]);
         
        %Define Objective Functions
        prob = optimproblem;
        prob.Objective =  - min_p / max_p; % question
        
        %Define constraints
        prob.Constraints.normalisation = maxx-minx == 2;
        prob.Constraints.a1 = a(1) * a(1) >= 0.25;
        
        %Tolerance
        %OptimalityTolerance = 1e-9;
        options = optimoptions(...
               @fmincon, ...
               'ConstraintTolerance', 1e-9, ...
               'OptimalityTolerance', 1e-9, ...
               'StepTolerance', 1e-10, ...
               'MaxFunctionEvaluations', 500, ...
               'Algorithm', 'active-set', ...
               'Display','iter'); %,'PlotFcn','optimplotfval');
        
        %Solve the problem
        
        for iter = 1:1:10
            x0.a = samples(iter,:);
            %show(prob);
            [sol,fval,exitflag,output] = solve(prob,x0,'Options', options);
            p_val = evaluate(min_p,sol);
            
            %Export Data to EXCEL
            if  p_val < 1e-6
                dat(table_index,:) = {zeta,omega/omega_0,sol.a};   
                table_index = table_index + 1;
                break;
            end
        end
        
    end
end
%mkdir('C:\Users\m1352\Desktop\Master Thesis background\dataframe\dataset\',num2str(14))
writetable(dat,'C:\Users\m1352\Desktop\Master Thesis background\dataframe\NA_sample_activeset05.csv')




%%Self-Function Part
function [min_p,max_p,p,f,x,minx,maxx,F,G] = power(zeta,omega_0,a,Cnt,Snt,CCnt,m,time_step)
x = zeros(1,time_step);
dx = zeros(1,time_step);
ddx = zeros(1,time_step);

for aindex = 1:m
    x = x + a(aindex) * Cnt(aindex,:);
    dx = dx +a(aindex) * Snt(aindex,:);
    ddx = ddx + a(aindex) * CCnt(aindex,:);
end

f = ddx + 2 * zeta * omega_0 * dx + omega_0^2 * x;
p = f .* dx;
min_p = min(p);
max_p = max(p);

G = ddx + 2 * zeta * omega_0 * dx;
F = -omega_0^2 * x;

%added Equality Condition
maxx = max(x);
minx = min(x);

end

function [Cnt] = cosharm(m,omega,time_step,t)
    Cnt = zeros(m,time_step);
    
    for index_term = 1:m
        for index_timepoint = 1:time_step
            time = t(index_timepoint);
            Cnt(index_term,index_timepoint) = cos((2 *index_term -1) * omega * time);
        end
    end
end

function [Snt] = sinharm(m,omega,time_step,t)
    Snt = zeros(m,time_step);
    
    for index_term = 1:m
        for index_timepoint = 1:time_step
            time = t(index_timepoint);
            Snt(index_term,index_timepoint) = - (2 *index_term -1) * omega * sin((2 *index_term -1) * omega * time);
        end
    end
end

function [CCnt] = ccosharm(m,omega,time_step,t)
    CCnt = zeros(m,time_step);
    
    for index_term = 1:m
        for index_timepoint = 1:time_step
            time = t(index_timepoint);
            CCnt(index_term,index_timepoint) = - (2 *index_term -1)^2 * omega^2 * cos((2 *index_term -1) * omega * time);
        end
    end
end


function samples = discrete_gaussian_sampling(m, num_samples, sigma, range)
    % m: Dimension of the space
    % num_samples: Number of sampling points
    % sigma: Standard deviation of the Gaussian distribution
    % range: Allowed range for sampling values (e.g., [-1, 1])

    % Initialize variables
    samples = zeros(num_samples, m);
    
    % Iterate to get the required number of samples
    for i = 1:num_samples
        while true
            % Sample m-dimensional vector
            sample = normrnd(0, sigma, 1, m);
            % Check if the sample is within the specified range
            if all(sample >= range(1) & sample <= range(2))
                samples(i, :) = sample;
                break;
            end
        end
    end
end
