clear all;
clc;

%Objective Function is p = f * Dx

%f = DDX + 2 * zeta * omega_0 * DX + omega_0^2 * X;
%p = f .* DX;

%for iter = 1:1:100
%Create table to store results
%CODE for simple case
sz = [400 3];
varTypes = ["double","double","cell"];
varNames = ["zeta","omega_ratio","solution"];
dat = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
% sz = [400 6];
% varTypes = ["double","double","cell","cell","cell","cell"];
% varNames = ["zeta","omega_ratio","solution","power","F","G"];
% dat = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

%Initialization
%data = readtable('C:\Users\m1352\Desktop\Master Thesis background\dataframe\numerical_continuation_benchmark.csv');

m = 3; % m terms Case

omega_0 = 1;
table_index = 1;
for zeta= 0.01:0.01:1
    iterflag = 1;
    omega = omega_0;
    while omega < 2 % made an edit!   
        Period = 2 * pi / omega;
        time_step = 1000;
        t = linspace(0,Period,time_step);
        
        % Print
        fprintf(['zeta ' num2str(zeta) ' omega ' num2str(omega)])
        
        Cnt = cosharm(m,omega,time_step,t);
        Snt = sinharm(m,omega,time_step,t);
        CCnt = ccosharm(m,omega,time_step,t);
        
        digits(9);
        
        %Create Variables
        a = optimvar('a',m); % m-by-1 variable
        
        %Convert the function into optimization expression
         [min_p,max_p,p,f,x,minx,maxx,F,G] = fcn2optimexpr(@power,zeta,omega_0,a,Cnt,Snt,CCnt,m,time_step);
         funminp= @(p)min(p);
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
               'Algorithm', 'sqp', ...
               'Display','iter'); %,'PlotFcn','optimplotfval');
        
            %Solve the problem
            %x0.a = [a10,a20,a30];
             if iterflag ==1
                    x0.a = [1,0,0];

             else
                    x0.a = [a1,a2,a3];

             end
              
              %show(prob);
              [sol,fval,exitflag,output] = solve(prob,x0, 'Options', options);
              a1 = sol.a(1);
              a2 = sol.a(2);
              a3 = sol.a(3);
              p_val = evaluate(min_p,sol);
%               p_dat = evaluate(p,sol);
%               F_dat = evaluate(F,sol);
%               G_dat = evaluate(G,sol);
              
              %Export Data to EXCEL
              %Check out whether the algorithm gave us a nice solution
              %x dx ddx
    
              if  abs(p_val) <= 1e-9
                                        
                      %dat(table_index,:) = {zeta,omega/omega_0,sol.a,{p_dat},{F_dat},{G_dat}};
                      dat(table_index,:) = {zeta,omega/omega_0,sol.a};
                      table_index = table_index + 1;
                      iterflag = iterflag + 1;                     
              else
                      sol_test_FLAG = 0;
                      break;
              end
          omega = omega + 0.01;
    end
end
writetable(dat,'C:\Users\m1352\Desktop\Master Thesis background\dataframe\numerical_continuation_sqp.csv')

table_index = 1;

for zeta= 0.01:0.01:1
    iterflag = 1;
    omega = omega_0;
    while omega > 0 % made an edit!   
        Period = 2 * pi / omega;
        time_step = 1000;
        t = linspace(0,Period,time_step);
        
        % Print
        fprintf(['zeta ' num2str(zeta) ' omega ' num2str(omega)])
        
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
               'Algorithm', 'sqp', ...
               'Display','iter'); %,'PlotFcn','optimplotfval');
        
            %Solve the problem
            %x0.a = [a10,a20,a30];
             if iterflag ==1
                    x0.a = [1,0,0];

             else
                    x0.a = [a1,a2,a3];
             end
              
              %show(prob);
              [sol,fval,exitflag,output] = solve(prob,x0, 'Options', options);
              a1 = sol.a(1);
              a2 = sol.a(2);
              a3 = sol.a(3);
              p_val = evaluate(min_p,sol);
%               p_dat = evaluate(p,sol);
%               F_dat = evaluate(F,sol);
%               G_dat = evaluate(G,sol);
              
              %Export Data to EXCEL
              %Check out whether the algorithm gave us a nice solution
              %x dx ddx
   
              if  abs(p_val) <= 1e-9
                      sol_test_FLAG = 1;
                      
                      %dat(table_index,:) = {zeta,omega/omega_0,sol.a,{p_dat},{F_dat},{G_dat}};
                      dat(table_index,:) = {zeta,omega/omega_0,sol.a};
                      
                      table_index = table_index + 1;
                      iterflag = iterflag + 1;
              else
                      sol_test_FLAG = 0;
                      break;
              end
        omega = omega - 0.01;
    end
end

writetable(dat,'C:\Users\m1352\Desktop\Master Thesis background\dataframe\numerical_continuation_sqp2.csv')

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
