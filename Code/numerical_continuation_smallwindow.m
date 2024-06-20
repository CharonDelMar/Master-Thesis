clear all;
clc;

%Objective Function is p = f * Dx

%f = DDX + 2 * zeta * omega_0 * DX + omega_0^2 * X;
%p = f .* DX;

%for iter = 1:1:100
%Create table to store results
sz = [400 5];
varTypes = ["double","double","cell","double","int64"];
varNames = ["zeta","omega_ratio","solution","power","test_FLAG"];
dat = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

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
        t = linspace(0,2 * Period,time_step);
        
        % Print
        fprintf(['zeta ' num2str(zeta) ' omega ' num2str(omega)])
        
        Cnt = cosharm(m,omega,time_step,t);
        Snt = sinharm(m,omega,time_step,t);
        CCnt = ccosharm(m,omega,time_step,t);
        
        digits(9);
        
        %Create Variables
        a = optimvar('a',m); % m-by-1 variable
        
        %Convert the function into optimization expression
         [min_p,max_p,p,f,x,minx,maxx,F_minx,G_minx,F_maxx,G_maxx] = fcn2optimexpr(@power,zeta,omega_0,a,Cnt,Snt,CCnt,m,time_step);
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
                
              %Export Data to EXCEL
              %Check out whether the algorithm gave us a nice solution
              %x dx ddx
    
              if  fval < 1e-6
                      sol_test_FLAG = 1;
                      dat(table_index,:) = {zeta,omega/omega_0,sol.a,-fval,sol_test_FLAG};   
                      table_index = table_index + 1;
                      iterflag = iterflag + 1;
                      omega = omega + 0.01;
              else
                      sol_test_FLAG = 0;
                      break;
              end
    end
end
%mkdir('C:\Users\m1352\Desktop\Master Thesis background\dataframe\dataset\',num2str(14))
writetable(dat,'C:\Users\m1352\Desktop\Master Thesis background\dataframe\numerical_continuation_sqp_2.csv')
%end

for zeta= 0.01:0.01:1
    iterflag = 1;
    omega = omega_0;
    while omega < 2 % made an edit!   
        Period = 2 * pi / omega;
        time_step = 1000;
        t = linspace(0,2 * Period,time_step);
        
        % Print
        fprintf(['zeta ' num2str(zeta) ' omega ' num2str(omega)])
        
        Cnt = cosharm(m,omega,time_step,t);
        Snt = sinharm(m,omega,time_step,t);
        CCnt = ccosharm(m,omega,time_step,t);
        
        %a10 =  data.solution_1(iterflag);
        %a20 =  data.solution_2(iterflag);
        %a30 =  data.solution_3(iterflag);
            
        digits(9);
        
        %Create Variables
        a = optimvar('a',m); % m-by-1 variable
        
        %Convert the function into optimization expression
         [min_p,max_p,p,f,x,minx,maxx,F_minx,G_minx,F_maxx,G_maxx] = fcn2optimexpr(@power,zeta,omega_0,a,Cnt,Snt,CCnt,m,time_step);
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
           
           % Didnt work so well: default interior point
           % Worked well: sqp and active-set
           
           options = optimoptions(...
               @fmincon, ...
               'ConstraintTolerance', 1e-9, ...
               'OptimalityTolerance', 1e-9, ...
               'StepTolerance', 1e-10, ...
               'MaxFunctionEvaluations', 500, ...
               'Algorithm', 'active-set', ...
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
                
              %Export Data to EXCEL
              %Check out whether the algorithm gave us a nice solution
              %x dx ddx
    
              if  fval < 1e-6
                      sol_test_FLAG = 1;
                      dat(table_index,:) = {zeta,omega/omega_0,sol.a,-fval,sol_test_FLAG};   
                      table_index = table_index + 1;
                      iterflag = iterflag + 1;
                      omega = omega + 0.01;
              else
                      sol_test_FLAG = 0;
                      break;     
              end                     
    end
end
%mkdir('C:\Users\m1352\Desktop\Master Thesis background\dataframe\dataset\',num2str(14))
writetable(dat,'C:\Users\m1352\Desktop\Master Thesis background\dataframe\numerical_continuation_activeset_2.csv')
%end





%%Self-Function Part
function [min_p,max_p,p,f,x,minx,maxx,F_minx,G_minx,F_maxx,G_maxx] = power(zeta,omega_0,a,Cnt,Snt,CCnt,m,time_step)
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
