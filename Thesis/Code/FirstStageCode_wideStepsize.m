clear all;
clc;

%Objective Function is p = f * Dx

%f = DDX + 2 * zeta * omega_0 * DX + omega_0^2 * X;
%p = f .* DX;

%Create table to store results
sz = [400 5];
varTypes = ["double","double","cell","double","string"];
varNames = ["zeta","omega_ratio","solution","power","test_FLAG"];
dat = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

%Initialization
m = 5; % m terms Case

omega_0 = 1;
table_index = 1;
for zeta = 0.01 : 0.05 : 1.0
    for omega = 0.01 : 0.05 : 1.0
        
        Period = 2 * pi / omega;
        time_step = 1000;
        t = linspace(0,2 * Period,time_step);
        
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
        prob.Objective =  - min_p;
        
        %Define constraints
        prob.Constraints.normalisation = maxx-minx == 2;
        
        %Tolerance
        %OptimalityTolerance = 1e-9;
        options = optimoptions(@fmincon,'ConstraintTolerance', 1e-13, 'StepTolerance', 1e-13,'PlotFcn','optimplotfval');
        
        %Solve the problem
        %x0.a = [0.1 0.1 0.1 0.1 0.1 0.1];
        
        x0.a = 0.1*ones(1,m);
        %show(prob);
        [sol,fval,exitflag,output] = solve(prob,x0);
        
        %Export Data to EXCEL
        %Check out whether the algorithm gave us a nice solution
        %x dx ddx
        x_data = zeros(1,time_step);
        dx_data = zeros(1,time_step);
        ddx_data = zeros(1,time_step);
        
        for term_index = 1:m
            x_data = x_data + sol.a(term_index) * Cnt(term_index,:);
            dx_data = dx_data + sol.a(term_index) * Snt(term_index,:);
            ddx_data = ddx_data + sol.a(term_index) * CCnt(term_index,:);
        end
        
        %f p
        f_data = ddx_data + 2 * zeta * omega_0 * dx_data + omega_0^2 * x_data;
        p_data = f_data .* dx_data;
        test_p_data = min(p_data);
    
        if test_p_data >= - 1e-6
            sol_test_FLAG = 'true';
        else
            sol_test_FLAG = 'false';
        end
    
        dat(table_index,:) = {zeta,omega/omega_0,sol.a,test_p_data,sol_test_FLAG};   
        table_index = table_index + 1;
    end
end
writetable(dat,'zeta_omega_fiveTerms_testdat.csv')






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
            Cnt(index_term,index_timepoint) = cos((index_term) * omega * time);
        end
    end
end

function [Snt] = sinharm(m,omega,time_step,t)
    Snt = zeros(m,time_step);
    
    for index_term = 1:m
        for index_timepoint = 1:time_step
            time = t(index_timepoint);
            Snt(index_term,index_timepoint) = - (index_term) * omega * sin((index_term) * omega * time);
        end
    end
end

function [CCnt] = ccosharm(m,omega,time_step,t)
    CCnt = zeros(m,time_step);
    
    for index_term = 1:m
        for index_timepoint = 1:time_step
            time = t(index_timepoint);
            CCnt(index_term,index_timepoint) = - (index_term)^2 * omega^2 * cos((index_term) * omega * time);
        end
    end
end
