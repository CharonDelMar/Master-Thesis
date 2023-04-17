T = 2 * pi;

%sys1
syms t alpha t_guess Dx x(t) x_sub1(t) x_sub2(t) Dx_sub1(t) Dx_sub2(t)
x(t) = cos(t) + cos(2 * t);
Dx(t) =  - diff(x,t);
x_sub1(t) = cos(t);
Dx_sub1(t) = - diff(x_sub1,t);
x_sub2(t) = cos(2 * t);
Dx_sub2(t) = - diff(x_sub2,t);

t_inf = solve(Dx == 0, t);
t_inf_sub1 = solve(Dx_sub1 == 0,t);
t_inf_sub2 = solve(Dx_sub2 == 0,t);


figure;
fplot(@(t)x(t),[0,T],LineWidth=2,Color = '#B22222');
hold on;
fplot(@(t)Dx(t),[0,T],LineStyle="-",LineWidth=1.5,Color = '#3D59AB');
hold on;
fplot(@(t)x_sub1(t),[0,T],LineStyle="-.",Color='#385E0F');
hold on;
fplot(@(t)x_sub2(t),[0,T],LineStyle="-.",Color='#32CD32');
hold on;
plot(t_inf,0,'.', Markersize = 20, Color = 'k');
hold on;
plot(t_inf,x(t_inf),'.' ,Markersize = 20,Color = 'k');
hold on;
for i = 1:1:numel(t_inf)
    line([t_inf(i),t_inf(i)],[x(t_inf(i)),0],LineStyle = '--',Linewidth = 1.5 ,Color = '#734A12');
    hold on;
end
hold on;

plot(t_inf_sub1,0,'*', Markersize = 10, Color = '#385E0F');
hold on;
plot(t_inf_sub1,x_sub1(t_inf_sub1),'*' ,Markersize = 10,Color = '#385E0F');
hold on;
plot(t_inf_sub1 + pi,0,'*', Markersize = 10, Color = '#385E0F');
hold on;
plot(t_inf_sub1 + pi,x_sub1(t_inf_sub1 + pi),'*' ,Markersize = 10,Color = '#385E0F');
hold on;
for i = 0:1:numel(t_inf_sub1)
    line([t_inf_sub1 + pi * i ,t_inf_sub1 + pi * i],[x_sub1(t_inf_sub1 + pi * i),0],LineStyle = ':', Color = '#385E0F');
    hold on;
end
hold on;

plot(t_inf_sub2,0,'*', Markersize = 10, Color = '#32CD32');
hold on;
plot(t_inf_sub2,x_sub2(t_inf_sub2),'*' ,Markersize = 10,Color = '#32CD32');
hold on;
plot(t_inf_sub2 + pi/2,0,'*', Markersize = 10, Color = '#32CD32');
hold on;
plot(t_inf_sub2 + pi/2,x_sub2(t_inf_sub2 + pi/2),'*' ,Markersize = 10,Color = '#32CD32');
hold on;
for i = 0:1:numel(t_inf_sub2)
    line([t_inf_sub2 + pi/2 * i ,t_inf_sub2 + pi/2 * i],[x_sub2(t_inf_sub2 + pi/2 * i),0],LineStyle = ':',Color = '#32CD32');
    hold on;
end

hold off;
grid on;
set(gca,'XAxisLocation','origin'); %将x轴的位置设置在y=0处。
set(gca,'YAxisLocation','origin'); %将y轴的位置设置在x=0处。
lgd = legend('cos(t) + cos(2t)', 'derivative','cos(t)','cos(2t)');
title('Primal Function: cos(t) + cos(2t)')
