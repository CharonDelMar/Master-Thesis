%%spot figure plot
[dat]=readtable('/Users/congxiaozhang/Documents/Master Thesis/Thesis/Code/test data with wide stepsize/odd series/zeta_omega_threeTerms_testdat.csv');

m = 3;
figure;
image(dat.zeta,dat.omega_ratio);
colormap(gca,'parula');
xlabel('\xi','FontName','BIZ UDGothic','FontSize',8);
set(gca,'XAxisLocation','origin');
ylabel('\Omega/\Omega_0','FontName','BIZ UDGothic','FontSize',8,'FontWeight','bold');
sgtitle([m + "-Terms Case \xi - \Omega/\Omega_0 Odd Series;"], ...
    'HorizontalAlignmen','left', ...
    'FontName','Times Roman', ...
    'FontWeight','bold', ...
    'FontSize',12);
hold off;