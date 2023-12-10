%%spot figure plot
[dat]=readtable('/Users/congxiaozhang/Documents/Master Thesis/Thesis/Code/test data with wide stepsize/odd series/zeta_omega_twoTerms_testdat.csv');

m = 2;
figure;
for i = 1:1:400
    if isequal(dat.test_FLAG(i), {'false'}) == 1
        scatter(dat.zeta(i),dat.omega_ratio(i),[],[0.7686 0.6862 0.5098],'filled');
        hold on;
    else
        scatter(dat.zeta(i),dat.omega_ratio(i),[],[0.5725 0.4039 0.4666],'filled');
        hold on;
    end
end
xlabel('\xi','FontName','BIZ UDGothic','FontSize',8);
set(gca,'XAxisLocation','origin');
ylabel('\Omega/\Omega_0','FontName','BIZ UDGothic','FontSize',8,'FontWeight','bold');
sgtitle([m + "-Terms Case \xi - \Omega/\Omega_0 Odd Series;"], ...
    'HorizontalAlignmen','left', ...
    'FontName','Times Roman', ...
    'FontWeight','bold', ...
    'FontSize',12);
hold off;

