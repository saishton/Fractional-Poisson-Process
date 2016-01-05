function [] = randomDLD(startTime,X,Y)

%==BEST FIT CCDFs==%
fo_ml = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[1 1],'StartPoint',[0.5 0.5]);
ft_ml = fittype('mlf(beta,1,-gamma*x.^beta,6)','options',fo_ml);
[cf_ml,gof_ml] = fit(X',Y',ft_ml);
cv_ml = coeffvalues(cf_ml);

fo_gp = fitoptions('Method', 'NonlinearLeastSquares','Lower',[-inf -inf 0],'StartPoint',[0.5 0.5 0.5]);
ft_gp = fittype('gpcdf(x,k,sigma,theta,''upper'')','options',fo_gp);
[cf_gp,gof_gp] = fit(X',Y',ft_gp);
cv_gp = coeffvalues(cf_gp);

fo_wb = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'StartPoint',[0.5 0.5]);
ft_wb = fittype('wblcdf(x,a,b,''upper'')','options',fo_wb);
[cf_wb,gof_wb] = fit(X',Y',ft_wb);
cv_wb = coeffvalues(cf_wb);

%==CALCULATE CCDFs==%
beta = cv_ml(1);
gamma = cv_ml(2);
ccdf_ml = mlf(beta,1,-gamma*X.^beta,6);

k = cv_gp(1);
sigma = cv_gp(2);
theta = cv_gp(3);
ccdf_gp = gpcdf(X,k,sigma,theta,'upper');

a = cv_wb(1);
b = cv_wb(2);
ccdf_wb = wblcdf(X,a,b,'upper');

%==GoF Tests==%
test_data = sort(X)';

z_ml = ones(length(test_data),1)-mlf(beta,1,-gamma*test_data.^beta,6);
z_gp = gpcdf(test_data,k,sigma,theta);
z_wb = wblcdf(test_data,a,b);

stats_ml = testStatistics(test_data,z_ml);
stats_gp = testStatistics(test_data,z_gp);
stats_wb = testStatistics(test_data,z_wb);

stats_ml.Root_MSE = gof_ml.rmse;
stats_gp.Root_MSE = gof_gp.rmse;
stats_wb.Root_MSE = gof_wb.rmse;
stats_ml.R_Squared = gof_ml.rsquare;
stats_gp.R_Squared = gof_gp.rsquare;
stats_wb.R_Squared = gof_wb.rsquare;

%==PLOT AND SAVE==%
rDLD = figure();
hold on
plot(X,Y,'o');
plot(X,ccdf_ml);
plot(X,ccdf_gp);
plot(X,ccdf_wb);
set(gca,'XScale','log');
set(gca,'YScale','log');
axis([1E1,1E4,1E-5,1E0]);
xlabel('Contact Time (s)');
ylabel('CCDF');
legend('Random DLD','ML','Gen. Pareto','Weibull');
hold off

imagefilename = [startTime,'/randomDLD-img.png'];
print(imagefilename,'-dpng')
close(rDLD);

%==Save Relevant Data==%
struc_ml = struct('Stability',beta,'Scale',gamma);
struc_gp = struct('Shape',k,'Scale',sigma,'Location',theta);
struc_wb = struct('Scale',a,'Shape',b);

ML = struct('Parameters',struc_ml,'Statistics',stats_ml);
GP = struct('Parameters',struc_gp,'Statistics',stats_gp);
WB = struct('Parameters',struc_wb,'Statistics',stats_wb);

datafilename = [startTime,'/randomDLD-data.mat'];
save(datafilename,'ML','GP','WB')