function [] = plotMLcomparison(input_filename,structure)

% input = 'journal.pone.0136497.s002.CSV';
% structure = '%f %f %f %*s %*s';

iF = 'input';
oF = 'output';

startTime = datestr(now,'yyyymmddTHHMMSS');
clean_input = strrep(input_filename, '.', '');
dir_ref = [oF,'\',startTime,'_',clean_input];
mkdir(dir_ref);

input = ['../',iF,'/',input_filename];

%==CALCULATE CCDF FOR SAMPLE==%
fid = fopen(input);
rawdata = textscan(fid,structure,'Delimiter',',');
fclose(fid);

data = cell2mat(rawdata);
data(:,1) = data(:,1)-data(1,1);
lowestID = min(min(data(:,2)),min(data(:,3)));
data(:,2) = data(:,2)-lowestID+1;
data(:,3) = data(:,3)-lowestID+1;
number_rows = size(data,1);
parfor i=1:number_rows
    thisrow = data(i,:);
    col2 = thisrow(1,2);
    col3 = thisrow(1,3);
    if col2 > col3
        thisrow(1,2) = col3;
        thisrow(1,3) = col2;
        data(i,:) = thisrow;
    end
end

[~, order] = sort(data(:,3));
partsorteddata = data(order,:);

[~, order] = sort(partsorteddata(:,2));
sorteddata = partsorteddata(order,:);

times = zeros(1,number_rows);
j = 1;
k = 1;
numpeople = max(max(sorteddata(:,2)),max(sorteddata(:,3)));
total_interactions = 0;
interactions = zeros(1,numpeople);
step_vector = [20 0 0];
while j<number_rows+1
    ID1 = sorteddata(j,2);
    ID2 = sorteddata(j,3);
    interactions(ID1) = interactions(ID1)+1;
    interactions(ID2) = interactions(ID2)+1;
    total_interactions = total_interactions+1;
    contact_time = 20;
    current_row = sorteddata(j,:);
    if j == number_rows
        next_row = [0 0 0]; 
    else
        next_row = sorteddata(j+1,:);
    end
    while isequal(next_row,current_row+step_vector)
        contact_time = contact_time+20;
        j = j+1;
        current_row = sorteddata(j,:);
        if j == number_rows
            next_row = [0 0 0];
        else
            next_row = sorteddata(j+1,:);
        end
    end
    times(k) = contact_time;
    j = j+1;
    k = k+1;
end

times(k:end) = [];
interactions(interactions==0) = []; %Remove nodes where no interaction occurs
activityPot = interactions/total_interactions;
[F,X] = ecdf(times);
ccdf_data = 1-F;
max_time = X(end);

%==MODIFY TO REMOVE ERRORS==%
Xrem = [X(1);X(end-9:end)];
X = X(2:end-10);
ccdf_data = ccdf_data(2:end-10);

%==BEST FIT CCDFs==%
fo_ml = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[1 1],'StartPoint',[0.5 0.5]);
ft_ml = fittype('mlf(beta,1,-gamma*x.^beta,6)','options',fo_ml);
[cf_ml,gof_ml] = fit(X,ccdf_data,ft_ml);
cv_ml = coeffvalues(cf_ml);

fo_gp = fitoptions('Method', 'NonlinearLeastSquares','Lower',[-inf -inf 0],'StartPoint',[0.5 0.5 0.5]);
ft_gp = fittype('gpcdf(x,k,sigma,theta,''upper'')','options',fo_gp);
[cf_gp,gof_gp] = fit(X,ccdf_data,ft_gp);
cv_gp = coeffvalues(cf_gp);

fo_wb = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'StartPoint',[0.5 0.5]);
ft_wb = fittype('wblcdf(x,a,b,''upper'')','options',fo_wb);
[cf_wb,gof_wb] = fit(X,ccdf_data,ft_wb);
cv_wb = coeffvalues(cf_wb);

%==CALCULATE OTHER CCDFs==%
tau = mean(times);
ccdf_exp = exp(-X./tau);

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
dataMod = times(~ismember(times,Xrem));
test_data = sort(dataMod)';

z_ml = ones(length(test_data),1)-mlf(beta,1,-gamma*test_data.^beta,6);
z_gp = gpcdf(test_data,k,sigma,theta);
z_wb = wblcdf(test_data,a,b);

stats_ml = testStatistics(test_data,z_ml);
stats_gp = testStatistics(test_data,z_gp);
stats_wb = testStatistics(test_data,z_wb);

stats_ml.Chi_Squared = sum(rdivide((ccdf_data-ccdf_ml).^2,ccdf_ml));
stats_gp.Chi_Squared = sum(rdivide((ccdf_data-ccdf_gp).^2,ccdf_gp));
stats_wb.Chi_Squared = sum(rdivide((ccdf_data-ccdf_wb).^2,ccdf_wb));

stats_ml.Root_MSE = gof_ml.rmse;
stats_gp.Root_MSE = gof_gp.rmse;
stats_wb.Root_MSE = gof_wb.rmse;
stats_ml.R_Squared = gof_ml.rsquare;
stats_gp.R_Squared = gof_gp.rsquare;
stats_wb.R_Squared = gof_wb.rsquare;

%==PLOT CCDFs==%
CCDF_fig = figure();
hold on
plot(X,ccdf_data,'o')
plot(X,ccdf_exp)
plot(X,ccdf_ml)
plot(X,ccdf_gp)
plot(X,ccdf_wb)
set(gca,'XScale','log');
set(gca,'YScale','log');
axis([1E1,1E4,1E-5,1E0]);
xlabel('Contact Time (s)');
ylabel('CCDF');
legend('Data','Exp','ML','Gen. Pareto','Weibull');
hold off
imagefilename = [dir_ref,'/plotMLcomparison-CCDF_fig-img.png'];
print(imagefilename,'-dpng')
close(CCDF_fig);

%===AP Stuff===%
[F_AP,X_AP] = ecdf(activityPot);
ap_ccdf = 1-F_AP;

fo_ex = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[1]);
ft_ex = fittype('expcdf(x,lambda,''upper'')','options',fo_ex);
[cf_ex,gof_ex] = fit(X_AP,ap_ccdf,ft_ex);
cv_ex = coeffvalues(cf_ex);

fo_gm = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[inf inf],'StartPoint',[1 1]);
ft_gm = fittype('gamcdf(x,a,b,''upper'')','options',fo_gm);
[cf_gm,gof_gm] = fit(X_AP,ap_ccdf,ft_gm);
cv_gm = coeffvalues(cf_gm);

fo_rl = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[1]);
ft_rl = fittype('raylcdf(x,sigma,''upper'')','options',fo_rl);
[cf_rl,gof_rl] = fit(X_AP,ap_ccdf,ft_rl);
cv_rl = coeffvalues(cf_rl);

%==BEST FIT CCDFs==%
lambda1 = cv_ex(1);
ccdf_ex = expcdf(X_AP,lambda1,'upper');

a1 = cv_gm(1);
b1 = cv_gm(2);
ccdf_gm = gamcdf(X_AP,a1,b1,'upper');

sigma1 = cv_rl(1);
ccdf_rl = raylcdf(X_AP,sigma1,'upper');

%==GoF Testing==%
test_data = sort(activityPot)';

z_ex = expcdf(test_data,lambda1);
z_gm = gamcdf(test_data,a1,b1);
z_rl = raylcdf(test_data,sigma1);

stats_ex = testStatistics(test_data,z_ex);
stats_gm = testStatistics(test_data,z_gm);
stats_rl = testStatistics(test_data,z_rl);

stats_ex.Root_MSE = gof_ex.rmse;
stats_gm.Root_MSE = gof_gm.rmse;
stats_rl.Root_MSE = gof_rl.rmse;
stats_ex.R_Squared = gof_ex.rsquare;
stats_gm.R_Squared = gof_gm.rsquare;
stats_rl.R_Squared = gof_rl.rsquare;


AP_fig = figure();
hold on
plot(X_AP,ap_ccdf,'o');
plot(X_AP,ccdf_ex);
plot(X_AP,ccdf_gm);
plot(X_AP,ccdf_rl);
xlabel('Activity Potential');
ylabel('CCDF');
legend('Data','Exp','Gamma','Rayleigh');
hold off
apfilename = [dir_ref,'/plotMLcomparison-AP_fig-img.png'];
print(apfilename,'-dpng')
close(AP_fig);

%==Save Relevant Data==%
struc_ml = struct('Stability',beta,'Scale',gamma);
struc_gp = struct('Shape',k,'Scale',sigma,'Location',theta);
struc_wb = struct('Scale',a,'Shape',b);
struc_ex = struct('Rate',lambda1);
struc_gm = struct('Shape',a1,'Scale',b1);
struc_rl = struct('Scale',sigma1);

EX = struct('Parameters',struc_ex,'Statistics',stats_ml);
GM = struct('Parameters',struc_gm,'Statistics',stats_gm);
RL = struct('Parameters',struc_rl,'Statistics',stats_rl);
ML = struct('Parameters',struc_ml,'Statistics',stats_ml);
GP = struct('Parameters',struc_gp,'Statistics',stats_gp);
WB = struct('Parameters',struc_wb,'Statistics',stats_wb);

datafilename = [dir_ref,'/plotMLcomparison-data.mat'];
save(datafilename,'ML','GP','WB','EX','GM','RL')

degreeleveldists(sorteddata,max_time,dir_ref,times);
create_avi(data,dir_ref);