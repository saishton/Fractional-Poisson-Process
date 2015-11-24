function [] = plotMLcomparison(input,structure)

% input = 'journal.pone.0136497.s002.CSV';
% structure = '%f %f %f %*s %*s';

startTime = datestr(now,'yyyymmddTHHMMSS');
mkdir(startTime);

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

times(~times) = [];
activityPot = interactions/total_interactions;

[F,X] = ecdf(times);
ccdf_data = 1-F;
max_time = X(end);

%==MODIFY TO REMOVE ERRORS==%
X = X(2:end-10);
ccdf_data = ccdf_data(2:end-10);

%==BEST FIT CCDFs==%
fo_ml = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[1 1],'StartPoint',[0.5 0.5]);
ft_ml = fittype('mlf(beta,1,-gamma*x.^beta,6)','options',fo_ml);
cf_ml = fit(X,ccdf_data,ft_ml);
cv_ml = coeffvalues(cf_ml);

fo_gp = fitoptions('Method', 'NonlinearLeastSquares','Lower',[-inf -inf 0],'StartPoint',[0.5 0.5 0.5]);
ft_gp = fittype('gpcdf(x,k,sigma,theta,''upper'')','options',fo_gp);
cf_gp = fit(X,ccdf_data,ft_gp);
cv_gp = coeffvalues(cf_gp);

fo_wb = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'StartPoint',[0.5 0.5]);
ft_wb = fittype('wblcdf(x,a,b,''upper'')','options',fo_wb);
cf_wb = fit(X,ccdf_data,ft_wb);
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
imagefilename = [startTime,'/plotMLcomparison-CCDF_fig-img.png'];
print(imagefilename,'-dpng')
close(CCDF_fig);

AP_fig = figure();
histogram(activityPot,'Normalization','cdf');
xlabel('Activity Potential');
ylabel('Cumulative Density');
apfilename = [startTime,'/plotMLcomparison-AP_fig-img.png'];
print(apfilename,'-dpng')
close(AP_fig);

datafilename = [startTime,'/plotMLcomparison-data.mat'];
save(datafilename)

degreeleveldists(sorteddata,max_time,startTime);
create_avi(data,startTime);