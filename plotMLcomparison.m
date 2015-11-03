function [] = plotMLcomparison

%==CALCULATE CCDF FOR SAMPLE==%
input = 'journal.pone.0136497.s002.CSV';
fid = fopen(input);
rawdata = textscan(fid,'%f %f %f %*s %*s','Delimiter',',');
fclose(fid);

data = cell2mat(rawdata);
data(:,1) = data(:,1)-data(1,1);

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
while j<number_rows+1
    contact_time = 20;
    step_vector = [20 0 0];
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

[F,X] = ecdf(times);
ccdf_data = 1-F;

%==MODIFY TO REMOVE ERRORS==%
X = X(2:end-10);
ccdf_data = ccdf_data(2:end-10);

%==BEST FIT CCDFs==%
fo_ml = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[1 1],'StartPoint',[0.5 0.5]);
ft_ml = fittype('mlf(beta,1,-gamma*x.^beta,6)','options',fo_ml);
cf_ml = fit(X,ccdf_data,ft_ml);
cv_ml = coeffvalues(cf_ml);

ft_gp = fittype('gpcdf(x,xi,sigma,theta,''upper'')');
cf_gp = fit(X,ccdf_data,ft_gp);
cv_gp = coeffvalues(cf_gp);

%==CALCULATE OTHER CCDFs==%
tau = mean(times);
ccdf_exp = exp(-X./tau);
beta = cv_ml(1);
gamma = cv_ml(2);
ccdf_ml = mlf(beta,1,-gamma*X.^beta,6);
sigma = cv_gp(1);
theta = cv_gp(2);
xi = cv_gp(3);
ccdf_gp = gpcdf(X,xi,sigma,theta,'upper');

%==PLOT CCDFs==%
figure()
hold on
plot(X,ccdf_data,'o')
plot(X,ccdf_exp)
plot(X,ccdf_ml)
plot(X,ccdf_gp)
set(gca,'XScale','log');
set(gca,'YScale','log');
axis([1E1,1E4,1E-5,1E0]);
xlabel('Contact Time (s)');
ylabel('CCDF');
legend('Data','Exp','ML','Gen. Pareto')