function [] = degreeleveldists(data,max_time,startTime)

all_people = [data(:,2)' data(:,3)'];
uni_people = unique(all_people);
tot_people = size(uni_people,2);
timesteps = max_time/20;

X = linspace(20,max_time,timesteps);
CCDF = zeros(tot_people,timesteps);

parfor i=1:tot_people
    current_id = uni_people(i);
    idx = (data(:,2)==current_id|data(:,3)==current_id);
    current_data = data(idx,:);
    current_numberrows = size(current_data,1);
    times = zeros(1,current_numberrows);
    
    j = 1;
    q = 1;
    while j<current_numberrows+1
        contact_time = 20;
        step_vector = [20 0 0];
        current_row = current_data(j,:);
        if j == current_numberrows
            next_row = [0 0 0];
        else
            next_row = current_data(j+1,:);
        end
        while isequal(next_row,current_row+step_vector)
            contact_time = contact_time+20;
            j = j+1;
            current_row = current_data(j,:);
            if j == current_numberrows
                next_row = [0 0 0];
            else
                next_row = current_data(j+1,:);
            end
        end
        times(q) = contact_time;
        j = j+1;
        q = q+1;
    end

    times(~times) = [];
    
    numbertimes = size(times,2);
    for l=1:timesteps
        this_time = X(l);
        CCDF(i,l) = 1-(sum(times <= this_time)/numbertimes);
    end
end

X = X(1:end-10);
CCDF = CCDF(:,(1:end-10));
Y_average = mean(CCDF);

%==BEST FIT CCDFs==%
fo_ml = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[1 1],'StartPoint',[0.5 0.5]);
ft_ml = fittype('mlf(beta,1,-gamma*x.^beta,6)','options',fo_ml);
cf_ml = fit(X',Y_average',ft_ml);
cv_ml = coeffvalues(cf_ml);

fo_gp = fitoptions('Method', 'NonlinearLeastSquares','Lower',[-inf -inf 0],'StartPoint',[0.5 0.5 0.5]);
ft_gp = fittype('gpcdf(x,k,sigma,theta,''upper'')','options',fo_gp);
cf_gp = fit(X',Y_average',ft_gp);
cv_gp = coeffvalues(cf_gp);

fo_wb = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'StartPoint',[0.5 0.5]);
ft_wb = fittype('wblcdf(x,a,b,''upper'')','options',fo_wb);
cf_wb = fit(X',Y_average',ft_wb);
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

%==PLOT AND SAVE==%
figure()
hold on
plot(X,CCDF,'o','MarkerEdgeColor',[0.8 0.8 0.8])
p1 = plot(X,Y_average,'o');
c1 = plot(X,ccdf_ml);
c2 = plot(X,ccdf_gp);
c3 = plot(X,ccdf_wb);
set(gca,'XScale','log');
set(gca,'YScale','log');
axis([1E1,1E4,1E-5,1E0]);
xlabel('Contact Time (s)');
ylabel('CCDF');
legend([p1 c1 c2 c3],{'Average at Deg. Level','ML','Gen. Pareto','Weibull'});
hold off

imagefilename = [startTime,'/degreeleveldists-img.png'];
print(imagefilename,'-dpng')
close

datafilename = [startTime,'/degreeleveldists-data.mat'];
save(datafilename)

RN = randi([1 tot_people],1,1);
randCCDF = CCDF(RN,:);
randomDLD(startTime,X,randCCDF);