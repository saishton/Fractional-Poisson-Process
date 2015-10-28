function [] = countingprob(t,beta,intstep)
% Plots the distribution of the fractional Poisson process for given time
% (t), index (beta) and integral step (intstep) and compares it to a Monte
% Carlo simulation and a standard Poisson distribution.

if nargin < 2
    error('countingprob:TooFewInputs','Requires at least two input arguments.');
elseif nargin > 3
    error('countingprob:TooManyInputs','Accepts at most three input arguments.');
elseif isempty(intstep)
    du = 0.01;
elseif isscalar(intstep) && intstep>0
    du = intstep;
else
    error('countingprob:BadInput','"du" must be a positive scalar.')
end

h = waitbar(0,'Initialising...');
format long g
startTime = datestr(now,'yyyymmddTHHMMSS');

nu=1; %Skewness parameter
delta=0; %Location parameter

%===Run a Monte Carlo simulation and build histogram===%

gamma_t = 1;
N=10000;
num=zeros(1,N);

waitbar(0,h,'Simulating Data...')
parfor k=1:N
dt=mlrnd(beta,gamma_t,1,100000);
time=cumsum(dt);
num(k)=find(time>t,1)-1;
end

BigN = max(num); %Number of events

parfor j=1:BigN
    freq(j)=length(find(num==j-1))/N;
end

%===Estimate an analytical formula===%
prob=zeros(1,BigN); %Preallocate vector size for speed
probpoiss=zeros(1,BigN); %Preallocate vector size for speed
intSteps=round(1000/du)+1;

waitbar(0,h,'Calculating Analytical Estimates...')

for n=0
    prob(n+1)=mlf(beta,1,-t^beta,5);
    progress=(n+1)/BigN;
    waitbar(progress,h)
end  

for n=1
    int=zeros(1,intSteps+1);
    parfor i=0:intSteps
        stepu=i*du;
        int(i+1)=du*stblcdf(t,beta,nu,(stepu*cos(pi*beta/2))^(1/beta),delta)*...
            exp(-stepu)*(1-stepu);
    end
    prob(n+1)=sum(int);
    progress=(n+1)/BigN;
    waitbar(progress,h)
end

for n=2:(BigN-1)
    int=zeros(1,intSteps+1);
    log_noi=log(du)-log(factorial(n-1))-log(n);
    parfor i=0:intSteps
        stepu=i*du;
        if stepu>1
            logint=log(stblcdf(t,beta,nu,(stepu*cos(pi*beta/2))^(1/beta),delta))-...
                stepu+(log(stepu)*(n-1))+log(n-stepu)+log_noi;
            int(i+1)=real(exp(logint));
        else
            int(i+1)=du*stblcdf(t,beta,nu,(stepu*cos(pi*beta/2))^(1/beta),delta)*...
                exp(-stepu)*(stepu)^(n-1)/factorial(n-1)*((n-stepu)/n);
        end
    end
    prob(n+1)=sum(int);
    progress=(n+1)/BigN;
    waitbar(progress,h)
end

parfor n=0:(BigN-1) %Poisson for comparison
    probpoiss(n+1)=poisspdf(n,t);
end

%===Plot data & save to file===%
waitbar(1,h,'Saving Outputs...')
x=0:(BigN-1);
figure()
hold on
plot(x,prob,'o')
plot(x,probpoiss,'or')
plot(x,freq,'x')
hold off

close(h)
datafilename = [startTime,'-data.mat'];
imagefilename = [startTime,'-img.png'];
save(datafilename)
print(imagefilename,'-dpng')
end