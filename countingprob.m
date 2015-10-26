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

format long g

runStartTime = datestr(now,'yyyymmddTHHMMSS');

nu=1; %Skewness parameter
delta=0; %Location parameter

%===Run a Monte Carlo simulation and build histogram===%

gamma_t = 1;
N=10000;
num=zeros(1,N);

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

prob0 = mlf(beta,1,-t^beta,5); %Probability for n=0
prob = zeros(1,BigN-1); %Preallocate vector size for speed
probpoiss = zeros(1,BigN-1); %Preallocate vector size for speed

for n=1:(BigN-1) 
    int=0;
    intSteps=round(1000/du)+1;
    parfor i=0:intSteps
        int=int+du*stblcdf(t,beta,nu,(i*du*cos(pi*beta/2))^(1/beta),delta)*...
            exp(-i*du)*(i*du)^(n-1)/factorial(n-1)*((n-i*du)/n);
    end
    prob(n)=int;
    probpoiss(n) = poisspdf(n,t);
end

prob=[prob0 prob];
probpoiss=[poisspdf(0,t) probpoiss]; %Poisson for comparison

%===Plot data & save to file===%
x=0:(BigN-1);
plot(x,prob,'o')
hold on
plot(x,probpoiss,'or')
plot(x,freq,'x')

filename_image = ['graphical-',runStartTime,'.png'];
filename_vars = ['variables-',runStartTime,'.mat'];
print(filename_image,'-dpng')
save(filename_vars)
end
