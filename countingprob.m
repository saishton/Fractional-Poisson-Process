% The distribution of the fractional Poisson process

clear all
close all

beta=0.99; %index
nu=1; %skewness
% u=5;
% gamma = (u*cos(pi*beta/2))^(1/beta); %scale parameter
delta=0; %location parameter

t=1; %time

% Monte Carlo simulation

gamma_t = 1;

N=10000;

num=zeros(1,N);

for k=1:N
k
dt=mlrnd(beta,gamma_t,1,100000);

time=cumsum(dt);

index=find(time>t);

num(k)=min(index)-1;

end

%building the histogram

for j=1:max(num)
    freq(j)=length(find(num==j-1))/N;
end

% Analytical formula

prob0 = mlf(beta,1,-t^beta,5); %probability for n=0

for n=1:(max(num)-1) %number of events
    n
% 
% stable = stblcdf(t,beta,nu,gamma,delta);
% 
% plot(t,stable)

% fun = @(u) stblcdf(t,beta,nu,(u.*cos(pi*beta/2)).^(1/beta),delta)*(exp(-u).*u.^(n-1)/factorial(n-1) - exp(-u).*u.^n/factorial(n));
% prob = integral(fun,0.0001,Inf);

du=0.01;
int=0;

% for i=0:100000
%     int= int + du*stblcdf(t,beta,nu,(i*du*cos(pi*beta/2))^(1/beta),delta)*(exp(-i*du)*(i*du)^(n-1)/factorial(n-1) - exp(-i*du)*(i*du)^n/factorial(n));
% end

for i=0:100000
    int= int + du*stblcdf(t,beta,nu,(i*du*cos(pi*beta/2))^(1/beta),delta)*exp(-i*du)*(i*du)^(n-1)/factorial(n-1)*((n-i*du)/n);
end

prob(n)=int;
probpoiss(n) = poisspdf(n,t);
end

prob=[prob0 prob];save
probpoiss=[poisspdf(0,t) probpoiss]; %Poisson for comparison

% Plots
x=0:(max(num)-1); %x axis
plot(x,prob,'o')
hold on
plot(x,probpoiss,'or')
plot(x,freq,'x')

% for i=1:1000
% funct(i) = stblcdf(t,beta,nu,(i*du*cos(pi*beta/2))^(1/beta),delta)*(exp(-i*du)*(i*du)^(n-1)/factorial(n-1) - exp(-i*du)*(i*du)^n/factorial(n));
% end
% 
% for i=1:1000
%     funct2(i) = stblcdf(t,beta,nu,(i*du*cos(pi*beta/2))^(1/beta),delta);
% end
% 
% for i=1:1000
%     funct3(i) = (exp(-i*du)*(i*du)^(n-1)/factorial(n-1) - exp(-i*du)*(i*du)^n/factorial(n));
% end

% funct3 is the difference between two Poisson probabilities (lambda = 1)
% 
% for i=1:1000
%     funct4(i) = poisspdf(n-1,i*du) - poisspdf(n,i*du);
% end