clear all; clc; close all; 
tic
log_n = [4:0.1:7];
n = [1:1:10^4 10.^log_n];
n = ceil(n);
%n = 1000:1000:10000000; %Number of simulations : 

%Borehole function Input in array: 
%1. yield stress (R)                            :   Normal mu = 40000, sigma = 2000
%2. Young's modulus (E)                         :   Normal mu = 2.9*10^7, sigma = 1.45*10^6
%3. Horizontal load (X)                         :   Normal mu = 500, sigma = 100
%4. Vertical load (Y)                           :   Normal mu = 1000, sigma = 100
[a,b] = size(n);
%Displacement
mean_D = zeros(1,b);
var_D = zeros(1,b);
momen3_D = zeros(1,b);
momen4_D = zeros(1,b);

%Stress
mean_S = zeros(1,b);
var_S = zeros(1,b);
momen3_S = zeros(1,b);
momen4_S = zeros(1,b);

L = 100; D0 = 2.2535; w = 4; t = 2;
for m = 1:b
    R = zeros(1,n(m)); E = R; X = R; Y = R; D = zeros(1,n(m)); S = D;
    
    fprintf('n = %d\n',n(m));
    %Initialization and Calculation
    for i = 1:n(m)
        R(1,i) = random('Normal', 40000, 2000);
        E(1,i) = random('Normal', 2.9*10^7, 1.45*10^6);
        X(1,i) = random('Normal', 500, 100);
        Y(1,i) = random('Normal', 1000, 100);
        
        S(1,i) = 600*Y(1,i)/(w*t^2) + 600*X(1,i)/(w^2*t);
        D(1,i) = 4*L^3/(E(1,i)*w*t) *sqrt( (Y(1,i)/t^2)^2 + (X(1,i)/w^2)^2);
    end
    
   
    %Post processing : 
    mean_D(m) = mean(D);
    deviation2 = D - mean_D(m);
    var_D(m) = mean(deviation2.^2);
    momen3_D(m) = mean(deviation2.^3);
    momen4_D(m) = mean(deviation2.^4);
    
    mean_S(m) = mean(S);
    deviation2 = S - mean_S(m);
    var_S(m) = mean(deviation2.^2);
    momen3_S(m) = mean(deviation2.^3);
    momen4_S(m) = mean(deviation2.^4);
    
end
skewness_D = momen3_D./(sqrt(var_D).^3);
kurtosis_D = momen4_D./(sqrt(var_D).^4);

skewness_S = momen3_D./(sqrt(var_S).^3);
kurtosis_S = momen4_D./(sqrt(var_S).^4);

figure(b+1)
plot(log(n)/log(10),mean_D)

figure(b+2)
plot(log(n)/log(10),var_D)

figure(b+3)
plot(log(n)/log(10),skewness_D)

figure(b+4)
plot(log(n)/log(10),kurtosis_D)

figure(b+5)
plot(log(n)/log(10),mean_S)

figure(b+6)
plot(log(n)/log(10),var_S)

figure(b+7)
plot(log(n)/log(10),skewness_S)

figure(b+8)
plot(log(n)/log(10),kurtosis_S)

disp('Done')
toc

%{ 
Determing distribution function :

figure(1)
histogram(S)
histfit(S,10^4,'normal')
pd_Normal_S = fitdist(transpose(S),'normal')
figure(2)
histogram(S)
histfit(S,10^4,'lognormal')
pd_Lognormal_S = fitdist(transpose(S),'lognormal')
figure(3)
histogram(S)
histfit(S,10^4,'Weibull')
pd_Weibull_S = fitdist(transpose(S),'weibull')
figure(4)
histogram(S)
pd_Gamma_S = fitdist(transpose(S),'gamma')
histfit(S,10^4,'Gamma')
title('Gamma distribution fitting for S')
figure(3)
title('Weibull distribution fitting for S')
figure(2)
title('Lognormal distribution fitting for S')
figure(1)
title('Normal distribution fitting for S')
%}