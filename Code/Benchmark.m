%clear all; clc; close all; 

tic
log_n = [4:0.1:7];
n = [1:10000 10.^log_n]; %Number of simulations : 
n = ceil(n);

%Benchmark for Code: 
%1. radius (r)                                  :   Normal mu = 0, sigma = 1
%2. transmissivity (T)                          :   Exponential beta = 1
[a,b] = size(n);
mean_r = zeros(1,b); var_r = zeros(1,b);
momen3_r = zeros(1,b);
momen4_r = zeros(1,b);

mean_T = zeros(1,b); var_T = zeros(1,b);
momen3_T = zeros(1,b);
momen4_T = zeros(1,b);

for m = 1:b
    r = zeros(1,n(m)); T = r;
    
    fprintf('n = %d\n',n(m));
    %Initialization and Calculation
    for i = 1:n(m)
        r(1,i) = random('Normal', 0, 1);
        T(1,i) = random('Exponential', 1);
    end
    
   
    %Post processing : 
    mean_r(m) = mean(r);
    deviation = r - mean_r(m);
    var_r(m) = mean(deviation.^2);
    momen3_r(m) = mean(deviation.^3);
    momen4_r(m) = mean(deviation.^4);
    
    mean_T(m) = mean(T);
    deviation = T - mean_T(m);
    var_T(m) = mean(deviation.^2);
    momen3_T(m) = mean(deviation.^3);
    momen4_T(m) = mean(deviation.^4);
    
end
%Normal : SKewness = 0, Kurtosis = 3
skewness_r = momen3_r./(sqrt(var_r).^3);
kurtosis_r = momen4_r./(sqrt(var_r).^4);

figure(b+1)
plot(log(n)/log(10),mean_r)
title('Mean for Normal distribution');

figure(b+2)
plot(log(n)/log(10),var_r)
title('Variance for Normal distribution');

figure(b+3)
plot(log(n)/log(10),skewness_r)
title('Skewness for Normal distribution');

figure(b+4)
plot(log(n)/log(10),kurtosis_r)
title('Kurtosis for Normal distribution');

%Exponential : Skewness = 2, 
skewness_T = momen3_T./(sqrt(var_T).^3);
kurtosis_T = momen4_T./(sqrt(var_T).^4);

figure(b+5)
plot(log(n)/log(10),mean_T)
title('Mean for Exponential distribution');

figure(b+6)
plot(log(n)/log(10),var_T)
title('Variance for Exponential distribution');

figure(b+7)
plot(log(n)/log(10),skewness_T)
title('Skewness for Exponential distribution');

figure(b+8)
plot(log(n)/log(10),kurtosis_T)
title('Kurtosis for Exponential distribution');


disp('Done')
toc

