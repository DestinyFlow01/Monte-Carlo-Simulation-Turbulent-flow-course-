clear all; clc; close all; 

tic
n = [1:1:10^3]; %Number of simulations : 

%Benchmark for Code: 
%1. radius (r)                                  :   Normal mu = 0.1, sigma = 0.0161812
%2. transmissivity (T)                          :   Exponential mu = 0.1, sigma = 0.0161812
%3. potentiometric head of upper aquifer (Hu)   :   Uniform 990 to 1110
%4. transmissivity of lower aquifer (Tl)        :   Uniform 63.1 to 116
%5. potentiometric head of lower aquifer (Hl)   :   Uniform 700 to 820
%6. length of borehole (L)                      :   Uniform 1120 to 1680
%7. hydraulic conductivity of borehole (Kw)     :   Uniform 9855 to 12045
[a,b] = size(n);
mean_f = zeros(1,b);
var_f = zeros(1,b);
momen3_f = zeros(1,b);
momen4_f = zeros(1,b);

for m = 1:b
    rw = zeros(1,n(m)); r = rw; Tu = rw; Hu = rw; Tl = rw; Hl = rw; L = rw; Kw = rw; f = zeros(1,n(m));
    
    fprintf('n = %d\n',n(m));
    %Initialization and Calculation
    for i = 1:n(m)
        rw(1,i) = random('Normal', 0.1, 0.0161812);
        r(1,i) = random('Lognormal', 7.71, 1.0056);
        Tu(1,i) = random('Uniform', 63070, 115600);
        Hu(1,i) = random('Uniform', 990, 1110);
        Tl(1,i) = random('Uniform', 63.1, 116);
        Hl(1,i) = random('Uniform', 700, 820);
        L(1,i) = random('Uniform', 1120, 1680);
        Kw(1,i) = random('Uniform', 9855, 12045);
        numerator = 2*pi*Tu(1,i)*(Hu(1,i) - Hl(1,i));
        denominator = log(r(1,i)/rw(1,i))*( 1 + Tu(1,i)/Tl(1,i) + 2*L(1,i)*Tu(1,i)/(Kw(1,i)*rw(1,i)^2 * log(r(1,i)/rw(1,i))));
        f(1,i) = numerator/denominator;
    end
    
   
    %Post processing : 
    mean_f(m) = mean(f);
    deviation = f - mean_f(m);
    var_f(m) = mean(deviation.^2);
    momen3_f(m) = mean(deviation.^3);
    momen4_f(m) = mean(deviation.^4);
    
end
skewness_f = momen3_f./(sqrt(var_f).^3);
kurtosis_f = momen4_f./(sqrt(var_f).^4);

figure(b+1)
plot(log(n)/log(10),mean_f)

figure(b+2)
plot(log(n)/log(10),var_f)

figure(b+3)
plot(log(n)/log(10),skewness_f)

figure(b+4)
plot(log(n)/log(10),kurtosis_f)

disp('Done')
toc