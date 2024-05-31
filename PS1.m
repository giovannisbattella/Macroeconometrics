%Exercise 1
datab=xlsread("PS1.xls", 'FRED Graph', 'C14:C320');
data = datab/100;
plot(data(:,1))
title('Percentage change of GDP, US, Quarterly data');
xlabel('Periods');
ylabel('GDP growth');
%% Exercise 2_OLS to estimate C, phi and sigma squared
numLags = 1
laggeddata = [NaN(numLags,1); data(1:end-numLags)];
x = laggeddata;
y = data;

model = fitlm (x,y)
%find the variance
variance_of_residuals = model.MSE;
disp(['Variance of the residuals: ', num2str(variance_of_residuals)]);

%% Exercise 3_autocorrelation
autocorr(data, 10);
title('Autocorrelation Function (ACF) 10 Lags');
xlabel('Lag');
ylabel('Autocorrelation');
%% Exercise 4_obtain and plot the coefficient of the Wold representation of y
p = 1; %number of lags
b = AREstimation(y,p);
F_AR1=[b(2:end)'; eye(p-1) zeros(p-1,1)]; % Since the process is an AR1 the F matrix
% is just the value of the coefficient of the first lag of the AR process.

% Let´s represent an MA process of 20 periods.
J=10

for j=1:J
    MATemp=F_AR1^(j-1);
    MACoeff_AR1(j)=MATemp(1,1);
end 
   
figure (5)  
plot(MACoeff_AR1, '*-');
title('Coefficients of the Wold representation of 10 periods-AR1');
%% Exercise 5_AR(2)
numlags = 2
laggeddata2 = [NaN(numlags,1); data(1:end-numlags)];
X = [laggeddata laggeddata2];

model2 = fitlm (X,y)
%find the variance
variance_of_residuals_AR(2) = model2.MSE;
disp(['Variance of the residuals AR(2): ', num2str(variance_of_residuals)]);


%% Exercise 6
%find the root

psi1 = model2.Coefficients.Estimate(2);
psi2 = model2.Coefficients.Estimate(3);

rootsMA2 = abs(roots([psi2 psi1 1]));

if rootsMA2>1 
    disp('The process is causal since the absolute value of the roots is bigger than 1')
else
    disp('The process is not causal')
end

%is causal and stationary



%% Ex 8_wold coefficients
p = 2;
c = AREstimation(y,p);
F_AR2 = [c(2:end)';eye(p-1) zeros(p-1,1)];

for j=1:J
    MATemp = F_AR2^(j-1)
    MACoeff_AR2(j) = MATemp(1,1)
end

figure (6)  
plot(MACoeff_AR2, '*-');
title('Coefficients of the Wold representation of 10 periods-AR2');
%% Exercise 10

% We can check the roots of the characteristic MA-polynomial:
% The process is yt = eps_t + psi1 eps_t-1 + psi2 eps_t-2. Writing this in
% terms of the lag operator we obtain: yt = [1 + psi1*L + psi2*L^2]eps_t.
% Thus, the characteristic polynomial is psi(z) = (1+psi1*z+psi2*z^2). We
% need the roots of this polynomical to be outside the (complex) unit
% circle. Fundamentalness rules out the existence of z = 1 as a solution.

% The function "roots" compute the roots of a polynomial. We have to feed
% in the coefficients attached to the terms in the polynomial, starting
% with the one attached to the highest polynomial degree in descending
% order. 

% In this case we have psi2 attached to z^2, psi1 attached to z^1 and 1
% attached to z^0.
rootsMAll=abs(roots([2 1.2 1]));


if rootsMAll<1 % This way, the condition returns true only if all values of the logical vector are 1
    disp('The process is invertible')
else
    disp('The process is not invertible')   
end











