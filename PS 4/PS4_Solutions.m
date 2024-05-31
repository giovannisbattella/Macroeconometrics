clear; clc; close all;
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
rng(12345);

%% Load the data

% Fernald:
data = readtable('quarterly_tfp.xlsx', 'Sheet', 'quarterly');
hours = data.dhours(2:308);  % change in hours
prod = data.dLP(2:308);   % change in labor productivity
dates = datetime("1947-04-01"):calmonths(3):datetime("2023-10-01");

% FRED:
% data1 = readtable('techdata2.csv');
% hours_pc = diff(log(data1.HOANBS./data1.CNP16OV))*100;
% hours = diff(log(data1.HOANBS))*100;
% prod = diff(log(data1.OPHNFB))*100;
% dates = data1.DATE(2:end);

y = [prod, hours]; % Data from Q2:1947 to Q4:2023

figure;
subplot(2,1,1)
plot(dates, prod, 'LineWidth', 1.5)
recessionplot
axis tight
title('Labor productivity growth', 'FontSize', 14, 'Interpreter','Latex')

subplot(2,1,2)
plot(dates, hours, 'LineWidth', 1.5)
recessionplot
axis tight
title('Hours worked growth', 'FontSize', 14, 'Interpreter','Latex')


%% Estimate the VAR
% Determine the lag length for the VAR using AIC
[paic, pbic, phq] = aicbic(y,8,1);
p = 4;
c = 1;

% Estimate the VAR yt = c + A1y_t-1 + ... Apyt-p + et
[beta, residuals] = VAR(y,p,c);

% Compute the variance-covariance matrix of the reduced form errors et
[T, N] = size(y);
Sigma = (residuals' * residuals)./(T-1-p-N*p);

% Compute the Wold IRFs by inverting (I-A(L))
horizon = 40;
wold = woldirf(beta, c, p, horizon);

%% Compute the LR responses and identify technology shocks
% Gali assumes that only technology shocks can
% have a permanent effect on the level of technology. Second shock is the
% non-technology shock then.

C1 = sum(wold,3);
D1 = chol(C1*Sigma*C1', 'lower');
K = C1\D1;
% disp(K)

% Compute the LR responses
irf_bq = bqIRF(wold, K);

%% Bootstrap confidence bands
nboot1 = 1000;
nboot2 = 2000;
prc = 68;
cumulate = [1,2];
[bootbq, upper, lower, boot_beta] = ...
    bootstrapBQ_corrected(y,p,c,beta,residuals,nboot1,horizon,prc,cumulate,nboot2);

%% Plot the IRFs
varnames = {'LABPROD', 'HOURS'};
shocknames = {'Technology', 'Non-Technology'} ;
shocks = [1,2];
plotLR(irf_bq, varnames, shocknames, cumulate, shocks, upper, lower, prc)

%% Compute the FEVD
irf_cumu = cumsum(irf_bq,3);
fevd = variance_decomp(irf_cumu);

fevd(:,:,[1,5,9,21,41])

plot_vardec(fevd, varnames, shocknames)

%% Compute the historical decomposition
histedec_all = zeros(T-p,N,N);
ystar_all = zeros(T-p,N);

for series=1:N
    [histdec, ystar] = hist_decmp(y, beta, residuals, c, p, K, series);
    histedec_all(:,:,series) = histdec;
    ystar_all(:,series) = ystar;
end

%% Plot HD
dates_final = dates(p+1:end);

series = 1:1:N;
varnames = varnames(series);
plothistdec(histedec_all,ystar_all,shocks,series,dates_final,shocknames,varnames)

%% Compute the conditional correlations of Gali (1999)
[N,~,horizon] = size(irf_bq);

% Compute the numerator for each shock (the product of the MA coefficients)
numerator = zeros(N,1);
for i=1:N
    
    temp = irf_bq(1,i,:) .* irf_bq(2,i,:);
    numerator(i) = sum(temp,3);
    
end

% Compute the denominators 
irf_squared = irf_bq.^2;

denominator = zeros(N,1);
for i=1:N
    v1 = sum(irf_squared(1,i,:),3);
    v2 = sum(irf_squared(2,i,:),3);
    denominator(i) = sqrt(v1 * v2);  
end

% Conditional correlations
cond_corr = numerator ./ denominator;

% Unconditional correlations
uncond_corr = corr(y);

% Just like Gali (1999) we find that the unconditional correlation between
% hours worked and productivity is small which is at odds with standard RBC
% models. However, the conditional correlation, i.e. the part of the
% correlation between hours and productivity driven by technology is
% negative. I.e. when technology increases, hours should decrease. On the
% other hand, the correlation for the non-technology shock (in Gali's
% example monetary shocks) is positive, which RBC predicted for the
% technology part. Of course this criticism of RBC models was met with
% critiques that the shock identified is not a technology shock.

%% Repeat the exercise but with AS-AD shocks: idea is to identify an AS shock as 
% the only driver of LR GDP
clear; clc; close all;

data = readtable("as_ad_data.csv");
gdp = diff(log(data.GDPC1))*100;
infl = diff(log(data.GDPDEF))*100;
dates = data.DATE(2:end);

y = [gdp, infl];

% Estimate the VAR
% Determine the lag length for the VAR using AIC
p = aicbic(y,4,1);
c = 1;

% Estimate the VAR yt = c + A1y_t-1 + ... Apyt-p + et
[beta, residuals] = VAR(y,p,c);

% Compute the variance-covariance matrix of the reduced form errors et
[T, N] = size(y);
Sigma = (residuals' * residuals)./(T-1-p-N*p);

% Compute the Wold IRFs by inverting (I-A(L))
horizon = 40;
wold = woldirf(beta, c, p, horizon);

% Compute the LR responses and identify technology shocks
% Gali assumes that only technology shocks can
% have a permanent effect on the level of technology. Second shock is the
% non-technology shock then.

C1 = sum(wold,3);
D1 = chol(C1*Sigma*C1', 'lower');
K = C1\D1;
% disp(K)

% Compute the LR responses
irf_bq = bqIRF(wold, K);

% Bootstrap confidence bands
nboot1 = 1000;
nboot2 = 2000;
prc = 68;
cumulate = [1,2];
[bootbq, upper, lower, boot_beta] = ...
    bootstrapBQ_corrected(y,p,c,beta,residuals,nboot1,horizon,prc,cumulate,nboot2);

% Plot the IRFs
varnames = {'GDP', 'PRICES'};
shocknames = {'AS', 'AD'} ;
shocks = [1,2];
plotLR(irf_bq, varnames, shocknames, cumulate, shocks, upper, lower, prc)

% Compute the FEVD
irf_cumu = cumsum(irf_bq,3);
fevd = variance_decomp(irf_cumu);
plot_vardec(fevd, varnames, shocknames)

% Compute the historical decomposition
histedec_all = zeros(T-p,N,N);
ystar_all = zeros(T-p,N);

for series=1:N
    [histdec, ystar] = hist_decmp(y, beta, residuals, c, p, K, series);
    histedec_all(:,:,series) = histdec;
    ystar_all(:,series) = ystar;
end

% Plot HD
dates_final = dates(p+1:end);

series = 1:1:N;
varnames = varnames(series);
plothistdec(histedec_all,ystar_all,shocks,series,dates_final,shocknames,varnames)

