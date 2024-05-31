%% PS5: News Shocks

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
close all; clear; clc;

%% Load the data (goes from Q1:1959 to Q4:2022)

% FRED-QD
fred_qd = readtable('fred_qd.csv');
clc; % this is to remove the warning that the data format for the dates is unusual

% TFP data
fernald_data = readtable('quarterly_tfp.xlsx',Sheet='quarterly');
tfp = fernald_data.dtfp_util/4; % We divide by 4 because the quarterly changes are annualized in the sheet (diff(log(x)))*400 is used

% Create dates (in beginning of quarter notation)
dates_fred_qd = datetime('1959-01-01'):calmonths(3):datetime('2022-10-01');
dates_tfp = datetime('1947-01-01'):calmonths(3):datetime('2022-10-01');

% Adjust the samples
fred_data = [fred_qd.GDPC1(3:end), fred_qd.PCECC96(3:end), fred_qd.S_P500(3:end), fred_qd.GDPCTPI(3:end)];
fred_data = fred_data(find(dates_fred_qd == '1959-01-01'):find(dates_fred_qd == '2022-07-01'),:);

tfp = tfp(find(dates_tfp == '1959-01-01'):find(dates_tfp == '2022-07-01'));

%% 1. Estimate a VAR(4) for the growth rate of stock prices and TFP growth

% Deflate the SP500 series
sp500 = fred_data(:,3)./ fred_data(:,4);

% Turn to growth rates
sp500gr = diff(log(sp500))*100;
tfp_gr = tfp(2:end);
y = [tfp_gr, sp500gr]; % if estimated in growth rates
%y = [cumsum(tfp), log(sp500)*100]; % if estimated in log-levels
[T,N] = size(y);

% Estimate the VAR(4)
p = 4;
c = 1;
[beta, residuals] = VAR(y,p,c);
disp('the maximum eigenvalue is')
disp(max(abs(eig(companionMatrix(beta,c,p)))))

%% 2. Identify a news shock using contemporaneous restrictions
% Get the variance-covariance matrix
sigma = (residuals' * residuals)./(T-c-p-N*p);

% Compute Wold coefficients
horizon = 40;
wold = woldirf(beta,c,p,horizon);

% Do this with a VECM
% r = 1; % One co-integrating relationshiop
% vecm_output = vecm_rr(cumsum(y,2),p,c,r,horizon);
% beta_vecm = vecm_output.par_var;
% wold = vecm_output.wold;
% sigma = cov(vecm_output.eps);

% Compute the Cholesky factor of sigma
S = chol(sigma, 'lower');

% Estimate the Cholesky IRFs
cholirf = choleskyIRF(wold,S);

% Plot Point IRF
varnames = {'TFP', 'SP500'};
shockname = "News Shock";
shock = 2;
cumulate = []; 
plotchol(cholirf, varnames, shockname, cumulate, shock)

%% 3. Plot the IRFs and bootstrapped bands for the log level responses
% Bootstrap the confidence bands
prc = 68;
nboot = 2000;
cumulate = [1, 2]; 
%cumulate = []; % if estimated in log-levels
[bootchol, upper, lower, boot_beta] = ...
    bootstrapChol(y,p,c,beta,residuals,nboot,horizon,prc,cumulate);

% Plot Cholesky IRFs
plotchol(cholirf, varnames, shockname, cumulate, shock, upper, lower, prc)

%% 4. Compute and plot the FEVD
level_irf = cumsum(cholirf,3); % if estimated in growth rates
% level_irf = cholirf;
fevd_level = variance_decomp(level_irf, shock);

% Plot the FEVD
plot_vardec(fevd_level, varnames, shockname)

%% 5. Long-run identification
C1 = sum(wold,3); % If estimated in growth rates
%C1 = squeeze(wold(:,:,end));
K = C1\chol(C1*sigma*C1','lower');

% Compute the LR responses
irf_bq = bqIRF(wold, K);

% Bootstrap confidence bands
level_est = 0;
[bootbq, upper, lower, boot_beta] = ...
    bootstrapBQ(y,p,c,beta,residuals,nboot,horizon,prc,cumulate, level_est);

% Plot the IRFs
shocknames = {'Technology Shock', 'Non-technology Shock'};
plotLR(irf_bq, varnames, shocknames, cumulate, [1,2], upper, lower, prc)

%% Second Approach: fsolve

% Using fsolve for a system of non-linear equations K
options=optimset('Display', 'off', 'TolFun',.0000000001,'MaxIter',100000,'MaxFunEvals',100000); % more iterations

start = [1,-2,3,4];
sol = fsolve(@(k) mysystem(k,sigma,C1), start, options);
K1 = [sol(1), sol(2); sol(3), sol(4)];

%% Third Approach: Using the orthonormal representaiton yt = C(L)SHH'S^{-1}et

D1 = sum(cholirf,3);
start = [1,-2,3,4];
sol1 = fsolve(@(h) mysystem1(h,D1), start, options);
H = [sol1(1), sol1(2); sol1(3), sol1(4)];
K2 = S * H;


%% Compute the correlation between the SR identified news shock and the technology shock
u_SR = S\residuals';
u_LR = K\residuals';

corrs = corrcoef([u_SR(2,:)', u_LR(1,:)']);

disp('the correlation between the news shock and the tech shock is')
disp(corrs(1,2))

figure;
scatter(u_SR(2,:)',u_LR(1,:)');hold on;
lsline; hold off;
xlabel('Cholesky Shock')
ylabel('Long-run shock')

%% 6. Maximization identification

% Create the growth rates of GDP and Consumption
GDP_gr = diff(log(fred_data(:,1)))*100;
CONS_gr = diff(log(fred_data(:,2)))*100;
y = [tfp_gr, sp500gr, GDP_gr, CONS_gr];

% for level:
% y = [cumsum(tfp_gr), cumsum(sp500gr), cumsum(GDP_gr), cumsum(CONS_gr)]; 


[T,N] = size(y);

% Estimate the VAR(4)
[beta, residuals] = VAR(y,p,c);
disp('the maximum eigenvalue is')
disp(max(abs(eig(companionMatrix(beta,c,p)))))

% Get the variance-covariance matrix
sigma = (residuals' * residuals)./(T-c-p-N*p);

% Compute Wold coefficients
wold = woldirf(beta,c,p,horizon);

% Compute the Cholesky factor of sigma
S = chol(sigma, 'lower');

% Estimate the Cholesky IRFs
cholirf = choleskyIRF(wold,S);

%% Identification using maximum of LR response

% Using fminsearch
var = 1; %TFP is ordered first
startvalue=ones(3,1)./norm(ones(3,1));
h_est = fminsearch(@(h) lr_target(h,cholirf,var),startvalue,options);
% h_est = fminsearch(@(h) lr_target_level(h,cholirf,var),startvalue,options); % for level

h = [0;h_est./norm(h_est)];

%% Second approach: fmincon
h_est = fmincon(@(h) lr_target1(h,cholirf,var),startvalue,[],[],[],[],[],[],@constraint);
h1 = [0;h_est./norm(h_est)];

%% Third approach: Givens rotations
% Using the rotation matrices

startvalue=zeros(2,1);
theta_est=fminunc(@(theta) -h_vec_max(theta,cholirf),startvalue,options);

theta_13=theta_est(1);
theta_14=theta_est(2);

h2=[0;-cos(theta_13)*cos(theta_14); -sin(theta_13)*cos(theta_14); -sin(theta_14)];

%% Compute the point IRF
max_irf = partial_irf(cholirf,h);

%% 7. Bootstrap the confidence bands and plot the IRFs
cumulate = [1,2,3,4];
[bootMAX, upper, lower, boot_beta] = bootstrapMAX(y,p,c,beta,residuals,nboot,horizon,prc,cumulate);
% [bootMAX, upper, lower, boot_beta] = bootstrapMAX_level(y,p,c,beta,residuals,2000,horizon,95); % for level

%% Plot
% cumulate = []; % for level
varnames = {'TFP', 'SP500', 'GDP', 'CONS'};
plotirf_partial(max_irf,upper,lower,varnames,shockname,prc,cumulate)
% plotirf_partial(max_irf,upper,lower,varnames,shockname,prc) % for level

%% 8. Perform forecast error variance decomposition

% Compute the total variation
irf_chol_cumu = cumsum(cholirf,3);
tot_var = squeeze(cumsum(sum((irf_chol_cumu.^2),2),3));

% Cumulate the IRF from the Shock
irf_cumu = cumsum(max_irf,2);
fevd = cumsum((irf_cumu.^2),2)./tot_var;
plot_vardec(fevd,varnames, shockname)

% Total variation from Wold
wold_cumu = cumsum(wold,3);
tot_var_wold = zeros(N,horizon+1);
temp = 0;
for hh = 1:horizon+1
    temp = temp + wold_cumu(:,:,hh)*sigma*wold_cumu(:,:,hh)';
    tot_var_wold(:,hh) = diag(temp);
end



%% Givens rotations

syms theta_12 theta_13 theta_14 theta_23 theta_24 theta_34

H_12 = [cos(theta_12) sin(theta_12) 0 0;-sin(theta_12) cos(theta_12) 0 0;0 0 1 0;0 0 0 1];
H_13 = [cos(theta_13) 0 sin(theta_13) 0;0 1 0 0;-sin(theta_13) 0 cos(theta_13) 0;0 0 0 1];
H_14 = [cos(theta_14) 0 0 sin(theta_14);0 1 0 0;0 0 1 0;-sin(theta_14) 0 0 cos(theta_14)];
H_23 = [1 0 0 0;0 cos(theta_23) sin(theta_23) 0;0 -sin(theta_23) cos(theta_23) 0;0 0 0 1];
H_24 = [1 0 0 0;0 cos(theta_24) 0 sin(theta_24);0 0 1 0;0 -sin(theta_24) 0 cos(theta_24)];
H_34 = [1 0 0 0;0 1 0 0;0 0 cos(theta_34) sin(theta_34);0 0 -sin(theta_34) cos(theta_34)];

BigH = H_12*H_13*H_14*H_23*H_24*H_34;


