 %%% Forecasting
%%% Macroeconometrics

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
clear; clc; close all;

%% Forecasting an AR(1):
% First, using the FULL SAMPLE:

PCEC = readtable('PCEC.csv');
GDPDEF = readtable('GDPDEF.csv');
%y = diff(log(PCEC.PCEC))*100;
y = diff(log(GDPDEF.GDPDEF))*400;

%dates = PCEC.DATE(2:end);
dates = GDPDEF.DATE(2:end);

figure;
plot(dates, y, "LineWidth", 2)
recessionplot
yline(0, "--r")
%title('Personal Consumption Expenditures Growth')
title('Inflation')

% First, estimate AR(1) by OLS:
phi_ols=my_AR(y,1,0); % 1 lag, no constant (think how to enter the constant for your PS!)

%%
% One period ahead forecast:

y_T=y(end);
forecast_1=phi_ols*y_T;

%disp(['The one period ahead forecast of the growth rate of Personal Consumption Expenditures is :' num2str(forecast_1)])
disp(['The one period ahead forecast of inflation is :' num2str(forecast_1)])

% j periods ahead forecasts, j=1,2,...,10:
H = 10;
forecast_j=zeros(H,1);

for j=1:H
    forecast_j(j)=phi_ols^(j)*y_T;
end

% Alternatively, using the previous period ahead forecast:

forecast_j=zeros(H,1);
forecast_j(1)=phi_ols*y_T;

for j=2:H
    forecast_j(j)=phi_ols*forecast_j(j-1);
end

% Plot of the series + forecasts:

%dates = 1947.5:0.25:2023.25;
datesH = dates(1):calmonths(3):dates(end)+calmonths(3*H);

figure;
plot(datesH,[y;NaN(H,1)],'LineWidth',2), axis tight; hold on;
plot(datesH,[NaN(size(y,1)-1,1);y_T;forecast_j],'LineWidth',2);
line(get(gca,'Xlim'),[0 0],'Color',[0 0 0],'LineStyle','--','LineWidth',1.5); hold off;
yline(mean(y), '--r')
xlabel('Date','FontWeight','Bold')
ylabel('%','FontWeight','Bold')
% legend('PCEC Growth','Forecast')
% title('PCEC Growth')
legend('Inflation','Forecast')
title('Inflation')
set(gca,'Fontsize',14)

iterfor_2 = forecast_j(2);
iterfor_3 = forecast_j(3);

%% Direct forecasts:

% Two periods ahead:
y_t=y(3:end);
y_t_2=y(1:end-2);

ols_2=OLS(y_t,y_t_2,0);
dirfor_2= ols_2*y_t(end);

% Three periods ahead:
y_t=y(4:end);
y_t_3=y(1:end-3);
X=[ones(length(y_t_3)) y_t_3];
ols_3=OLS(y_t,y_t_3,0);
dirfor_3=ols_3*y_t(end);

table(dirfor_2,dirfor_3, iterfor_2, iterfor_3)

%% Forecasting an AR(3):

% First, estimate AR(3) by OLS:
phi_ols=my_AR(y,3,0); % 1 lag, no constant

% Construct F_hat:
F_hat=[phi_ols';eye(2) zeros(2,1)];

% j periods ahead forecasts, j=1,2,...,H:
forecast_j_3=zeros(H,1);
y_T_1=y(end-1);
y_T_2=y(end-2);
Z_T=[y_T,y_T_1,y_T_2]';

for j=1:H
    
    Floop=F_hat^(j)*Z_T;
    forecast_j_3(j)=Floop(1,1);
    
end

figure;
plot(datesH,[y;NaN(H,1)],'LineWidth',2), axis tight; hold on;
plot(datesH,[NaN(size(y,1)-1,1);y_T;forecast_j_3],'LineWidth',2);
line(get(gca,'Xlim'),[0 0],'Color',[0 0 0],'LineStyle','--','LineWidth',1.5); hold off;
xlabel('Date','FontWeight','Bold')
ylabel('%','FontWeight','Bold')
% legend('PCEC Growth','Forecast')
% title('PCEC Growth')
legend('Inflation','Forecast')
title('Inflation')
set(gca,'Fontsize',14)

%% Simple pseudo out-of-sample one period ahead forecasts - AR(1):

% First of all, we have to estimate the parameters using half of the
% sample. In this case, as we have T=307, I will use T_0=154. Hence:

T = length(y);  % Total sample
T_0= ceil(T/2); % Training sample
T_1 = T - T_0;   % Testing sample

forecasts1=zeros(T_1,1);

% Notice that here I make a loop that repeats the steps until the end of
% the sample:

for k=1:T_1
    
    phi_ols=my_AR(y(1:T_0 + k -1),1,0);
    forecasts1(k,1)=phi_ols*y(T_0+k-1); % 1 period ahead   

end

% Plot:
y_hat1=[NaN(T_0-1,1);y(T_0); forecasts1];

figure;
plot(dates,y,'LineWidth',2), axis tight; hold on;
plot(dates,y_hat1,'LineWidth',2);
line(get(gca,'Xlim'),[0 0],'Color',[0 0 0],'LineStyle','--','LineWidth',1.5); hold off;
xlabel('Date','FontWeight','Bold')
ylabel('%','FontWeight','Bold')
legend('PCEC Growth','Forecast')
title('PCEC Growth')
set(gca,'Fontsize',14)

%% Computing 1 period ahead mean squared errors
w1 = y(T_0+1:end)-forecasts1(:,1);
MSE1=mean(w1.^2);

disp('The MSE associated to the 1 ahead pseudo out of sample forecast using half of the sample are given by :')
table(MSE1)




