%%% Quantitative and Statistical Methods II
%%% Problem Set 3

%% Importing the data:
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
clc; clear;

%% Read in the data
originaldata = readtable('newdata_PS3.csv');

% We use only data between Q2 1954 and Q4 2022
data = originaldata(31:304,2:6);
dates = originaldata.DATE(31:304);

% Convert into matrix format
data = table2array(data);
UNRATE = data(:,3);
GS10 = data(:,4);
GDPDEF = data(:,1);
GDPC1 = data(:,2);
FEDFUNDS = data(:,5);
%% 1) 

Spread=GS10-FEDFUNDS;
GDPC1_gr=diff(log(GDPC1))*100;
GDPDEF_gr=diff(log(GDPDEF))*100;

% Plot the data
dates_data = dates;
figure;
subplot(2,2,1)
axis tight;
hold on;
plot(dates_data,Spread)
plot(dates_data,GS10)
plot(dates_data,FEDFUNDS)
yline(0, '--r')
hold off;
recessionplot
title('Interest Rates')
legend({'Spread', 'GS10', 'FFR'})

subplot(2,2,2)
plot(dates_data(2:end),GDPC1_gr)
axis tight;
recessionplot
title('Real GDP growth')

subplot(2,2,3)
plot(dates_data(2:end),GDPDEF_gr)
axis tight;
recessionplot
title('GDP Deflator Growth')

subplot(2,2,4)
plot(dates_data,UNRATE)
axis tight;
recessionplot
title('Unemployment Rate')


%% 2) 

% Estimate the reduced-form VAR(4) with finaldata:

finaldata=[GDPC1_gr,UNRATE(2:end),GDPDEF_gr, FEDFUNDS(2:end),Spread(2:end)]; 

% We start from the second observation for the rates since we lost one
% observation by taking diff(log(.)) of GDPC1 and GDPDEF.

% First of all, we have to create the (T-p) x n and (T-p) x np+1 matrices 
% of the SUR representation, i.e. Y and X respectively:

[T,N]=size(finaldata); % # of variables
p=4; % 4 lags
c=1; % including constant

[pi_hat,err]=VAR(finaldata,p,c); % VAR estimation

%% 3)

% Wold representation impulse responses:

% Populate the companion form matrix
BigA=[pi_hat(2:end,:)'; eye(N*p-N) zeros(N*p-N,N)]; % BigA companion form, np x np matrix
BigA1=companionMatrix(pi_hat, c, p);

% Check stability
ev = abs(eig(BigA));
evmax = ['The maximum eigenvalue is ', num2str(max(ev)), '.'];
disp(evmax)

hor=40;

C=zeros(N,N,hor+1);

for j=1:hor+1
    BigC=BigA^(j-1);
    C(:,:,j)=BigC(1:N,1:N); % Impulse response functions of the Wold representation
end

% Alternatively, use the function woldiirf
CC = woldirf(pi_hat, c, p, hor);


%% Bootstrap the Wold CI
nboot1 = 1000;
nboot2 = 2000;
prc = 68;
cumulate = [];
[bootwold, upper, lower, boot_beta] = ...
    bootstrapWold_corrected(finaldata,p,c,pi_hat,err,nboot1,hor,prc,cumulate,nboot2);

%%
% Plot the IRFs:
colorBNDS=[0.8 0.8 0.8];
VARnames={'GDP';'UNR';'INF';'FFR';'SPR'};
Shocknames={'$u_1 \uparrow$';'$u_2 \uparrow$';'$u_3 \uparrow$';'$u_4 \uparrow$';'$u_5 \uparrow$'};

figure;
counter = 0;
for k=1:N % variable k
    for j=1:N % shock j
    counter = counter + 1;
    subplot(N,N,counter)
    fill([0:hor fliplr(0:hor)]' ,[squeeze(upper(k,j,:)); flipud(squeeze(lower(k,j,:)))],...
        colorBNDS,'EdgeColor','None'); hold on;
    plot(0:hor,squeeze(C(k,j,:)),'LineWidth',3.5,'Color','k'); hold on;
    set(gca,'FontSize',16)
    line(get(gca,'Xlim'),[0 0],'Color',[1 0 0],'LineStyle','-','LineWidth',1.5); hold off;
    
    if j ==1
        ylabel(strcat(VARnames{k}), 'FontSize', 18);
    end
    
    if k == 1
        title(strcat(Shocknames{j}), 'FontSize', 22,'Interpreter','Latex');
    end
    xlim([0 hor]);
    end
end

legend({'68% confidence bands','IRF'},'FontSize',16,'Orientation','Horizontal')

%% 4) and 5)

% Cholesky:

omega=(err'*err)./(T - N*p-1-p); % Estimate of omega
S=chol(omega,'lower'); % Cholesky factorization, lower triangular matrix

D=zeros(N,N,hor+1);

for i=1:hor+1
    D(:,:,i)=C(:,:,i)*S; % Cholesky wold respesentation
end

% Alternatively use the function choleskyIRF
DD = choleskyIRF(C, S);

DDcumu = DD;
cumu = [1, 3];
DDcumu(cumu,:,:) = cumsum(DDcumu(cumu,:,:),3);


%% Bootstrap the confidence bands
[bootchol, upper, lower, boot_beta] = ...
    bootstrapChol_corrected(finaldata,p,c,pi_hat,err,nboot1,hor,prc,cumu,nboot2);

%% Plot the IRFs:
Dplot = DDcumu;

Shocknames={'';'';'';'Monetary Policy';''};
VARnames={'GDP';'UNR';'PRICES';'FFR';'SPR'};

figure;
counter = 0;
for k=1:N % variable k
    for j=1:N % shock j
    counter = counter + 1;
    subplot(N,N,counter)
    fill([0:hor fliplr(0:hor)]' ,[squeeze(upper(k,j,:)); flipud(squeeze(lower(k,j,:)))],...
        colorBNDS,'EdgeColor','None'); hold on;
    plot(0:hor,squeeze(Dplot(k,j,:)),'LineWidth',3.5,'Color','k'); hold on;
    set(gca,'FontSize',16)
    line(get(gca,'Xlim'),[0 0],'Color',[1 0 0],'LineStyle','-','LineWidth',1.5); hold off;
    
    if j ==1
        ylabel(strcat(VARnames{k}), 'FontSize', 18);
    end
    
    if k == 1
        title(strcat(Shocknames{j}), 'FontSize', 18);
    end
    xlim([0 hor-1]);
    end
end

legend({'68% confidence bands','IRF'},'FontSize',16,'Orientation','Horizontal')


