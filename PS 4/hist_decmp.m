function [histdec, ystar] = hist_decmp(y, beta, residuals, c, p, K, series)

% Function to compute historical decomposition of series of choice 

% Inputs:   y           = TxN matrix of original data
%           p           = integer VAR lag order
%           c           = 1 if constant required
%           beta        = (Np+1 x N) matrix of estimated coefficients (Np x N) if
%                          no constant is included
%           residuals   = (T-p) x N matrix of OLS residuals from VAR(p)
%                          estimation

% Outputs:  histdec     = T-p x N matrix of contributions of shocks to series
%           ystar       = demeaned series

[T, N] = size(y);

% Strucutral shocks
struct_shock = (K\residuals')';

% Demean the data
ystar = y(p+1:end,series) - mean(y(p+1:end,series),1);

% Compute the MA coefficients of the entire sample
wold_long = woldirf(beta, c, p, T-p);
ma_coeff = zeros(N,N,T-p);

for i=1:T-p
    ma_coeff(:,:,i) = wold_long(:,:,i) * K;
end

% Compute the historical decomposition for the series of choice
histdec = zeros(T-p, N);

for t=1:T-p % for each time period
    for j=1:N % for each shock                     
        histdec(t,j) =  squeeze(ma_coeff(series,j,1:t))' * flipud(struct_shock(1:t,j));           
    end
end

end