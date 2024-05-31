function [ynext, draws] = bootstrapVAR(y,p,c,beta,residuals)

% Function to bootstrap a VAR using the resampling technique

% Inputs:   y           = TxN matrix of original data
%           p           = integer VAR lag order
%           c           = 1 if constant required
%           beta        = (Np+1 x N) matrix of estimated coefficients (Np x N) if
%                          no constant is included
%           residuals   = (T-p) x N matrix of OLS residuals from VAR(p)
%                          estimation

% Outputs:  ynext       = T x N matrix of bootstrapped time series

[T, N] = size(y);
yinit = y(1:p, :);

if c==1
    const = beta(1,:);
    pi = beta(2:end,:);
else
    const = 0;  
    pi = beta;
end

ynext=zeros(T,N);
ynext(1:p,:) = yinit;
yinit=reshape(flipud(yinit)',1,[]); % Initial data to feed into model.
draws = zeros(T-p,1);
for i=1:T-p 
  draws(i) = randi(T-p);
  ynext(p+i,:)= const + yinit*pi + residuals(draws(i),:);  % Draw a random integer to jointly select from residuals
  yinit = reshape(flipud(ynext(i+1:p+i,:))',1,[]);  % New values to feed into model
end


end