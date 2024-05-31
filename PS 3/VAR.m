function [beta, residuals]  = VAR(y,p,c)

% Function to estimate a VAR(p) with or without constant using OLS
% Inputs:   y = T x N matrix of endogeneous variables
%           p = VAR lag order
%           c = 1 if constant required
% Outputs:  beta = (Np+1 x N) matrix of estimated coefficients (Np x N) if
%                   no constant is included
%           residuals = (T-p x N) matrix of OLS residuals

[T, ~] = size(y);
yfinal = y(p+1:T,:);
if c == 1
    X = [ones(T-p,1), lagmakerMatrix(y,p)];
else
    X = lagmakerMatrix(y,p);
end

beta = (X'*X)\X'*yfinal;
residuals = yfinal - X*beta;

end