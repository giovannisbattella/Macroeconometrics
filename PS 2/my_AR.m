function [beta, sigma2, tstat, res, se] = my_AR(y,p,c)

    yf = y(p+1:end);
    X = lagmakerMatrix(y,p);

    if c == 1
        X = [ones(length(yf),1), X];
    end

    K = size(X,2);
    T = size(yf,1);

    % OLS
    beta = (X'*X)\X'*yf;

    % Residuals
    res = yf - X*beta;

    % Sigma2
    sigma2 = (res'*res)/(T-K);

    % Standard errors
    se = sqrt(diag(sigma2.*inv(X'*X)));

    % T-statistics
    tstat = abs(beta./se);

end