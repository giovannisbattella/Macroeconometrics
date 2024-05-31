function [aic, bic, hq] = aicbic(y , pmax, c)

% Function to determine the lag order of a VAR(p) according to the Akaike
% and Bayesian information criterion and Hannan-Quinn criterion

% Inputs:   y = T x N matrix of endogeneous variables
%           pmax = integer maximum VAR lag order
%           c = 1 if constant required
% Outputs:  aic = integer lag length recommended by AIC
%           bic = integer lag length recommended by BIC
%           hq  = integer lag length recommended by HQ

    [T,N] = size(y);
    for p = 1:pmax
        [beta, e] = VAR(y, p, c);
        a = log(det((e'*e)/T));
        b = (p*N^2 + N)/T;
        aic(p) = a + 2*b;
        bic(p) = a + log(T)*b;
        hq(p) = a + 2*b*log(log(T));
    end
    [minaic,aic] = min(aic);
    [minbic,bic] = min(bic);
    [minhq,hq] = min(hq);

end