function [bootbq, upper, lower, boot_beta] = bootstrapBQ_corrected(y,p,c,beta,residuals,nboot1,horizon,prc,cumulate,nboot2,scaling)

% Function to compute bootstrapped Cholesky IRFs using the resampling method
% with Kilian correction

% Inputs:   y           = TxN matrix of original data
%           p           = integer VAR lag order
%           c           = 1 if constant required
%           beta        = (Np+1 x N) matrix of estimated coefficients (Np x N) if
%                          no constant is included
%           residuals   = (T-p) x N matrix of OLS residuals from VAR(p)
%                          estimation
%           nboot1      = integer number of bootstrap iterations for bias
%                         computations
%           nboot2      = integer number of bootstrap iterations for bias
%                         corrected bootstrap
%           horizon     = integer horizon for the IRFs
%           prc         = integer between 0 and 100 to select size of bands
%           scaling = 2 x 1 vector where the first argument is the variable
%                     to be used for scaleing and the second is the shock size

% Outputs:  bootbq    = N x N x horizon + 1 x nboot array of
%                         bootstrapped BQ LR IRFs
%           upper       = N x N x horizon + 1 array of upper percentiles
%           lower       = N x N x horizon + 1 array of lower percentiles
%           boot_beta   = N x Np+1 x nboot array of bootstrapped VAR
%                         coefficient estimates

[T, N] = size(y);
boot_beta = zeros(N, size(beta,1), nboot1);

for b=1:nboot1
    
    % Bootstrap a new data set
    varboot = bootstrapVAR(y,p,c,beta,residuals);
    
    % Compute the new VAR coefficients and save
    [betaloop, ~] = VAR(varboot, p, 1);
    boot_beta(:,:,b) = betaloop';
    
end

[Beta, ~] = remove_bias(beta,c,p,boot_beta);

if nargin > 10
    [bootbq, upper, lower, boot_beta] = bootstrapBQ(y,p,c,Beta',residuals,nboot2,horizon,prc,cumulate,scaling);
else
    [bootbq, upper, lower, boot_beta] = bootstrapBQ(y,p,c,Beta',residuals,nboot2,horizon,prc,cumulate);
end

end