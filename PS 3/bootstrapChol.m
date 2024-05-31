function [bootchol, upper, lower, boot_beta] = bootstrapChol(y,p,c,beta,residuals,nboot,horizon,prc,cumulate,scaling)

% Function to compute bootstrapped Cholesky IRFs using the resampling method

% Inputs:   y           = TxN matrix of original data
%           p           = integer VAR lag order
%           c           = 1 if constant required
%           beta        = (Np+1 x N) matrix of estimated coefficients (Np x N) if
%                          no constant is included
%           residuals   = (T-p) x N matrix of OLS residuals from VAR(p)
%                          estimation
%           nboot       = integer number of bootstrap iterations
%           horizon     = integer horizon for the IRFs
%           prc         = integer between 0 and 100 to select size of bands
%           scaling = 2 x 1 vector where the first argument is the variable
%                     to be used for scaleing and the second is the shock size

% Outputs:  bootchol    = N x N x horizon + 1 x nboot array of
%                         bootstrapped Cholesky IRFs
%           upper       = N x N x horizon + 1 array of upper percentiles
%           lower       = N x N x horizon + 1 array of lower percentiles
%           boot_beta   = N x Np+1 x nboot array of bootstrapped VAR
%                         coefficient estimates

[T, N] = size(y);
bootchol = zeros(N,N,horizon+1,nboot);
boot_beta = zeros(N, size(beta,1), nboot);

for b=1:nboot
    
    % Bootstrap a new data set
    varboot = bootstrapVAR(y,p,c,beta,residuals);
    
    % Compute the new VAR coefficients and save
    [betaloop, err_loop] = VAR(varboot, p, 1);
    boot_beta(:,:,b) = betaloop';
    
    % Compute the wold IRF and save
    wold_loop = woldirf(betaloop,c,p,horizon);
    
    % Compute Cholesky factor
    sigma_loop = (err_loop' * err_loop) ./ (T - 1 - p - N*p);    
    S_loop = chol(sigma_loop, 'lower');
    
    % Compute the Cholsky IRFs
    if nargin > 9
        cholirf_loop = choleskyIRF(wold_loop, S_loop, scaling);
    else
        cholirf_loop = choleskyIRF(wold_loop, S_loop);
    end

    % Cumulate where necessary
    for i=1:N
        cholirf_loop(cumulate,i,:) = cumsum((cholirf_loop(cumulate,i,:)),3);
    end
    
    bootchol(:,:,:,b) = cholirf_loop;
    
end


up = (50 + prc/2);
low = (50 - prc/2);

% Extract the desired percentiles from the bootstrap distribution of the 
% Wold IRFs
upper = prctile(bootchol,up,4);
lower = prctile(bootchol,low,4);

end