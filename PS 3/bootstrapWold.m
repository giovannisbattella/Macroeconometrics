function [bootwold, upper, lower, boot_beta] = bootstrapWold(y,p,c,beta,residuals,nboot,horizon,prc, cumulate)

% Function to compute bootstrapped Wold IRFs using the resampling method

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

% Outputs:  bootwold    = N x N x horizon + 1 x nboot array of
%                         bootstrapped Wold IRFs
%           upper       = N x N x horizon + 1 array of upper percentiles
%           lower       = N x N x horizon + 1 array of lower percentiles
%           boot_beta   = N x Np+1 x nboot array of bootstrapped VAR
%                         coefficient estimates

[~, N] = size(y);
bootwold = zeros(N,N,horizon+1,nboot);
boot_beta = zeros(N, size(beta,1), nboot);

for b=1:nboot
    
    % Bootstrap a new data set
    varboot = bootstrapVAR(y,p,c,beta,residuals);
    
    % Compute the new VAR coefficients and save
    [betaloop, ~] = VAR(varboot, p, 1);
    boot_beta(:,:,b) = betaloop';
    
    % Compute the wold IRF and save
    wold_loop = woldirf(betaloop,c,p,horizon);
    
    % Cumulate where necessary
    for i=1:N
        wold_loop(cumulate,i,:) = cumsum(squeeze(wold_loop(cumulate,i,:)),2);
    end
    
    bootwold(:,:,:,b) = wold_loop;
    
end


up = (50 + prc/2);
low = (50 - prc/2);

% Extract the desired percentiles from the bootstrap distribution of the 
% Wold IRFs
upper = prctile(bootwold,up,4);
lower = prctile(bootwold,low,4);

end