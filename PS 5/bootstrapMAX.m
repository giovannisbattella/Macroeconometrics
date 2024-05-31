function [bootMAX, upper, lower, boot_beta] = bootstrapMAX(y,p,c,beta,residuals,nboot,horizon,prc,cumulate)

[T, N] = size(y);
bootMAX = zeros(N,horizon+1,nboot);
boot_beta = zeros(N, size(beta,1), nboot);
startvalue=ones(3,1)./norm(ones(3,1));
options=optimset('Display', 'off', 'TolFun',.0000000001,'MaxIter',100000,'MaxFunEvals',100000); % more iterations

for b=1:nboot
    
    % Bootstrap a new data set
    varboot = bootstrapVAR(y,p,c,beta,residuals);
    
    % Compute the new VAR coefficients and save
    [betaloop, err_loop] = VAR(varboot, p, 1);
    boot_beta(:,:,b) = betaloop';
    
    % Compute the wold IRF
    wold_loop = woldirf(betaloop,c,p,horizon);
    
    % Compute Cholesky
    sigma_loop = (err_loop' * err_loop) ./ (T - 1 - p - N*p);
    S_loop = chol(sigma_loop, 'lower');
    chol_loop = choleskyIRF(wold_loop,S_loop);
    
    % Compute LR impact vector h 
    var = 1; %TFP is ordered first
    h_est_loop = fminsearch(@(h) lr_target(h,chol_loop,var),startvalue,options);
    h2_loop = [0;h_est_loop./norm(h_est_loop)];
    
    % Compute the IRFs
    max_effect_irfs_loop = partial_irf(chol_loop, h2_loop);

    % Cumulate where necessary
    max_effect_irfs_loop(cumulate,:) = cumsum((max_effect_irfs_loop(cumulate,:)),2);
    
    % Store
    bootMAX(:,:,b) = max_effect_irfs_loop;
    
end


up = (50 + prc/2);
low = (50 - prc/2);

% Extract the desired percentiles from the bootstrap distribution of the 
% IRFs
upper = prctile(bootMAX,up,3);
lower = prctile(bootMAX,low,3);

end