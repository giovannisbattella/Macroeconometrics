function fevd = variance_decomp(irf, shock)

% Function to compute forecast error variance decomposition. If only one
% shock is of interest, then input its position.

% Inputs:   irf = (N x N x horizon) array of impulse responses 
%           shock = if a specific shock is required   

% Output:  fevd = N x N x horizon array of forecast error shares

[N,~,horizon] = size(irf);

% Square the IRFs
irf = irf.^2;

% Compute the total variation for each variable by summing over shocks
% (dimension 2) first and then summing over horizons (dimension 3).
% We use the cumsum() command to get the total variation at each horizon
total_variation = squeeze(cumsum(sum(irf,2),3));

% Now compute the variation of each individual shock over the horizon and
% divide it by the total variation

if nargin > 1 % For a specific shock
    shock_variation = cumsum(squeeze(irf(:,shock,:)),2);
    fevd = shock_variation ./ total_variation;
else % For all shocks
    fevd = zeros(N,N,horizon);
    for s=1:N
        shock_variation = cumsum(squeeze(irf(:,s,:)),2);
        fevd(:,s,:) = shock_variation ./ total_variation;
    end

end


