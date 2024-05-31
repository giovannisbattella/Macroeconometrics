function [bqirf] = bqIRF(wold, K, scaling)

% Function to compute the point estimate of the IRF of a VAR identified
% using BQ long-run restrictions

% Inputs:   wold    = (N x N x horizon +1) array of Wold IRFs
%           K       = N x N lower triangular matrix Cholesky factor of LR
%           scaling = 2 x 1 vector where the first argument is the variable
%                     to be used for scaleing and the second is the shock size

% Outputs:  bqirf = N x N x horizon + 1 array of LR identified IRFs

[N,~, horizon] = size(wold);
bqirf = zeros(N,N,horizon);

for h=1:horizon
    
    bqirf(:,:,h) = wold(:,:,h) * K;
    
end

if nargin > 2
    bqirf = bqirf ./ (bqirf(scaling(1),scaling(1),1) * 1/scaling(2));
end

end