function [cholirf] = choleskyIRF(wold, S, scaling)

% Function to compute the point estimate of the IRF of a VAR identified
% using Cholesky

% Inputs:   wold    = (N x N x horizon +1) array of Wold IRFs
%           S       = N x N lower triangular matrix Cholesky factor
%           scaling = 2 x 1 vector where the first argument is the variable
%                     to be used for scaleing and the second is the shock size

% Outputs:  chol = N x N x horizon + 1 array of Cholesky identified IRFs

[N,~, horizon] = size(wold);
cholirf = zeros(N,N,horizon);

for h=1:horizon
    
    cholirf(:,:,h) = wold(:,:,h) * S;
    
end

if nargin > 2
    cholirf = cholirf ./ (cholirf(scaling(1),scaling(1),1) * 1/scaling(2));
end

end