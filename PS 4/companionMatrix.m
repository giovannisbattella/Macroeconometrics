function [comp, N] = companionMatrix(beta, c, p)

% Function to create the F matrix of the companion form associated with a
% VAR(P)

% Inputs:   beta = (Np+1 x N) matrix of estimated coefficients (Np x N) if
%                   no constant is included
%           c = 1 if constant required
%           p = VAR lag order

% Outputs:  comp = Np x Np matrix F of the companion form
%           N = number of dependent variables in the VAR

if c==1
    N = (size(beta,1) - 1) / p;
    beta = beta(2:end,:);
else
    N = size(beta,1) / p;
end

comp = zeros(N*p,N*p);

comp(1:N,:) = beta';
comp(N+1:end,1:N*(p-1)) = eye(N*(p-1));


end