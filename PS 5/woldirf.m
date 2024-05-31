function irfwold = woldirf(beta, c, p, horizon)

% Function to compute the matrices of the Wold representation of a
% stationary VAR(p)

% Inputs:   beta = (Np+1 x N) matrix of estimated coefficients (Np x N) if
%                   no constant is included
%           c = 1 if constant required
%           p = VAR lag order
%           horizon = how many of the Wold matrices (+1) will be computed

% Outputs:  irfwold = N x N x horizon + 1 array of Wold coefficient
%                     matrices

[BigA, N] = companionMatrix(beta,c,p);

irfwold = zeros(N,N,horizon + 1);

for h=1:horizon+1
    
    temp = BigA^(h-1);
    irfwold(:,:,h) = temp(1:N, 1:N);
    
end

end
