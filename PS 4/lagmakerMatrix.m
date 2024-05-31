function x = lagmakerMatrix(y,p)

% Function to create the matrix of regressors according to the 
% SUR representaiton for a VAR

% Inputs:   y = T x N matrix of endogeneous variables
%           p = VAR lag order
% Outputs:  x = T-p x Np matrix of lagged dependent variables as regressors

[T, N] = size(y);

x = zeros(T-p, N*p);
counter = 0;
for i=1:p
    for j=1:N
        counter = counter + 1;
        x(:,counter) = y(p+1-i:T-i,j);      
    end
end



end