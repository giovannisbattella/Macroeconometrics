function [output] = vecm_rr(y,p,c,r,h)

    % First-difference the data in y to get (1-L)yt
    yy = diff(y);

    % Create X_t-1, the p-1 lags of (1-L)yt
    X = lagmakerMatrix(yy, p-1);

    % Regress (1-L)yt onto X_t-1 and get the residuals R0
    yt = yy(p:end,:);
    R0 = yt - X*((X'*X)\X'*yt);

    % Regress Y_t-1 onto X_t-1 and get the residuals R1
    R1 = y(1:end-p,:) - X*((X'*X)\X'*y(1:end-p,:));

    % Compute the matrices S_ij
    [T,N] = size(R0);
    S_00 = R0'*R0/T;
    S_11 = R1'*R1/T;
    S_01 = R0'*R1/T;
    S_10 = R1'*R0/T;

    % Compute the r orthonormal eigenvectors
    M = S_11^(-0.5)*S_10*inv(S_00)*S_01*S_11^(-0.5);
    [v,D] = eigs(M,[], r);

    % Recover beta(tilde) and alpha(tilde) to get Pi
    beta = S_11^(-1/2)*v;
    alpha = S_01*beta*inv(beta'*S_11*beta); % can get this from OLS too
    pi = alpha*beta';

    % Get the gamma parameters from OLS regressions 
    lhs = yy(p:end,:);
    if c == 0
        rhs = [y(p:end-1,:)*beta, X];    
        psi = (rhs'*rhs)\rhs'*lhs;
        gamma = psi(r+1:end,:)';
        % alpha1 = psi(1:r,:);   % same as above
        eps = lhs - rhs*psi;
    else
        rhs = [ones(T,1), y(p:end-1,:)*beta, X];    
        psi = (rhs'*rhs)\rhs'*lhs;
        gamma = psi(r+2:end,:)';
        % alpha1 = psi(2:r+1,:);   % same as above
        eps = lhs - rhs*psi;
    end
    % pi1 = alpha*beta';    % same as above

    % Recover the VAR parameters
    A1 = pi + eye(N) + gamma(:,1:N);
    Ap = -gamma(:, end-N+1:end);

    if p > 2

        Ai = zeros(N,N,p-2);
        for i = 2:p-1
            Ai(:,:,i-1) = gamma(:,1+N*(i-1):i*N) -  gamma(:, 1+N*(i-2):(i-1)*N);
        end

    end
    par_var = [A1, reshape(Ai,N,(p-2)*N), Ap]';
    wold = woldirf(par_var,0,p,h);

    output.wold = wold;
    output.par_var = par_var;
    output.beta = beta;
    output.alpha = alpha;
    output.eps = eps;

end
