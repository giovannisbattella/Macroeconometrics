function [sys] = mysystem(k,Omega,C1)

    K = [k(1), k(2); k(3), k(4)];
    
    n = size(Omega,2);
    first = inv(K)*Omega*inv(K)' - eye(n); % 2 by 2 matrix
    f1 = first(1,1);
    f2 = first(1,2);
    f3 = first(2,2);
    second = C1(1,:)*K(:,2); % scalar
    f4 = second;
    sys = [f1;f2;f3;f4];

end
