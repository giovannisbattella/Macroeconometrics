function [sys] = mysystem1(h,D1)

H = [h(1), h(2); h(3), h(4)];
n = size(H,2);
first = H'*H - eye(n);
f1 = first(1,1);
f2 = first(1,2);
f3 = first(2,2);
second = D1(1,:)*H(:,2); % scalar
f4 = second;
sys = [f1;f2;f3;f4];
