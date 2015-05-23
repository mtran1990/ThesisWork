% lapjv test
% Linear Assignment Problem JV

% add path
addpath('../util/lapjv');

% infinity example
x = magic(5);
x(2) = Inf;
disp(x);
[rowsol,cost] = lapjv(magic(5));
disp(rowsol); % 3 2 1 5 4
disp(cost);   %15

% % non-square assigment
% n=100;
% A=1./randn(n);
% tic
% [a,b]=lapjv(A);
% toc % about 0.2 sec
% A1=[A zeros(n,1)+max(max(A))];
% tic
% [a1,b1]=lapjv(A1);
% toc % about 0.01 sec. The nonsquare one can be done faster!
% %check results
% disp(norm(a-a1))
% disp(b-b)

% rm path
rmpath('../util/lapjv');