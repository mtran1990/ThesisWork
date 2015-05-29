% Gauss-Newton Test

addpath(genpath('./util'));

%% test 1
% loc = [0 0; 3 0; 2 2];
% R   = [2 1 2]';
% 
% [beta,iter] = gaussNewton(loc,R);
% disp(beta);
% disp(iter);

%% test 2
loc = [0 0; 2 0];
R   = [1.5 1.5]';

[beta,iter] = gaussNewton(loc,R,[1;0]);
disp(beta);
disp(iter);