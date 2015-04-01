% Gauss-Newton algorithm
% assuming "r_i" is the distance error

function [beta,iter] = gaussNewton(loc, R, beta, threshAbs, maxit)
% inputs:
% loc   : locations of the sensors that produced the R vector, the size
%         should be nx2, where n is the number of measurements
% R     : range values of the measurements, the size should be nx1
% beta  : optional, a starting guess for the algorithm

if(nargin<3), beta = [0; 0]; end
if(nargin<4), threshAbs = 1e-4; end
if(nargin<5), maxit = 100; end

iter = 0;

while(1)
    
    bOld = beta;
    iter = iter+1;
    
    J = calcJacobian(loc,beta);
    r = calcR(loc,beta,R);
    
    beta = beta-(J.'*J)\J.'*r;
    
    eAbs = sum(abs(beta-bOld));
    
    if(eAbs<threshAbs || iter>maxit)
        break;
    end
            
end

if(sum(isnan(beta)))
    error('Possible error in results');
end

end

function J = calcJacobian(loc,beta)

[m,~] = size(loc);

xy = ones(m,1)*beta.';

J = 2*(xy-loc);

end

function r = calcR(loc,beta,R)

[m,~] = size(loc);

xy = ones(m,1)*beta.';

r = sum((xy-loc).^2,2)-R.^2;

end