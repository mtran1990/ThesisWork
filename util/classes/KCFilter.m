% Kalman Consensus Filter
classdef KCFilter < handle
    %% properties
    properties (SetAccess = private)
        F
        H
        Q
        R
        sigA2
        xu
        P
        xp
        u
        U
    end
    
    %% methods
    methods (Access = public)
        
        function obj = KCFilter(params)
            dt = params.dt;
            F_ = [1 dt; 0 1];
            obj.F = [F_ zeros(2,2); zeros(2,2) F_];

            obj.H = [1 0 0 0; 0 0 1 0];
            
            
            obj.sigA2 = params.sigA2;
            sigM2 = params.M2;
            
            Q_ = [dt^4/4 dt^3/2; dt^3/2 dt^2]*obj.sigA2;
            obj.Q = [Q_ zeros(2); zeros(2) Q_];
            
            obj.R = [sigM2 0; 0 sigM2];
            
            obj.xp = params.x0;
            obj.P = params.P0;
        end
        
        function addMeasurement(obj, z, Pd, Pg, lam)
            
            % can't calculate weights without measurements
            if(nargin == 4)
                % get a matrix of betas to elementwise multiply with
                beta = calculateWeights(z,Pd,Pg,lam);
                beta = ones(2,1)*beta;

                % don't need the first element it seems (corresponds to false
                % alarm)
                % need to sum across rows
                zVec = sum(beta(:,2:end).*z,2);

                obj.u = obj.H.'*obj.R\zVec;
                obj.U = obj.H.'*obj.R\obj.H;
            else
                obj.u = 0;
                obj.U = 0;
            end
            
        end
        
        function iterFilter(obj, xj, u, U)
            
            % probably need to rewrite this part
            % thinking of returning a boolean, where true is returned when
            % the track was updated with measurements and false is returned
            % otherwise
            % might be able to use that to help update the track states
            
            eps = 0.25;
            y = sum(u,2)+obj.u;
            S = sum(U,3)+obj.U;
            
            % update step
            M = (inv(obj.P)+S);
            obj.xu = obj.xp + M\(y-S*obj.xp)+eps*M\sum(xj-obj.xp);
            
            % prediction step
            obj.P = obj.F*M*obj.F.'+obj.Q;
            obj.xp = obj.F*obj.xu;            
            
        end
        
    end
    
    %% private methods
    methods (Access = private)
        
        function Beta = calculateWeights(obj, z, Pd, Pg, lam)
            % z  : 2xN vector of N measurements
            % Pd : Probability of Detection
            % Pg : Probability of measurement lying inside the gate
            % lam: Mean of the poisson distribution of false alarms            
            
            zTild = z - obj.H*obj.xp;
            
            % innovation covariance (covariance of zTild?)
            innCov = obj.H*obj.P*obj.H'+obj.R;
            
            % dimension of target
            dim = 2;
            
            b = (2*pi)^(dim/2)*lam*sqrt(det(innCov))*(1-Pd*Pg)/Pd;
            
            % N+1 length vector representing the probability that a
            % measurement is associated with this track
            % N: # of measurements
            % first element is the probability of being associated with no
            % tracks
            N = size(z,2);
            Beta = zeros(1,N+1);
            
            % might be a way to do this using vectors / matrices instead
            den = b;
            for k = 1:N
                den = den+exp(-zTild(:,k).'*(innCov\zTild(:,k))/2);
            end
            
            Beta(1) = b/den;
            
            for k = 2:N+1
                
                num = exp(-zTild(:,k).'*(innCov\zTild(:,k))/2);
                Beta(k) = num/den;
                
            end
            
            
        end
        
        
    end
end
        
        
    