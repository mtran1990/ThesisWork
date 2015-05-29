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
        xdiff
        
        beta0
        pTild
    end
    
    %% methods
    methods (Access = public)
        
        function obj = KCFilter(params)
            dt = params.dt;
            F_ = [1 dt; 0 1];
            obj.F = [F_ zeros(2,2); zeros(2,2) F_];

            obj.H = [1 0 0 0; 0 0 1 0];
            
            
            obj.sigA2 = params.sigA2;
            sigM2 = params.sigM2;
            
            Q_ = [dt^4/4 dt^3/2; dt^3/2 dt^2]*obj.sigA2;
            obj.Q = [Q_ zeros(2); zeros(2) Q_];
            
            obj.R = [sigM2 0; 0 sigM2];
            
            obj.xp = params.x0;
            obj.P = params.P0;
        end
        
        function addMeasurement(obj, z, Pd, Pg, lam)
            
            % can't calculate weights without measurements
            if(nargin == 5)
                
                beta = obj.calculateWeights(z,Pd,Pg,lam);
                obj.beta0 = beta(1);
                obj.pTild = obj.calcPTild(beta,z);
                
                % get a matrix of betas to elementwise multiply with
                beta_tmp = ones(2,1)*beta;

                % don't need the first element it seems (corresponds to
                % false alarm)
                % need to sum across rows
                zVec = sum(beta_tmp(:,2:end).*z,2);

                obj.u = obj.H.'*(obj.R\zVec);
                obj.U = obj.H.'*(obj.R\obj.H);
                
                if(obj.u(1) < 0)
                    error('here');
                end
                
            else
                obj.u = zeros(4,1);
                obj.U = zeros(4,4);
            end
            
            % right now, I'm still calling this function even when a track
            % has no new measurements in order to reset the values of u, U,
            % and xdiff
            obj.xdiff = zeros(4,1);
        end
        
        function addTrack(obj,T2)
            if(~isempty(obj.u))
                % adds the info from another track
                obj.u = obj.u+T2.u;
                obj.U = obj.U+T2.U;

                obj.xdiff = obj.xdiff + (obj.xp-T2.xp);
            end
        end
        
        function updated = iterFilter(obj)
            
            % updated: true if filter updated with measurements
            %          false if filter updated with no measurements
            %          empty if filter did not update
            
            % if u is empty, track hasn't updated yet, so don't do anything
            % yet
            if(~isempty(obj.u))
            
                % if u is a column vectors of zeros, then the track didn't
                % receive any measurements
                if(isequal(zeros(4,1),obj.u))
                    % update step
                    xu_ = obj.xp;

                    % prediction step
                    % NOTE: not sure if this part is right, need to double
                    % check
                    obj.P = obj.F*obj.P*obj.F.'+obj.Q;
                    obj.xp = obj.F*xu_;
                    updated = false;
                else            
                    eps = 0.25;
                    y = obj.u;
                    S = obj.U;

                    S0 = (1-obj.beta0)*S;

                    % update step
                    Mtild = (inv(obj.P)+S);
                    K = Mtild\(obj.H).'/obj.R;
                    M = obj.beta0*obj.P+(1-obj.beta0)*inv(Mtild)+K*obj.pTild*K.';

                    xu_ = obj.xp + Mtild\(y-S0*obj.xp)+...
                        eps*M/(1+norm(M,'fro'))*obj.xdiff;                

                    % prediction step
                    P = obj.F*M*obj.F.'+obj.Q;
                    obj.xp = obj.F*xu_;
                    updated = true;
                    
                    if(P(1) < 0)
                        error('here');
                    end
                    
                    obj.P = P;
                    
                end

                obj.updateState(xu_);
            
            else
                
                updated = [];
            
            end
            
        end
        
        function cost = computeCost(obj,T2)
            % computes the cost between the current track and the given
            % track (using Mahalanobis distance)
            S = (obj.P+T2.P);
            xdiff = obj.xp-T2.xp;
            c = xdiff.'*S*xdiff;
            
            cost = c+log(det(S));
        end
        
        function S = getInnCovariance(obj)
            
            S = obj.H*obj.P*obj.H.'+obj.R;
            
        end
        
    end
    
    %% private methods
    methods (Access = private)
        
        function Beta = calculateWeights(obj, z, Pd, Pg, lam)
            % z  : 2xN vector of N measurements
            % Pd : Probability of Detection
            % Pg : Probability of measurement lying inside the gate
            % lam: Mean of the poisson distribution of false alarms            
            
            % need to convert tmp to be the same size as zTild
            tmp = obj.H*obj.xp;
            tmp = tmp*ones(1,size(z,2));
            
            zTild = z - tmp;
            
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
            
            for k = 1:N
                
                num = exp(-zTild(:,k).'*(innCov\zTild(:,k))/2);
                % first element in beta is for false alarms
                Beta(k+1) = num/den;
                
            end
            
            
        end
        
        function pTild = calcPTild(obj,beta,z)
            
            [M,N] = size(z);
            pTild = zeros(M);
            
            beta_tmp = ones(2,1)*beta(2:end);
            
            zi = sum(beta_tmp.*z,2);
            er = zi-obj.H*obj.xp;
            mat2 = er*er.';
            
            for k = 1:N
                
                zTild = z(:,k) - obj.H*obj.xp;                
                mat = zTild*zTild.';
                
                pTild = pTild+beta(k+1)*mat;
            end
            
            pTild = pTild-mat2;
            
        end
        
        function updateState(obj,xu)
            
            obj.xu = [obj.xu xu];
            
        end
    end
end
        
        
    