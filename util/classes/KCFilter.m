% Kalman Consensus Filter
classdef KCFilter < handle
    %% properties
    properties (SetAccess = private)
        % Model Parameters
        % modelParams:
        % F: state transition matrix
        % H: state vector to measurement matrix (observation model)
        % Q: model noise covariance
        % R: measurement noise covariance
        modelParams

        % State Parameters
        % stateParams:
        % xu: 4xN   updated state vector
        % xp: 4xN   predicted state
        % P : 4x4xN covariance matrix
        % t : 1xN   time of update
        stateParams
        
        % information vector / matrix
        % the rest are temporary variables used when updating the filter
        u
        U
        xdiff
        beta0
        pTild
                
        % Simulation Parameters
        % sParams:
        % dim   :    dimension of the simulation (e.g. 2D)
        % maxV  :    Max Velocity
        % dt    :    time between measurements
        % sigA2 :    std. dev. of acceleration in model
        % sigM2 :    variance of measurements
        % n0    :    noise density
        % Kap   :    gate size parameter        
        % Pd    :    Probability of detection
        % Pg    :    Probability of a measurement lying inside a gate
        % gamG  :    Gate theshold, calculated from Pg
        % Bfa   :    Density of false alarm measurements
        % V     :    Volume of the simulation        
        sParams
    end
    
    %% methods
    methods (Access = public)
        
        function obj = KCFilter(sParams,x0,P0)
            
            obj.sParams = sParams;            
            obj.initModel;            
            obj.initState(x0,P0);

        end
        
        function addMeasurement(obj, z)
            
            % can't calculate weights without measurements
            if(nargin == 2)
                
                % get the model parameters
                [~,H,~,R] = obj.getModelParams;
                
                beta = obj.calculateWeights(z);
                obj.beta0 = beta(1);
                obj.pTild = obj.calcPTild(beta,z);
                
                % get a matrix of betas to elementwise multiply with
                beta_tmp = ones(2,1)*beta;

                % don't need the first element it seems (corresponds to
                % false alarm)
                % need to sum across rows
                zVec = sum(beta_tmp(:,2:end).*z,2);

                obj.u = H.'*(R\zVec);
                obj.U = H.'*(R\H);
                
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

                % get the state variables
                [~,xp ,~] = obj.getState;
                [~,xp2,~] = T2.getState;
                
                obj.xdiff = obj.xdiff + (xp-xp2);
            end
        end
        
        function updated = iterFilter(obj,t)
            
            % updated: true if filter updated with measurements
            %          false if filter updated with no measurements
            %          empty if filter did not update
            
            % if u is empty, track hasn't updated yet, so don't do anything
            % yet
            if(~isempty(obj.u))
            
                % get the model parameters
                [F,H,Q,R] = obj.getModelParams;
                
                % get the state variables
                [~,xp,P] = obj.getState;
                
                % if u is a column vectors of zeros, then the track didn't
                % receive any measurements
                if(isequal(zeros(4,1),obj.u))
                    % update step
                    xu = xp;

                    % prediction step
                    % NOTE: not sure if this part is right, need to double
                    % check
                    P = F*P*F.'+Q;
                    xp = F*xu;
                    updated = false;
                else            
                    eps = 0.01;
                    y = obj.u;
                    S = obj.U;

                    S0 = (1-obj.beta0)*S;

                    % update step
                    Mtild = (inv(P)+S);
                    K = Mtild\H.'/R;
                    M = obj.beta0*P+...
                        (1-obj.beta0)*eye(size(Mtild))/Mtild+...
                        K*obj.pTild*K.';

                    xu = xp + Mtild\(y-S0*xp)+...
                        eps*M/(1+norm(M,'fro'))*obj.xdiff;                

                    % prediction step
                    P = F*M*F.'+Q;
                    xp = F*xu;
                    updated = true;
                    
                end

                obj.updateState(xu,xp,P,t);
            
            else
                
                updated = [];
            
            end
            
        end
        
        function cost = computeCost(obj,T2)
            % computes the cost between the current track and the given
            % track (using Mahalanobis distance)
            
            % get the state variables
            [~,xp ,P ]  = obj.getState;
            [~,xp2,P2] = T2.getState;
            
            S = (P+P2);
            xdiff_ = xp-xp2;
            c = xdiff_.'*S*xdiff_;
            
            cost = c+log(det(S));
        end        
        
        function inside = checkGate(obj,y)
            
            % get the model parameters
            [~,H,~,~] = obj.getModelParams;
            
            % get the state variables
            [~,xp,~]  = obj.getState;
            
            S = obj.getInnCovariance;
            
            U_ = obj.factorS(S);
            
            yhat = H*xp;
            
            testStat = norm( (U_.'\(y-yhat)), 2);
            
            if testStat<= obj.sParams.gamG
                inside = true;
            else
                inside = false;
            end
            
        end
        
        function U = factorS(~,S)
            % make sure S_ is symmetric
            S_ = 0.5*(S+S.');
            
            U = cholcov(S_);
        end
        
        function U = getEllipsoidMat(obj,P)
            
            % get the model parameters
            [~,H,~,R] = obj.getModelParams;
            
            S = H*P*H.'+R;
            
            U = obj.factorS(S);
            
        end
        
        function [xu,xp,P,t] = getState(obj,getAll)
            
            % usually only the most recent values are needed
            % if getAll is true, then, pull out all the state parameters
            if(nargin < 2 || ~getAll)
                if(~isempty(obj.stateParams.xu))
                    xu = obj.stateParams.xu(:,end);
                    t  = obj.stateParams.t(end);
                else
                    xu = [];
                    t  = [];
                end
                xp = obj.stateParams.xp(:,end);
                P  = obj.stateParams.P(:,:,end);
                
            else
                xu = obj.stateParams.xu;
                xp = obj.stateParams.xp;
                P  = obj.stateParams.P;
                t  = obj.stateParams.t;
            end
        end
        
    end
    
    %% private methods
    methods (Access = private)
        
        function initModel(obj)
            
            dt = obj.sParams.dt;
            F_ = [1 dt; 0 1];
            F = [F_ zeros(2,2); zeros(2,2) F_];

            H = [1 0 0 0; 0 0 1 0];
            
            
            sigA2 = obj.sParams.sigA2;
            sigM2 = obj.sParams.sigM2;
            
            Q_ = [dt^4/4 dt^3/2; dt^3/2 dt^2]*sigA2;
            Q = [Q_ zeros(2); zeros(2) Q_];
            
            R = [sigM2 0; 0 sigM2];            
            
            obj.modelParams = struct('F',F,'H',H,'Q',Q,'R',R);
            
        end
        
        function initState(obj,x0,P0)
            
            obj.stateParams = struct('xu',[],'xp',x0,'P',P0,'t',[]);
            
        end
        
        function [F,H,Q,R] = getModelParams(obj)
            
            F = obj.modelParams.F;
            H = obj.modelParams.H;
            Q = obj.modelParams.Q;
            R = obj.modelParams.R;
            
        end                
        
        function Beta = calculateWeights(obj, z)
            % z  : 2xN vector of N measurements
            % Pd : Probability of Detection
            % Pg : Probability of measurement lying inside the gate
            % lam: Mean of the poisson distribution of false alarms            
            
            % get the model parameters
            [~,H,~,R] = obj.getModelParams;
            
            % get the state variables
            [~,xp,P] = obj.getState;
            
            % get the necessary simulation parameters
            Pd = obj.sParams.Pd;
            Pg = obj.sParams.Pg;
            lam = obj.sParams.Bfa*obj.sParams.V;
            
            % need to convert tmp to be the same size as zTild
            tmp = H*xp;
            tmp = tmp*ones(1,size(z,2));
            
            zTild = z - tmp;
            
            % innovation covariance (covariance of zTild?)
            innCov = H*P*H.'+R;
            
            % dimension of target
            dim = obj.sParams.dim;
            
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
            
            % get the model parameters
            [~,H,~,~] = obj.getModelParams;
            
            % get the state variables
            [~,xp,~] = obj.getState;
            
            [M,N] = size(z);
            tmp = zeros(M);
            
            y = H*xp;
            
            beta_tmp = ones(2,1)*beta(2:end);
            
            % Note: there was an error in the paper by Sandell I think, I'm
            % using the equation from Fortmann instead here
            % er needs to be the same size as z, so multiply "y" by 1xN
            % matrix
            er = z-y*ones(1,N);
            zi = sum(beta_tmp.*er,2);
            mat2 = zi*zi.';
            
            for k = 1:N
                
                zTild = z(:,k) - y;                
                mat = zTild*zTild.';
                
                tmp = tmp+beta(k+1)*mat;
            end
            
            pTild = tmp-mat2;
            
        end
        
        function updateState(obj,xu,xp,P,t)
            
            obj.stateParams.xu(:,end+1) = xu;
            obj.stateParams.xp(:,end+1) = xp;
            obj.stateParams.P(:,:,end+1) = P;
            obj.stateParams.t(:,end+1) = t;
            
        end
        
        function S = getInnCovariance(obj)
            
            % get the model parameters
            [~,H,~,R] = obj.getModelParams;
            
            % get the state variables
            [~,~,P] = obj.getState;
            
            % S is the covariance of the innovation vector
            % H*P_k|k-1*H.'+R?
            S = H*P*H.'+R;
            
        end
        
    end
end
        
        
    