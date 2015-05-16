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
        
        function addMeasurement(obj, z)
            
            obj.u = obj.H.'*obj.R\z;
            obj.U = obj.H.'*obj.R\obj.H;
            
        end
        
        function iterFilter(obj, xj, u, U)
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
    
end
        
        
    