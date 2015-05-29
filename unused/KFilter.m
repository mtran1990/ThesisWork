classdef KFilter < handle
    %% properties
    properties (SetAccess = private)
        F
        H
        Q
        R
        xu
        Pu
    end
    
    %% methods
    methods (Access = public)
        
        function obj = KFilter(F,H,x0,P0,Q,R)
            if(nargin > 0)
                obj.F = F;
                obj.H = H;
                obj.xu = x0;
                obj.Pu = P0;
                obj.Q = Q;
                obj.R = R;
            else
                dt = 0.5;
                F = [1 dt; 0 1];
                obj.F = [F zeros(2,2); zeros(2,2) F];
            end
        end
        
        function [xu, Pu] = iterFilter(obj, meas)
            % prediction step
            xp = obj.F*obj.xu;
            Pp = obj.F*obj.Pu*obj.F.'+obj.Q;
            % update step
            y = meas-obj.H*xp;
            S = obj.H*Pp*obj.H.'+obj.R;
            K = Pp*obj.H.'/S;
            xu = xp+K*y;
            Pu = (eye(size(Pp))-K*obj.H)*Pp;
            obj.xu = xu;
            obj.Pu = Pu;
        end
        
    end
    
end