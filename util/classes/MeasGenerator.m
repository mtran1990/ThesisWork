classdef MeasGenerator < handle
    %% properties
    properties (SetAccess = private)
        Bfa
        Pd
        n0
    end
    
    methods (Access = public)        
        % constructor
        function obj = MeasGenerator(params)
            obj.Bfa = params.Bfa;
            obj.Pd = params.Pd;
            obj.n0 = params.n0;
        end
        
        function out = getMeasurements(tgtLoc, mapDim)

            clutter = getClutter(mapDim);
            meas = genMeas(tgtLoc);
            
            out = [clutter meas];
            
        end
        
    end
    
    methods (Access = private)
        
        function clutter = genClutter(obj, mapDim)
            
            x = mapDim(1);
            y = mapDim(2);
            
            m = poissrnd(obj.Bfa*(x*y));
            
            clutter = -[x;y]/2+[x 0; 0 y]*rand(2,m);
            
        end
        
        function meas = genMeas(obj, tgtLoc)
            % assume tgtLoc is 2xn, where n is the # of tgts
            [~,n] = size(tgtLoc);
            idx = rand(1,n)>obj.Pd;            
            meas = tgtLoc(:,idx);
        end
        
    end    
    
end