classdef MeasGenerator < handle
    %% properties
    properties (SetAccess = private)
        Bfa
        V
        Pd
        n0
        x
        y
    end
    
    methods (Access = public)        
        % constructor
        function obj = MeasGenerator(params,mapDims)
            obj.Bfa = params.Bfa;
            obj.V = params.V;
            obj.Pd = params.Pd;
            obj.n0 = params.n0;
            obj.x = mapDims(1);
            obj.y = mapDims(2);
        end
        
        function out = getMeasurements(obj, tgtLoc)

            clutter = obj.genClutter;
            meas = obj.genMeas(tgtLoc);
            
            out = [clutter meas];            
            
        end
        
    end
    
    methods (Access = private)
        
        function clutter = genClutter(obj)            
            
            m = poissrnd(obj.Bfa*obj.V);
            
            clutter = -[obj.x;obj.y]/2*ones(1,m)+...
                [obj.x 0; 0 obj.y]*rand(2,m);
            
%             if(m>0)
%                 disp('check');
%             end
            
        end
        
        function meas = genMeas(obj, tgtLoc)
            % assume tgtLoc is 2xn, where n is the # of tgts
            [~,n] = size(tgtLoc);
            idx = rand(1,n)<obj.Pd;            
            meas = tgtLoc(:,idx)+obj.n0^2*randn(size(tgtLoc(:,idx)));
        end
        
    end    
    
end