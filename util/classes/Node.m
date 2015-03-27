
classdef Node < handle
    %% properties
    properties (SetAccess = private)
        loc
        tgtDist
    end
    
    %% methods
    methods (Access = public)
        % constructor
        function obj = Node(x,y)
            
            if(nargin>0)
                obj.loc = [x; y];
            else
                obj.loc = [0; 0];
            end
        end
        
        function out = calcTgtDist(obj, tgtLoc)
            
            obj.tgtDist = sqrt(sum(tgtLoc-obj.loc).^2);
            out = obj.tgtDist;
            
        end
        
    end
    
end