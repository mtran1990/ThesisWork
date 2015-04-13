
classdef Node < handle
    %% properties
    properties (SetAccess = private)
        loc
        tgtDist
        tgtMeas
        allRanges
        sensorLoc
        rawLoc
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
            
            obj.tgtDist = sqrt(sum((tgtLoc-obj.loc).^2));
            out = obj.tgtDist;
            
        end
        
        function setTgts(obj, tgtMeas)
            obj.tgtMeas = tgtMeas;                        
        end
        
        function receiveMeas(obj, nodes)                        
            % function to fill allRanges with the measurements from the
            % current node and all its neighbors
            % takes in a list of nodes who are the current node's neighbors            
            n = length(nodes);            
            m = length(obj.tgtMeas);
            
            sensorCount = zeros(1,n);
            
            MAXMEAS = 100;
            
            while(m>MAXMEAS)
                MAXMEAS = MAXMEAS*2;
            end            
            
            obj.allRanges = zeros(1,MAXMEAS);
            obj.allRanges(1:m) = obj.tgtMeas;
            
            count = m;
            
            for k = 1:n
                tmp = nodes(k).tgtMeas;
                numMeas = length(tmp);
                
                if(count+numMeas>MAXMEAS)
                    while(count+numMeas>MAXMEAS)
                        MAXMEAS = MAXMEAS*2;
                    end
                    newArr = zeros(1,MAXMEAS);
                    newArr(1:count) = obj.allRanges;
                    obj.allRanges = newArr;
                end                    
                
                obj.allRanges(count+1:count+numMeas) = tmp;
                
                count = count+numMeas;
                
                sensorCount(k) = count;                
            end                        
            
            % trim the array
            obj.allRanges(count+1:end) = [];
            
            sensorLoc_ = zeros(2,length(obj.allRanges));
            sensorLoc_(:,1:m) = obj.loc*ones(1,m);
            
            count = m;
            
            for k = 1:n                
                endIdx = sensorCount(k);
                nodeLoc = nodes(k).loc;
                sensorLoc_(:,count+1:endIdx) = nodeLoc*ones(1,(endIdx-count));
                count = endIdx;
            end
            
            % sensorLoc is a 2x(length(allRanges))
            % column values are the location of the sensor that produced
            % the associated measurement
            obj.sensorLoc = sensorLoc_;
            
        end
        
        function estimateRawLoc(obj)
            % take allRanges and associate them with tracks
            % for now, assume all measurements come from the same target
            
            % need at least 3 measurements for algorithm to work
            if(length(obj.allRanges)>2)
                obj.rawLoc = gaussNewton(obj.sensorLoc', obj.allRanges');
            end
            
        end
        
    end
    
end