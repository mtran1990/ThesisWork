
classdef Node < handle
    %% properties
    properties (SetAccess = private)
        loc
        tgtDist
        TrackManager
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
            obj.TrackManager.addMeasurement(tgtMeas);
        end
        
        function [validTracks,tList] = getValidTracks(obj)
            % validTracks: 1xT index vector, where T is the number
            %              of valid tracks for this node, and each element
            %              represents that track's index in the full list
            % tList      : list of valid tracks in the trackList
            
            % (confirmed and tentative only for now)
            [validTracks, tList] = obj.TrackManager.getValidTracks;
            
        end
        
        function matchTracks(obj, nList)
            % nList: a list of nodes connected with the current node
            for n = nList
               
                [~,rTracks] = n.getValidTracks;
                obj.TrackManager.matchTracks(rTracks);
                
            end                        
            
        end
        
        function updateTracks(obj,t)
            
            obj.TrackManager.updateTrackStates(t);
            
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
        
        function initTrackManager(obj,tracker,sParams)
            
            if(isempty(obj.TrackManager))
                obj.TrackManager = MNinitiator(tracker,sParams);
            end
            
        end
        
        function h = plotTracks(obj,ax)
            
            h = obj.TrackManager.plotTracks(ax);
            
        end
        
        function N = getNumTracks(obj)
            
            N = size(obj.TrackManager.trackStates,2);
            
        end
        
        function info = getTrackInfo(obj,track)
            
            [xu,xp,P,t] = ...
                obj.TrackManager.trackList(track).getState(true);
            
            info = struct('xu',xu,'xp',xp,'P',P,'t',t);
            
        end
        
    end
    
end