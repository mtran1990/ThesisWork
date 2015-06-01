
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
        
        function initTrackManager(obj,tracker,sParams)
            
            if(isempty(obj.TrackManager))
                obj.TrackManager = MNinitiator(tracker,sParams);
            end
            
        end
        
        function h = plotTracks(obj,ax,track,t)
            
            if(nargin<3)
                h = obj.TrackManager.plotTracks(ax);
            else
                h = obj.TrackManager.plotTracks(ax,track,t);
            end
            
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