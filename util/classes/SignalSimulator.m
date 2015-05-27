
classdef SignalSimulator < handle
    %% properties
    properties (SetAccess = private)
        
        tParams
        adjMat
        tgtPath
        nodeList
        mapDims
        
        mGenerator
    end
    
    %% public methods
    methods (Access = public)
        % constructor
        function obj = SignalSimulator(op, params)
            
            switch op
                case 0
                    if(nargin == 2)
                        obj.tgtPath = params.tgtPath;
                        obj.nodeList = params.nodeList;
                        obj.adjMat = params.adjMat;
                        obj.mapDims = params.mapDims;
                        obj.tParams = params.tParams;
                        
                        obj.mGenerator = ...
                            MeasGenerator(params.sParams,params.mapDims);
                        
                        % since TrackManager of nodes aren't defined by
                        % default, need to initialize them
                        obj.initNodeTM(params.nodeList,...
                            params.tracker,params.sParams);
                    else
                        error('Incorrect number of arguments');
                    end
                case 1
                    if(nargin == 2)
                        obj.tgtPath = arg1;
                        % arg2: # of nodes
                        % use a Gaussian Distribution to set node locations
                        % calc mapDims from max node distances
                        
                    else
                        error('Incorrect number of arguments');
                    end                    
                otherwise
                    error('Invalid use of SignalSimulator constructor');
            end
        end
        
        function showMap(obj, figNum)
            figure(figNum);
            axis([-obj.mapDims(1)/2 obj.mapDims(1)/2 ...
                -obj.mapDims(2)/2 obj.mapDims(2)/2]);
            grid on;
            xlabel('Distance (m)');
            ylabel('Distance (m)');
            axis square;
            hold on;
            
            for a = obj.nodeList
            
                plot(a.loc(1), a.loc(2), 'ro');
                
            end
            hold off;
            
        end
        
        function showConnections(obj, figNum)
            
            obj.showMap(figNum);
            hold on;
            
            n = length(obj.nodeList);            
            for k = 1:n
                
                node = obj.nodeList(k);
                
                for a = obj.nodeList(obj.adjMat(k,:))
                    
                    xy = [node.loc a.loc];
                    plot(xy(1,:),xy(2,:),'k');
                    
                end
                
                
            end
            
            hold off;
            
        end
        
        function plotPath(obj, figNum, t)
            
            figure(figNum);
            axis([-obj.mapDims(1)/2 obj.mapDims(1)/2 ...
                -obj.mapDims(2)/2 obj.mapDims(2)/2]);
            grid on;
            xlabel('Distance (m)');
            ylabel('Distance (m)');
            axis square;
            hold on;
            
            loc = obj.tgtPath(t);
            
            plot(loc(1,:),loc(2,:),'g');
            
            hold off;
            
        end
        
        function printDist(obj)
            
            tgtLoc = obj.getTgtLoc(0);
                        
            fprintf('Node Distances to Target\n');
            for k = 1:length(obj.nodeList)
                tmp = obj.nodeList(k);
                fprintf('Node %2d: %.2fm\n',k,tmp.calcTgtDist(tgtLoc));
            end
        
        end
        
        function simulate(obj)
            
            while(~obj.isDone)
                obj.iterSim;
            end
            
        end
        
        function iterSim(obj)
            
            % runs one iteration of the simulation
            obj.advanceTime;
            
            loc = obj.getTgtLocNow;
            obj.measTgts(loc);
            obj.exchangeMeas;
            obj.updateNodes;
        end        
        
        function getRawEstimates(obj)
            % tells each node to use the pooled measurements to give an
            % intial estimate of the target location
            n = length(obj.nodeList);
            
            for k=1:n
                node = obj.nodeList(k);
                node.estimateRawLoc;
            end
            
        end
        
    end
        
    methods (Access = private)
        
        function initNodeTM(~,nodes,tracker,sParams)
            
            for node = nodes
                node.initTrackManager(tracker,sParams);
            end
            
        end
        
        function [loc] = getTgtLoc(obj, t)
            
            loc = obj.tgtPath(t);
            
        end
        
        function measTgts(obj, tgtLoc)
            
            for a = obj.nodeList
                meas = obj.mGenerator.getMeasurements(tgtLoc);
                a.setTgts(meas);
            end
            
        end
        
        function exchangeMeas(obj)
            % tells each node to exchange its measurements with its
            % neighbors
            n = length(obj.nodeList);
            
            for k=1:n
                node = obj.nodeList(k);
                
                neighbors = obj.nodeList(obj.adjMat(k,:));
                node.matchTracks(neighbors);
                
            end
            
        end
        
        function updateNodes(obj)
            
            for n = obj.nodeList                
                n.updateTracks;                
            end
            
        end
        
        function advanceTime(obj)
            obj.tParams.now = obj.tParams.now+obj.tParams.dt;
        end
        
        function loc = getTgtLocNow(obj)
           
            loc = obj.tgtPath(obj.tParams.now);
            
        end
        
        function done = isDone(obj)
            done = obj.tParams.now > obj.tParams.end;
        end
    end
    
    
end