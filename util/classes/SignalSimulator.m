
classdef SignalSimulator
    %% properties
    properties (SetAccess = private)
        
        % time
        adjMat
        tgtPath
        nodeList
        mapDims
        
    end
    
    %% public methods
    methods (Access = public)
        % constructor
        function obj = SignalSimulator(op, arg1, arg2, arg3, arg4)
            
            switch op
                case 0
                    if(nargin == 5)
                        obj.tgtPath = arg1;
                        obj.nodeList = arg2;
                        obj.adjMat = arg3;
                        obj.mapDims = arg4;
                    else
                        error('Incorrect number of arguments');
                    end
                case 1
                    if(nargin == 3)
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
            
        
        function measTgts(obj, tgtLoc)
            
            for a = obj.nodeList
                a.measureTgts(tgtLoc);
            end
            
        end
        
        function exchangeMeas(obj)
            % tells each node to exchange its measurements with its
            % neighbors
            n = length(obj.nodeList);
            
            for k=1:n
                node = obj.nodeList(k);
                
                neighbors = obj.nodeList(obj.adjMat(k,:));
                node.receiveMeas(neighbors);
                
            end
            
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
        
        function [loc] = getTgtLoc(obj, t)
            
            loc = obj.tgtPath(t);
            
        end
        
    end
    
    
end