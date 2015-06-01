
classdef SignalSimulator < handle
    %% properties
    properties (SetAccess = private)
        
        tParams
        adjMat
        tgtPath
        nodeList
        mapDims
        
        mGenerator
        
        % gui info
        % gui:
        % nodeHandles: handles for the location of the nodes on the map
        % connHandles: handles for the lines between nodes
        % poinHandles: handles for the point view
        gui
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
                        
                        obj.initGUIstruct;
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
        
        function showMap(obj, ax)
            
            if(isempty(obj.gui.nodeHandles))
                N = length(obj.nodeList);
                h = zeros(1,N);
                hold on;
                for k = 1:N
                    
                    a = obj.nodeList(k);
                    h(k) = plot(ax,a.loc(1), a.loc(2), 'ro');
                    
                end
                hold off;
                
                obj.gui.nodeHandles = h;
                
            else
                % else toggle the visibility
                str = get(obj.gui.nodeHandles(1),'visible');
                
                if(strcmp(str,'on'))
                    option = 'off';
                else
                    option = 'on';
                end
                
                set(obj.gui.nodeHandles,'visible',option);
                
            end
            
        end
        
        function showConnections(obj, ax)
            
            if(isempty(obj.gui.connHandles))

                N = sum(obj.adjMat(:));
                h = zeros(1,N);
                i = 1;
                
                hold on;
                n = length(obj.nodeList);
                for k = 1:n
                    
                    node = obj.nodeList(k);
                    
                    for a = obj.nodeList(obj.adjMat(k,:))
                        
                        xy = [node.loc a.loc];
                        h(i) = plot(ax,xy(1,:),xy(2,:),'k');
                        i = i+1;
                    end
                    
                    
                end
                hold off;
                
                obj.gui.connHandles = h;
                
            else
                % else toggle the visibility
                str = get(obj.gui.connHandles(1),'visible');
                
                if(strcmp(str,'on'))
                    option = 'off';
                else
                    option = 'on';
                end
                
                set(obj.gui.connHandles,'visible',option);
                
            end
            
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
            t = obj.advanceTime;
            
            loc = obj.getTgtLocNow;
            obj.measTgts(loc);
            obj.exchangeMeas;
            obj.updateNodes(t);
        end        
        
        function initMap(obj,ax)
           
            axis(ax,[-obj.mapDims(1)/2 obj.mapDims(1)/2 ...
                -obj.mapDims(2)/2 obj.mapDims(2)/2]);
            grid on;
            xlabel('Distance (m)');
            ylabel('Distance (m)');
            axis square;
            
        end
        
        function showTracks(obj,ax,node,track,t)
            % two ways to call the function, either plot the entire track
            % or piece by piece
            
            if(nargin == 3)
                
                h = obj.gui.trackHandles{node};
                
                if(isempty(h))
                    obj.gui.trackHandles{node} = ...
                        obj.nodeList(node).plotTracks(ax);
                else
                    str = get(h(1),'visible');
                    
                    if(strcmp(str,'on'))
                        option = 'off';
                    else
                        option = 'on';
                    end
                    
                    set(h,'visible',option);
                end
                
            elseif(nargin > 3)
                
                if(~isempty(obj.gui.poinHandles))
                    delete(obj.gui.poinHandles);
                    obj.gui.poinHandles = [];
                end
                obj.gui.poinHandles = ...
                    obj.nodeList(node).plotTracks(ax,track,t);
                
            end
                
            
        end
        
        function hidePlots(obj)
            
            fields = fieldnames(obj.gui);
            
            for k = 1:numel(fields)
                
                h = obj.gui.(fields{k});
                if(~isempty(h))
                    if(iscell(h))

                        for j = 1:numel(h)

                            h2 = h{j};
                            set(h2,'visible','off');

                        end

                    else
                        set(h,'visible','off');
                    end
                end
                
            end
            
        end
        
        function N = getNumNodes(obj)
            N = length(obj.nodeList);
        end
        
        function N = getNumTracks(obj,idx)
            
            N = obj.nodeList(idx).getNumTracks;
            
        end
        
        function info = getTrackInfo(obj,node,track)
            
            info = obj.nodeList(node).getTrackInfo(track);
            
        end
        
        function loc = getNodeLoc(obj,num)
            
            if(nargin<2)
                loc = [obj.nodeList.loc];
            else
                loc = obj.nodeList(num).loc;
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
        
        function initNodeTM(~,nodes,tracker,sParams)
            
            for node = nodes
                node.initTrackManager(tracker,sParams);
            end
            
        end
        
        function initGUIstruct(obj)
            
            obj.gui = struct('nodeHandles',[],'connHandles',[],...
                'trackHandles',[],'poinHandles',[]);
            obj.gui.trackHandles = cell(1,obj.getNumNodes);
            
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
        
        function updateNodes(obj,t)
            
            for n = obj.nodeList                
                n.updateTracks(t);                
            end
            
        end
        
        function t = advanceTime(obj)
            % t: the new current time
            t = obj.tParams.now+obj.tParams.dt;
            obj.tParams.now = t;
        end
        
        function loc = getTgtLocNow(obj)
           
            loc = obj.tgtPath(obj.tParams.now);
            
        end
        
        function done = isDone(obj)
            done = obj.tParams.now > obj.tParams.end;
        end
        
    end
    
    
end