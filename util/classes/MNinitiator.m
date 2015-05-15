classdef MNinitiator < handle
    %% properties
    properties (SetAccess = private)
        % To confirm a track, first N1 measurements must be inside the gate
        % and at least M2 of the next N2 measurements must be inside the
        % gate
        N1 
        M2
        N2
        
        % Simulation Parameters        
        maxV
        dt
        
        % list of all tracks
        trackList
    end
    
    %% methods
    methods (Access = public)
        
        function obj = MNinitiator(N1,M2,N2,maxV)
            obj.N1 = N1;
            obj.M2 = M2;
            obj.N2 = N2;
            obj.maxV = maxV;
        end
        
        function addMeasurement(obj, z)
            
            if(isempty(obj.trackList))
                % assuming N measurements are arranged as 2xN
                for k = 1:size(z,2)
                    initTrack(z(:,k));
                end
            end
            
            
        end
        
    end
    
    %% private methods
    methods (Access = private)
        
        function initTrack(obj,z)
            k = length(obj.trackList)+1;
            
            % state vector is [x x' y y']
            x0 = [z(1); 0; z(2); 0];
            
            obj.trackList(k) = KCFilter(obj.dt,x0,P0);
            
        end
        
    end
    
end