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
        % sParams:
        % maxV  :    Max Velocity
        % dt    :    time between measurements
        % sigA2 :    std. dev. of acceleration in model
        % sigM2 :    variance of measurements
        % Kap   :    gate size parameter
        sParams
        
        % list of all tracks
        trackList
        
        % 4xN vector tracking the state of N tracks
        % [s m mbar n]'
        % s    : deleted, tentative, or confirmed state
        % m    : # of measurements collected since creation
        % mbar : # of missed measurements
        % n    : age of initiator
        trackStates
    end
    
    %% methods
    methods (Access = public)
        
        function obj = MNinitiator(N1,M2,N2,sParams)
            obj.N1 = N1;
            obj.M2 = M2;
            obj.N2 = N2;
            obj.sParams = sParams;
        end
        
        function addMeasurement(obj, z)
            
            if(isempty(obj.trackList))
                % assuming N measurements are arranged as 2xN
                N = size(z,2);
                for k = 1:N
                    initTrack(z(:,k));
                end
                
                obj.trackStates = ones(4,N);
                obj.trackStates(3,:) = 0;
                
            else
                % need to gate the measurements
                
            end
            
            
        end
        
    end
    
    %% private methods
    methods (Access = private)
        
        function initTrack(obj,z)
            k = length(obj.trackList)+1;
            
            % state vector is [x x' y y']
            x0 = [z(1); 0; z(2); 0];
            
            P0 = obj.genP0mat;
            
            params = genKFparams(x0,P0);
            
            obj.trackList(k) = KCFilter(params);
            
        end       
        
        function P0 = genP0mat(obj)
            
            % following the format given in Ozkan, Lecture 2 for single
            % point initiation
            sigA2 = obj.sParams.sigA2;
            Kap = obj.sParams.Kap;
            maxV = obj.sParams.maxV;
            
            P0 = diag([sigA2 maxV/Kap sigA2 maxV/Kap]);
        end
        
        function params = genKFparams(obj,x0,P0)
            
            params = obj.sParams;
            params = rmfield(params,{'maxV','Kap'});
            params.x0 = x0;
            params.P0 = P0;
            
        end
        
    end
        
end