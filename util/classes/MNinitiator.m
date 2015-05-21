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
        % dim   :    dimension of the simulation (e.g. 2D)
        % maxV  :    Max Velocity
        % dt    :    time between measurements
        % sigA2 :    std. dev. of acceleration in model
        % sigM2 :    variance of measurements
        % n0    :    noise density
        % Kap   :    gate size parameter        
        % Pd    :    Probability of detection
        % Pg    :    Probability of a measurement lying inside a gate
        % gamG  :    Gate theshold, calculated from Pg
        % Bfa   :    Density of false alarm measurements
        % V     :    Volume of the simulation
        sParams
        
        % list of all tracks
        trackList
        
        % 4xN vector tracking the state of N tracks
        % [s m mbar n]'
        % s    : deleted, tentative, or confirmed state (0, 1, or 2)
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
            
            obj.gamG = chi2inv(sParams.Pg,sParams.dim);
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
                [vMat, tmpList] = getValidationMat(z);
                
                % cluster measurements in vMat?
                [measVec, tgtsVec, numClust, measOutClus] = ...
                    obj.clusterTgts(vMat);
                
                % TODO: update tracks based on clusters
                %       create new tracks from rest of measurements
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
        
        function updateTracks(obj,measVec,tgtsVec,tList,z,numClust)
            
            for k = 1:numClust
                
                zClust = z(measVec{k});
                tgts = tList(tgtsVec{k});
                
            end
            
            
            
        end
        
        function P0 = genP0mat(obj)
            
            % following the format given in Ozkan, Lecture 2 for single
            % point initiation
            sigA2 = obj.sParams.sigA2;
            Kap = obj.sParams.Kap;
            maxV = obj.sParams.maxV;
            
            sigV2 = (maxV/Kap)^2;
            
            P0 = diag([sigA2 sigV2 sigA2 sigV2]);
        end
        
        function params = genKFparams(obj,x0,P0)
            
            params = obj.sParams;
            params.x0 = x0;
            params.P0 = P0;
            
        end
        
        function [vMat, tmpList] = getValidationMat(obj,z)
            
            % N: # of measurements
            N = size(z,2);
            
            % M: # of tracks (confirmed and tentative for now)
            validTracks = obj.trackStates(1,:) ~= 0;
            M = sum(validTracks);
            
            vMat = zeros(N,M);
            
            tmpList = obj.trackList(validTracks);
            
            for i=1:N
                for j=1:M
                    
                    S = tmpList(j).P;
                    yhat = tmpList(j).H*tmpList(j).xp;
                    
                    y = z(:,i);
                    
                    if(insideGate(S,y,yhat))
                        vMat(i,j) = 1;
                    end
                    
                end
            end
            
        end
        
        % might consider trying to optimize these for loops later
        %{
        function yhat = getEstimatesFromTrack(~,list)
            
            % assuming all H's are the same (all sensors the same)
            H = list(1).H;
            
            yhat = H*[list.xp];
            
            % Uncomment if H's are different
            %{
            M = length(list);
            yhat = zeros(2,M);
            for k = 1:M
                
                yhat(:,k) = list(k).H*list(k).xp;
                
            end
            %}
            
        end
        %}
        
        function inside = insideGate(obj,S,y,yhat)
            
            U = cholcov(S);
            
            testStat = norm( (U.'\(y-yhat)), 2);
            
            if testStat<= obj.sParams.gamG
                inside = true;
            else
                inside = false;
            end
            
        end
        
        % returns a cell array of tgts separated into clusters
        function [measInClus, tgtsInClus, numClust, measOutClus, vCopy] ...
                = clusterTgts(~, vMat)
            
            vCopy = vMat;
            
            [m,n] = size(vMat);
            measInClus = num2cell(1:m);
                                    
            % steps: remove null rows and columns
            % use OR operator to cluster
            zeroRows = sum(vCopy,2)==0;
            % measOutClus = find(zeroRows == 1);
            measOutClus = zeroRows;
            measInClus(zeroRows) = [];
                        
            zeroCols = sum(vCopy,1)==0;
            
            vCopy(zeroRows,:) = [];
            vCopy(:,zeroCols) = [];
            
            % starting from first target, find other measurements in the
            % measurement gate, then combine those rows
            for k = 1:size(vCopy,2)
                
                % search through other measurements in the same column
                % (target)
                idx = find(vCopy(:,k)==1);
                
                % if there is more than one index, then there are multiple
                % measurements in the measurement gate
                if(length(idx)>1)
                    
                    % i1: the first measurement in the current
                    % target's measurement gate
                    i1 = idx(1);
                    
                    % loop over all other measurements
                    for j = 2:length(idx)
                        
                        inew = idx(j);
                        % combine rows
                        vCopy(i1,:) = or(vCopy(i1,:),vCopy(inew,:));
                        
                        % add measurements to current cluster
                        measInClus{i1} = [measInClus{i1} ...
                            measInClus{inew}];
                        
                        % delete old row in vCopy and measInClus
                        vCopy(inew,:) = [];
                        measInClus(inew) = [];
                    end
                    
                end
                                
            end
            
            % pull out targets in each cluster
            % we have n targets, but we have removed columns of all zeros
            tgtNums = 1:n;
            tgtNums = tgtNums(~zeroCols);
            
            numClust = size(vCopy,1);
            
            % in the final vCopy, each row is a cluster
            % 1's indicate that the target of the associated column is part
            % of that cluster
            tgtsInClus = cell(1,numClust);            
            for k = 1:numClust
                tgtsInClus{k} = tgtNums(logical(vCopy(k,:)));
            end
            
        end
        
    end
        
end