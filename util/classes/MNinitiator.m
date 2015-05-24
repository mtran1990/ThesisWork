classdef MNinitiator < handle
    %% properties
    properties (SetAccess = private)
        % To confirm a track, first N1 measurements must be inside the gate
        % and at least M2 of the next N2 measurements must be inside the
        % gate
        
        % Tracker Parameters
        % tracker:
        % tType : Type of tracker (only MN right now)
        % N1   : see above
        % M2   : see above
        % N2   : see above
        tracker
        
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
        
        function obj = MNinitiator(tracker,sParams)
            obj.tracker = tracker;
            obj.sParams = sParams;
        end
        
        function addMeasurement(obj, z)
            
            if(isempty(obj.trackList))
                
                % if there are no tracks currently, create new tracks from
                % the measurements
                obj.createNewTracks(z);
                
            else
                
                obj.clustAndUpdateTracks(z);
                
            end
                        
        end
        
        function [validTracks,tList] = getValidTracks(obj)
            % validTracks: 1xT logical vector, where T is the total number
            %              of current tracks for this node
            % tList      : list of valid tracks in the trackList
            
            % (confirmed and tentative only for now)
            validTracks = obj.trackStates(1,:) ~= 0;
            tList = obj.trackList(validTracks);
            
        end
        
        function matchTracks(obj,rTracks)
            % rTracks: received tracks from another node
            
            % find the cost matrix
            [cMat, Na, Nb] = obj.createCostMatrix(rTracks);
            
            % run the JV algorithm
            rowsol = lapjv(cMat);
            
            % combine track data
            obj.combineTracks(rowsol,rTracks,Na,Nb);
        end
        
        function updateTrackStates(obj)
            
            % looking at all the tracks, there are three categories
            % 1) Not Valid Tracks: These aren't updated at all
            % 2) Valid Tracks with measurements inside their gates: these
            %    are updated
            % 3) Valid Tracks with no measurements inside their gates
            [validTracks,tList] = getValidTracks;
            
            N = length(tList);
            for k = 1:N
                
                
                
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
                
        function params = genKFparams(obj,x0,P0)
            
            params = obj.sParams;
            params.x0 = x0;
            params.P0 = P0;
            
        end
        
        function createNewTracks(obj,z)
            % assuming N measurements are arranged as 2xN
            N = size(z,2);
            for k = 1:N
                initTrack(z(:,k));
            end
            
            if(isempty(obj.trackStates))
                obj.trackStates = ones(4,N);
                obj.trackStates(3,:) = 0;
            else
                tmp = ones(4,N);
                tmp(3,:) = 0;
                obj.trackStates = [obj.trackStates tmp];
            end
            
        end
        
        function clustAndUpdateTracks(obj,z)
            
            [validTracks,tList] = getValidTracks;
            
            % need to gate the measurements
            vMat = obj.getValidationMat(z,tList);

            % cluster measurements in vMat?
            [measVec, tgtsVec, numClust, measOutClus, tgtsOutClus] = ...
                obj.clusterTgts(vMat);

            % update tracks in clusters
            obj.updateTracks(measVec,tgtsVec,tList,z,numClust);
            
            % TODO: update track states
%             obj.updateTrackStates(validTracks,tgtsOutClus);
            
            % create new tracks from leftover measurements
            zOut = z(measOutClus);            
            obj.createNewTracks(zOut);
            
        end
        
        function updateTracks(obj,measVec,tgtsVec,tList,z,numClust,...
                tgtsOutClus)
            
            Pd = obj.sParams.Pd;
            Pg = obj.sParams.Pg;
            lam = obj.sParams.Bfa*obj.sParams.V;
            
            % update targets in clusters
            for k = 1:numClust
                
                zClust = z(measVec{k});
                tgts = tList(tgtsVec{k});
                
                for t = tgts                    
                    t.addMeasurement(zClust,Pd,Pg,lam);                    
                end
                
            end
            
            % update targets with no measurements
            for t = tList(tgtsOutClus)
                t.addMeasurement;
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
        
        function [vMat] = getValidationMat(~,z,tList)
            
            % N: # of measurements
            N = size(z,2);
            
            % M: # of tracks (confirmed and tentative only for now)
            M = length(tList);                        
            
            vMat = zeros(N,M);
            
            for i=1:N
                for j=1:M
                    
                    S = tList(j).P;
                    yhat = tList(j).H*tList(j).xp;
                    
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
        function [measInClus, tgtsInClus, numClust, measOutClus, ...
                tgtsOutClus, vCopy] = clusterTgts(~, vMat)
            % Outputs:
            % measInClus : 1xC cell array, where M is the number of
            %              clusters, contains indicies of measurements in
            %              the cluster
            % tgtsInClus : 1xC cell array, where M is the number of
            %              clusters, contains indicies of targets in the
            %              cluster
            % numClust   : C, the number of clusters
            % measOutClus: 1xN logical vector, where N is the number of
            %              total measurements, takes on values of 1 where
            %              measurements are not in any validation gate, and
            %              0 when the associated measure is inside a gate
            % tgtsOutClus: 1xM logical vector, where M is the number of
            %              valid targets, takes on values of 1 where
            %              targets have no measurements inside their
            %              validation gates, and 0 when they have a
            %              measurement inside their validation gate            
            
            vCopy = vMat;
            
            [m,n] = size(vMat);
            measInClus = num2cell(1:m);
                                    
            % steps: remove null rows and columns
            % use OR operator to cluster
            zeroRows = sum(vCopy,2)==0;            
            measOutClus = zeroRows;
            measInClus(zeroRows) = [];
                        
            zeroCols = sum(vCopy,1)==0;
            tgtsOutClus = zeroCols;
            
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
        
        function [cMat, Na, Nb] = createCostMatrix(obj, rTracks)
            % rTracks: received tracks from another node
            
            [~,currTracks] = obj.getValidTracks;
            
            Na = length(currTracks);
            Nb = length(rTracks);
            
            costs = zeros(Na,Nb);
            
            for a = 1:Na
                for b = 1:Nb
                    
                    costs(a,b) = currTracks(a).computeCost(rTracks(b));
                    
                end
            end            
            
            % calculate costs for false alarms, right now all the sensors
            % are assumed to be the same, with the same false alarm
            % densities
            g = obj.calcFAcost;
            
            % augment the matrix
            if(Na<Nb)
                
                aug = inf(Na);
                aug(eye(Na)) = g;
                
                cMat = [costs aug];
            else
                
                aug = inf(Nb);
                aug(eye(Nb)) = g;
                
                cMat = [costs; aug];
            end
            
        end
        
        function g = calcFAcost(obj)
            
            Pd  = obj.sParams.Pd;
            Bfa = obj.sParams.Bfa;
            V   = obj.sParams.V;
            dim = obj.sParams.dim;
            
            Pnta = V*Pd*(1-Pd)+Bfa;
            Pntb = V*Pd*(1-Pd)+Bfa;
            
            g = 2*log(V*Pd^2/( (2*pi)^(dim/2)*Pnta*Pntb ));
            
        end
        
        function combineTracks(obj,rowsol,rTracks,Na,Nb)
            
            % rowsol: pairing given by JV algorithm
            % Na    : # of tracks in current node
            % Nb    : # of tracks in other node
            
            % function explanation
            %{
            % this function basically takes the pairings in rowsol and
            % combines the information in the tracks
            
            % I'll try to explain the Na<Nb case here
            % In this case, we've augmented our cost matrix to be a 
            % (Na)x(Nb+Na) size matrix
            % rowsol will have length Na
            
            % The first element of rowsol corresponds to the first track in
            % currTracks below (which is a list of all the valid tracks).
            
            % The second element of rowsol refers to the second track in
            % the same list and so on.
            
            % The value of rowsol determines which track from the other
            % node JV found a pairing with
            
            % If this value is larger than Nb, then we know it's pointing
            % to one of the false alarms
            
            % Otherwise, we use this information to fuse the tracks
            
            % Note: There might be better way to do this instead of two if
            %       statements, I'll try to think of one later
            %}
            
            [~,currTracks] = obj.getValidTracks;
            
            % two cases
            if(Na<Nb)
                
                % length of rowsol is the same as the smaller value, i.e.
                % Na in this case
                for k = 1:Na
                    
                    % if the value in rowsol is less than or equal to Nb,
                    % then it isn't pointing to a false alarm, and points
                    % to a valid track from the other node
                    if(rowsol(k)<=Nb)
                        
                        T1 = currTracks(k);
                        T2 = rTracks(rowsol(k));                        
                        
                        T1.addTrack(T2);
                        
                    end
                end
                
            else
                
                % length of rowsol is Nb
                for k = 1:Nb
                    
                    if(rowsol(k)<=Na)
                        
                        T1 = currTracks(rowsol(k));                        
                        T2 = rTracks(k);
                        
                        T1.addTrack(T2);
                        
                    end
                    
                end
                
            end
            
        end
    end
        
end