classdef MNinitiator < handle
    %% properties
    properties (SetAccess = private)
        % To confirm a track, first N1 measurements must be inside the gate
        % and at least M2 of the next N2 measurements must be inside the
        % gate
        % Confirmed tracks will be deleted if they are not updated within Nd steps
        
        % Tracker Parameters
        % tracker:
        % tType : Type of tracker (only MN right now)
        % N1   : see above
        % M2   : see above
        % N2   : see above
        % Nd   : see above
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
        
        % 5xN vector tracking the state of N tracks
        % [s m mbar n nd]'
        % s    : deleted, tentative, or confirmed state (0, 1, or 2)
        % m    : # of measurements collected since creation
        % mbar : # of missed measurements        
        % n    : age of initiator
        % nd   : consecutively missed measurements
        % note: should I change this to a struct array instead?        
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
            % validTracks: 1xT index vector, where T is the number
            %              of valid tracks for this node, and each element
            %              represents that track's index in the full list
            % tList      : list of valid tracks in the trackList
            
            % can't return values unless tracks have been initialized
            if(~isempty(obj.trackList))
                % confirmed and tentative only for now
                % also, track must have been updated at least once, i.e.
                % age is 1 or greater
                validTracks = find(obj.trackStates(1,:) ~= 0);
                tList = obj.trackList(validTracks);
            else
                validTracks = [];
                tList = [];
            end
        end
        
        function matchTracks(obj,rTracks)
            % rTracks: received tracks from another node
            
            % right now can only run if this node has tracks to match to
            % consider copying tracks from other node if this is the case?
            % also, can only run if there is at least one valid track
            
            validTracks = obj.getValidTracks;
            
            if(~(isempty(obj.trackList) || isempty(rTracks) || ...
                    isempty(validTracks) ))
            
                % find the cost matrix
                [cMat, Na, Nb] = obj.createCostMatrix(rTracks);

                if(isempty(cMat))
                    error('here');
                end
                
                % run the JV algorithm
                rowsol = lapjv(cMat);

                % combine track data
                obj.combineTracks(rowsol,rTracks,Na,Nb);
            end
        end
        
        function updateTrackStates(obj)
            
            % looking at all the tracks, there are three categories
            % 1) Not Valid Tracks: These aren't updated at all
            % 2) Valid Tracks with measurements inside their gates: these
            %    are updated
            % 3) Valid Tracks with no measurements inside their gates
            [validTracks,tList] = obj.getValidTracks;
            
            N = length(tList);
            for k = 1:N
                updated = tList(k).iterFilter;
                obj.iterState(updated,validTracks(k));
            end
            
        end
        
        function h = plotTracks(obj,ax,hand)
            
            if(isempty(hand))
                N = size(obj.trackList,2);
                maxColor = 10;
                map = jet(maxColor);
                h = zeros(1,N);

                hold on;
                for k = 1:N
                    if(isempty(obj.trackList(k).xu))
                        x = obj.trackList(k).xp(1,:);
                        y = obj.trackList(k).xp(3,:);
                    else
                        x = obj.trackList(k).xu(1,:);
                        y = obj.trackList(k).xu(3,:);
                    end
                    idx = 1+mod(k,maxColor);
                    h(k) = plot(ax,x,y,'color',map(idx,:));
                end
                hold off;
            else
                
                str = get(hand(1),'visible');
                
                if(strcmp(str,'on'))
                    option = 'off';
                else
                    option = 'on';
                end
                
                set(hand,'visible',option);
                h = hand;
            end
        end
        
    end
    
    %% private methods
    methods (Access = private)
        
        function initTrack(obj,z)            
            
            % state vector is [x x' y y']
            x0 = [z(1); 0; z(2); 0];            
            P0 = obj.genP0mat;            
            params = obj.genKFparams(x0,P0);
            
            if(isempty(obj.trackList))
                obj.trackList = KCFilter(params);
            else
                obj.trackList = [obj.trackList KCFilter(params)];
            end            
            
        end
                
        function params = genKFparams(obj,x0,P0)
            
            params = obj.sParams;
            params.x0 = x0;
            params.P0 = P0;
            
        end
        
        function createNewTracks(obj,z)
            
            % can only create tracks if z is non-empty
            if(~isempty(z))
                % assuming N measurements are arranged as 2xN
                N = size(z,2);
                for k = 1:N
                    obj.initTrack(z(:,k));
                end

                if(isempty(obj.trackStates))
                    obj.trackStates = zeros(5,N);
                    obj.trackStates(1,:) = 1;
                else
                    tmp = zeros(5,N);
                    tmp(1,:) = 1;
                    obj.trackStates = [obj.trackStates tmp];
                end
            end
            
        end
        
        function clustAndUpdateTracks(obj,z)
            
            [validTracks,tList] = obj.getValidTracks;
            
            % need to gate the measurements
            vMat = obj.getValidationMat(z,tList);

            % cluster measurements in vMat?
            [measVec, tgtsVec, numClust, measOutClus, tgtsOutClus] = ...
                obj.clusterTgts(vMat);

            % update tracks in clusters
            obj.updateTracks(measVec,tgtsVec,tList,z,numClust,tgtsOutClus);            
            
            % create new tracks from leftover measurements
            zOut = z(:,measOutClus);
            obj.createNewTracks(zOut);
            
        end
        
        function updateTracks(obj,measVec,tgtsVec,tList,z,numClust,...
                tgtsOutClus)
            
            Pd = obj.sParams.Pd;
            Pg = obj.sParams.Pg;
            lam = obj.sParams.Bfa*obj.sParams.V;
            
            % update targets in clusters
            for k = 1:numClust
                
                zClust = z(:,measVec{k});
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
        
        function [vMat] = getValidationMat(obj,z,tList)
            
            % N: # of measurements
            N = size(z,2);
            
            % M: # of tracks (confirmed and tentative only for now)
            M = length(tList);                        
            
            vMat = zeros(N,M);
            
            for i=1:N
                for j=1:M
                    
                    % S is the covariance of the innovation vector
                    % H*P_k|k-1*H.'+R?
%                     S = tList(j).P;
                    S = tList(j).getInnCovariance;
                    yhat = tList(j).H*tList(j).xp;
                    
                    y = z(:,i);
                    
                    if(obj.insideGate(S,y,yhat))
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
            
%             % first and third elements correspond to the measurement
%             % covariances
%             S_ = S([1 3],[1 3]);
            
            % make sure S_ is symmetric
            S_ = 0.5*(S+S.');
            
            U = cholcov(S_);
            
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
                aug(logical(eye(Na))) = g;
                
                cMat = [costs aug];
            else
                
                aug = inf(Nb);
                aug(logical(eye(Nb))) = g;
                
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
        
        function iterState(obj,updated,idx)
            
            N1 = obj.tracker.N1;
            M2 = obj.tracker.M2;
            N2 = obj.tracker.N2;
            Nd = obj.tracker.Nd;
            
            % increase the age of the track by 1
            obj.trackStates(4,idx) = obj.trackStates(4,idx)+1;
            
            % if updated, increase the number of measurements collected and
            % reset the number of consecutively missed measurements
            % else, increase the number of missed measurements and add one
            % to the number of consecutively missed measurements
            if(isempty(updated) || updated)
                obj.trackStates(2,idx) = obj.trackStates(2,idx)+1;
                obj.trackStates(5,idx) = 0;
            else
                obj.trackStates(3,idx) = obj.trackStates(3,idx)+1;
                obj.trackStates(5,idx) = obj.trackStates(5,idx)+1;
            end
            
            % if we've reached the Nd missed measurements, then delete the
            % track
            if(obj.trackStates(5,idx)>= Nd)
                obj.trackStates(1,idx) = 0;
            else
                % we need N1 consecutive measurements to move from tentative
                % to confirmed
                if(obj.trackStates(1,idx) <= N1)

                    % if we missed any of the first N1 measurements, then
                    % delete the track
                    if(obj.trackStates(3,idx)>0)
                        obj.trackStates(1,idx) = 0;
                    end

                else
                    if(obj.trackStates(3,idx)>(M2+N2))
                        obj.trackStates(1,idx) = 0;
                    elseif(obj.trackStates(2,idx)>=(N1+M2))
                        obj.trackStates(1,idx) = 2;
                    else
                        obj.trackStates(1,idx) = 1;
                    end
                end
            end
            
        end
    end
        
end