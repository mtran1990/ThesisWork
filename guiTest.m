%% GUI Test

clear;

utilDir = genpath('./util');
addpath(utilDir);

A = 1;
freq = 0.1;

% tgtPath = @(t)(circPath(A,freq,t));

% two targets
tgtPath = @(t)([circPath(A,freq,t) circPath(3*A,freq,-t)]);
%% Node Tests
a = Node;
b = Node(1,2);
c = Node(0,2);
d = Node(1,0);

nodeList = [a b c d];

adjMat = logical([0 0 1 1; 0 0 1 0; 1 1 0 1; 1 0 1 0]);

% one node
% nodeList = a;
% 
% adjMat = false;

params = genParams(tgtPath, nodeList, adjMat);

sim = SignalSimulator(0,params);

%% Measurement Tests
sim.simulate;

createMainFig(1,sim);


%% Cleanup
rmpath(utilDir);