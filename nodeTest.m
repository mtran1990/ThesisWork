
%% Thesis Simulation Tests

clear;

utilDir = genpath('./util');
addpath(utilDir);

A = 1.5;
freq = 0.2;

tgtPath = @(t)(circPath(A,freq,t));

%% Node Tests
a = Node;
b = Node(1,2);
c = Node(0,2);
d = Node(1,0);

nodeList = [a b c d];

adjMat = logical([0 0 1 1; 0 0 1 0; 1 1 0 1; 1 0 1 0]);

params = genParams(tgtPath, nodeList, adjMat);

sim = SignalSimulator(0,params);

sim.showMap(1);
sim.printDist;
t = linspace(0,5,100);
sim.plotPath(1,t);
sim.showConnections(1);

%% Measurement Tests
tgtLoc = tgtPath(0);

sim.measTgts(tgtLoc);
% sim.exchangeMeas;
% sim.getRawEstimates;

rmpath(utilDir);