
%% Thesis Simulation Tests

clear;

addpath(genpath('./util'));

a = Node;
b = Node(1,2);
c = Node(0,2);
d = Node(1,0);

nodeList = [a b c d];

adjMat = [0 0 1 1; 0 0 1 0; 1 1 0 1; 1 0 1 0];

mapDims = [12 12];

A = 1.5;
freq = 0.2;

tgtPath = @(t)(circPath(A,freq,t));

sim = SignalSimulator(0,tgtPath,nodeList,adjMat,mapDims);

sim.showMap(1);
sim.printDist;
t = linspace(0,5,100);
sim.plotPath(1,t);