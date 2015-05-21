function params = genParams(tgtPath, nodeList, adjMat)

mapDims = [12 12];

% Simulation Parameters
% sParams:
% dim   :    dimension of the simulation (e.g. 2D)
% maxV  :    Max Velocity
% dt    :    time between measurements
% sigA2 :    variance of acceleration in model
% sigM2 :    variance of measurements
% n0    :    noise density
% Kap   :    gate size parameter        
% Pd    :    Probability of detection
% Pg    :    Probability of a measurement lying inside a gate
% gamG  :    Gate theshold, calculated from Pg
% Bfa   :    Density of false alarm measurements
% V     :    Volume of the simulation
dim   = 2;
maxV  = 10;
dt    = 0.1;
sigA2 = 0.1;
sigM2 = 0.1;
n0    = 0.5;
Kap   = 3;
Pd    = 0.95;
Pg    = 0.95;
gamG  = chi2inv(Pg,dim);
Bfa   = 1E-5;
V     = mapDims(1)*mapDims(2);

sParams = struct('dim',dim,'maxV',maxV,'dt',dt,'sigA2',sigA2,...
    'sigM2',sigM2,'n0',n0,'Kap',Kap,'Pd',Pd,'Pg',Pg,'gamG',gamG,...
    'Bfa',Bfa,'V',V);

params = struct('tgtPath', tgtPath, 'nodeList', nodeList, ...
    'adjMat', adjMat, 'mapDims', mapDims, 'sParams', sParams);

end