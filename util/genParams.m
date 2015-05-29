function params = genParams(tgtPath, nodeList, adjMat)

mapDims = [12 12];

% Time Parameters
% tParams
% start: simulation start time
% end  : simulation end time
% dt   : time between measurements (Note: also used below)
% now  : current time
s = 0;
e = 10;
dt = 0.1;
n = s;

tParams = struct('start',s,'end',e,'dt',dt,'now',n);

% Right now I'm only working on the M-N tracker
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

tType = 'MN';
N1    = 2;
M2    = 2;
N2    = 3;
Nd    = 5;

tracker = struct('tType',tType,'N1',N1,'M2',M2,'N2',N2,'Nd',Nd);

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
dt    = dt;
sigA2 = 0.5;
sigM2 = 0.1;
n0    = 0.1;
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
    'adjMat', adjMat, 'mapDims', mapDims, 'tParams',tParams,...
    'sParams', sParams, 'tracker',tracker);

end