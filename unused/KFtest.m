%% Kalman Filter Test

utilpath = genpath('./util');
addpath(utilpath);

A = 1.5;
freq = 0.2;

tgtPath = @(t)(circPath(A,freq,t));

% set up filter
dt = 0.1;
F = [1 dt; 0 1];
F = [F zeros(2,2); zeros(2,2) F];

H = [1 0 0 0; 0 0 1 0];

sigA = 0.1;
Q = [dt^4/4 dt^3/2; dt^3/2 dt^2]*sigA;
Q = [Q zeros(2); zeros(2) Q];

sigM2 = 0.01;
R = [sigM2 0; 0 sigM2];

x0 = [1.5 0 0 0]';
P0 = 10*eye(4);

KF = KFilter(F,H,x0,P0,Q,R);

% set up measurements

t = 0:dt:5;
n = length(t);

sigN = 0.1;

xPath = tgtPath(t);
xMeas = xPath+sqrt(sigN)*randn(2,n);

x = zeros(4,n);

for k = 1:n
   
    x(:,k) = KF.iterFilter(xMeas(:,k));
        
end

figure(1);
plot(xMeas(1,:),xMeas(2,:),x(1,:),x(3,:),xPath(1,:),xPath(2,:));
xlabel('Distance (m)');
ylabel('Distance (m)');
grid on;
legend('Measured Position','KF output');
axis([-3 3 -3 3]);
axis square;

rmpath(utilpath);