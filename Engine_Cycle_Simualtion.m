
%  __________ Engine Cycle Simulation with Wiebe Heat Release and Iterative Matching

clear; 
clc;
close all

%% __________ Fuel & mixture properties

qc = 44.510e3;      % Fuel LHV [kJ/kg fuel]
A  = 15.27;         % Air/fuel ratio
F  = 1/A;           % Fuel/air ratio
as = 11.25;         % Stoichiometric coefficient
Mf = 101.21;        % Fuel molecular weight [kg/kmol]

%%  __________  Initial guesses / tolerances 

ts = -19;           % Start of combustion [deg CA]
td = 40;            % Combustion duration [deg CA] 
Te = 800;           % Exhaust gas temp guess [K]
Pi = 0.2;           % Intake pressure [bar]
Pe = 1.1;           % Exhaust/back pressure [bar]
f1 = 0.01;          % Residual-burn fraction at IVC (mixing factor)
ae = 5e-5;          % Convergence tolerance (relative)
maxIter = 10000;    % Safety cap
r1 = 1; 
r2 = 1; 
r3 = 1;

%%  __________ Engine / geometry parameters 

N     = 2000;       % rpm
a     = 5;          % Wiebe a
g     = 1.24;       % gamma (const)
mfmep = 1.2;        % Friction MEP [bar]
bmep  = 3;          % Target brake MEP [bar] 
nWie  = 3;          % Wiebe exponent
rc    = 10;         % Compression ratio
Mair  = 28.85;      % kg/kmol
Ru    = 8314.47;    % J/(kmol·K)

VD  = 1.4e-3;       % m^3 (total)
vd  = VD/4;         % m^3 per cylinder
b   = (VD/pi)^(1/3);
s   = b;
Rg  = (2*0.18)/s;   % rod-to-crank ratio
vc  = vd/(rc-1);
vl  = vc*rc;

%%  __________ Iterative convergence loop 

iter = 0;
stuck = '';
while (r1>ae || r2>ae || r3>ae) && iter < maxIter
    iter = iter + 1;

    [W, imep2, Qin, eta100, errs, traces] = run_cycle( ...
        a,g,bmep,nWie,rc,Mair,Ru,vd,Rg,vl, ...
        F,qc,as,Mf,ts,td,Te,Pi,Pe,f1);

    r1 = errs.r1; r2 = errs.r2; r3 = errs.r3;

    % Signed, damped updates toward targets 
    kp = 1e-3;   % bar step gain for IMEP via Pi
    kf = 1e-3;   % fraction gain toward f2
    kt = 1e-3;   % temperature gain toward T5

    % Move Pi to make IMEP match bmep*1e5 (proxy); clamp to [0.1, 1.5] bar

    Pi = Pi + kp * (bmep*1e5 - errs.imep2) / 1e5;    % scale to bar
    Pi = min(max(Pi, 0.1), 1.5);

    % Move f1 toward its target f2; clamp
    f1 = f1 + kf * (errs.f2 - f1);
    f1 = min(max(f1, 0.0), 0.2);

    % Move Te toward target T5; clamp (raise upper bound so we don’t saturate)
    Te = Te + kt * (errs.T5 - Te);
    Te = min(max(Te, 500), 2000);

    if iter==maxIter
        [~, idxMax] = max([r1 r2 r3]);
        labels = {'r1','r2','r3'};
        stuck = labels{idxMax};
    end

end

%  __________  Results 
theta = traces.theta; P = traces.P; V1 = traces.V1; T = traces.T;

if iter >= maxIter && (r1>ae || r2>ae || r3>ae)
    warning('Did not converge: maxIter reached. Stuck on %s = %.4g (tol %.4g).', stuck, max([r1 r2 r3]), ae);
end

disp('--- Final Results ---');
fprintf('Iterations                    = %d\n', iter);
fprintf('Indicated work per cycle [J]  = %.3f\n', traces.W);
fprintf('Indicated efficiency   [%%]    = %.2f\n', eta100);
fprintf('Converged Pi [bar]            = %.4f\n', Pi);
fprintf('Converged Pe [bar]            = %.4f\n', Pe);
fprintf('Converged Te [K]              = %.1f\n', Te);
fprintf('Converged f1 [-]              = %.5f\n', f1);
fprintf('Errors: r1=%.3g, r2=%.3g, r3=%.3g (tol %.3g)\n', r1, r2, r3, ae);

%%  __________  Plots (pressure full cycle; temp classic slice; PV curve)

figure('Name','Engine Cycle');
subplot(2,2,1)
plot(theta, P/1e5, 'r'); grid on;
title('Gasoline');
xlabel('crank angle [deg]');
ylabel('Pressure [bar]');

subplot(2,2,2)
plot(theta, T, 'r'); grid on;   % classic view
xlabel('crank angle [deg]');
ylabel('Temperature [K]');

subplot(2,2,3:4)
plot(V1, P/1e5, 'r'); grid on;
xlabel('Volume [m^3]');
ylabel('Pressure [bar]');

%  __________  Inner function 
function [W, imep2, Qin, eta100, errs, tr] = run_cycle( ...
    a,g,bmep,nWie,rc,Mair,Ru,vd,Rg,vl, ...
    F,qc,as,Mf,ts,td,Te,Pi,Pe,f1)

% Mixture properties
ys   = 1/(1+4.76*as);
yair = 1 - ys;
M    = yair*Mair + ys*Mf;
R    = Ru/M;
Ti   = 35 + 273.15;  % K

% Energy input
qin  = F/(1+F)*qc;               % kJ/kg mix
T1   = (1-f1)*Ti + f1*Te*(1-(1-Pi/Pe)*((g-1)/g));
m1   = (Pi*1e5*vl)/(R*T1);
Qin  = m1*(1-f1)*qin*1e3;        % J

% Crank-angle grid
theta = -360:1:360;
nSteps = numel(theta);

V1 = zeros(1,nSteps);
V2 = zeros(1,nSteps);
xb = zeros(1,nSteps);
DQ = zeros(1,nSteps);
P  = zeros(1,nSteps);
dP = zeros(1,nSteps);

% Cycle integration
k = 0;
for iDeg = theta
    k = k+1;
    o = deg2rad(iDeg);

    % Kinematics
    V1(k) = vd/(rc-1) + (vd/2)*(Rg+1 - cos(o) - sqrt(Rg^2 - (sin(o))^2));
    V2(k) = (vd/2)*sin(o)*(1 + cos(o)/sqrt(Rg^2-(sin(o))^2));

    % Wiebe burn fraction with td
    if iDeg < ts
        xb(k) = 0;
    elseif iDeg > (ts+td)
        xb(k) = 1;
    else
        xb(k) = 1 - exp(-a*((iDeg-ts)/td)^nWie);
    end

    % Heat release rate with td
    if iDeg >= ts && iDeg <= ts+td
        DQ(k) = nWie*a*Qin/deg2rad(td) * (1-xb(k)) * ...
            ((deg2rad(iDeg-ts)/deg2rad(td))^(nWie-1));
    else
        DQ(k) = 0;
    end

    % Pressure update
    if iDeg < -180
        P(k) = Pi*1e5;
    elseif iDeg > 180
        P(k) = Pe*1e5;
    elseif k > 1
        P(k) = P(k-1) + dP(k-1)*deg2rad(1);
    end

    dP(k) = (-g*P(k)/V1(k))*V2(k) + ((g-1)/V1(k))*DQ(k);
end

% Work via PV integral (trapz across compression+expansion window)
mask = (theta>=-180 & theta<=180);
W = trapz(V1(mask), P(mask));   % J
imep2 = W / vd;
eta100 = (W/Qin)*100;

% Temperatures (full cycle for physics; plotting can slice later)
T = P .* V1 / (m1*R);

% Errors (original logic)
P4 = P(theta==180);
T4 = T(theta==180);
T5 = T4*(Pe*1e5/P4)^((g-1)/g);
f2 = (1/rc)*((Pe*1e5/P4))^(1/g);

r1 = abs(1 - (imep2/(bmep*1e5)));  % IMEP match
r2 = abs(1 - (f1/f2));             % burned fraction match proxy
r3 = abs(1 - T5/Te);               % exhaust T match

errs = struct('r1',r1,'r2',r2,'r3',r3);

% add targets so the outer loop can move toward them:
errs.imep2 = imep2;   % actual IMEP
errs.f2    = f2;      % target f1
errs.T5    = T5;      % target Te

tr = struct('theta',theta,'P',P,'V1',V1,'T',T,'W',W);
end
