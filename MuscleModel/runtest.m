% Compare calcMuscleForces() and calcUmbergerProbe() with Umberger's results
% (see Figure 1A-D in Umberger et al., 2003).
%
% Uchida.1611021431

clear all; close all; clc; format long;

% Set all parameters used by calcMuscleForces() and calcUmbergerProbe().
maxFiberVel = 12; %fiber-lengths per second
param_rFT = 0.2; %proportion of fast-twitch muscle fibers
param_Arel = 0.1 + 0.4*param_rFT;
param_Brel = param_Arel*maxFiberVel;
params = struct('width',0.80, 'Lceopt',0.055, 'Arel',param_Arel, ...
                'Brel',param_Brel, 'Fmax',3127, 'rFT',param_rFT, ...
                'VceMax_LceoptsPerSecond',param_Brel/param_Arel, ...
                'muscleMass',0.805, 'scalingFactorS',1.0, ...
                'versionNumber',2003);
metabolicParams = [param_rFT, 0.055, param_Brel/param_Arel, 0.805, 1.0];
VCEmax_mps = params.VceMax_LceoptsPerSecond * params.Lceopt; %[m/s]

% Position and velocity of muscle endpoint.
NPTS = 401;
p = params.Lceopt*ones(NPTS,1);
vnorm = linspace(-1,1,NPTS)';
v = vnorm*VCEmax_mps;

% Excitation and activation.
u = ones(NPTS,1);
a = u;

% Calculate muscle forces F(t) and Fiso(t).
idx_F = 1;
idx_Fiso = 2;
muscleForces = NaN(NPTS,2);
for i = 1:NPTS
    muscleForces(i,:) = calcMuscleForces(a(i), p(i), v(i), params);
end

%% original probe
% Calculate heat rates.
idx_UmbA = 1;
idx_UmbM = 2;
idx_UmbSL = 3;
idx_UmbW = 4;
idx_UmbTotal = 5;
heatRates = NaN(NPTS,5);
for i = 1:NPTS
	heatRates(i,:) = calcUmbergerProbe(p(i), v(i), muscleForces(i,idx_F), ...
                     muscleForces(i,idx_Fiso), u(i), a(i), params);
end

% Shortening plots (match Figure 1A-D in Umberger et al., 2003).
idx0 = 1;
idx1 = (NPTS+1)/2;

% Normalized velocity from -1 to 0 --> shortening velocity from 1 to 0.
shorteningVel = -v(idx0:idx1) / VCEmax_mps;
% Normalized CE force [N/N].
force = muscleForces(idx0:idx1,idx_F) / params.Fmax;
% Mechanical power [W/kg].
mech_power = heatRates(idx0:idx1,idx_UmbW);
% Total rate of energy liberation [W/kg].
total_energy_rate = heatRates(idx0:idx1,idx_UmbTotal);
% Mechanical efficiency (in [0,1]).
efficiency = mech_power ./ total_energy_rate;

figure(1);
hold on
subplot(4,1,1); plot(shorteningVel,force); ylabel('Force (F/Fmax)');
subplot(4,1,2); plot(shorteningVel,mech_power); ylabel('Power (W/kg)');
subplot(4,1,3); plot(shorteningVel,total_energy_rate); ylabel('Energy Rate (W/kg)');
subplot(4,1,4); plot(shorteningVel,efficiency); ylabel('Efficiency');
xlabel('Shortening Speed (V/Vmax)');

%% smooth probe
% Calculate heat rates.
idx_UmbA = 1;
idx_UmbW = 2;
idx_UmbTotal = 3;
heatRates = NaN(NPTS,3);
for i = 1:NPTS
	heatRates(i,:) = calcUmbergerCost2010(p(i), v(i), muscleForces(i,idx_F), ...
                         muscleForces(i,idx_Fiso), u(i), a(i), metabolicParams);
end

% Shortening plots (match Figure 1A-D in Umberger et al., 2003).
idx0 = 1;
idx1 = (NPTS+1)/2;

% Normalized velocity from -1 to 0 --> shortening velocity from 1 to 0.
shorteningVel = -v(idx0:idx1) / VCEmax_mps;
% Normalized CE force [N/N].
force = muscleForces(idx0:idx1,idx_F) / params.Fmax;
% Mechanical power [W/kg].
mech_power = heatRates(idx0:idx1,idx_UmbW);
% Total rate of energy liberation [W/kg].
total_energy_rate = heatRates(idx0:idx1,idx_UmbTotal);
% Mechanical efficiency (in [0,1]).
efficiency = mech_power ./ total_energy_rate;

figure(2);
hold on
subplot(4,1,1); plot(shorteningVel,force); ylabel('Force (F/Fmax)');
subplot(4,1,2); plot(shorteningVel,mech_power); ylabel('Power (W/kg)');
subplot(4,1,3); plot(shorteningVel,total_energy_rate); ylabel('Energy Rate (W/kg)');
subplot(4,1,4); plot(shorteningVel,efficiency); ylabel('Efficiency');
xlabel('Shortening Speed (V/Vmax)');
