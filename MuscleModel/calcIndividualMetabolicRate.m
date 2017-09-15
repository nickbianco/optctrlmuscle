function [norm_average_muscle_energy_rate] = calcIndividualMetabolicRate(model, mat)
import org.opensim.modeling.*

time = mat.Time;
numColPoints = length(time);

DatStore = mat.DatStore;
OptInfo = mat.OptInfo;
MuscleNames = mat.MuscleNames;
auxdata = mat.OptInfo.result.setup.auxdata;

numDOFs = DatStore.nDOF;
numMuscles = DatStore.nMuscles;

control = OptInfo.result.solution.phase.control;
state = OptInfo.result.solution.phase.state;

% Get controls
e       = control(:,1:numMuscles); 
%    e(e<0)=0;
if any(any(e < 0))
    if any(any(e < -0.1))
        e(e < -0.1)
        error('VERY negative excitation!');
    else
        warning('Slightly negative excitation...clipping at 0.');
        e(e < 0) = 0;
    end
end
if any(any(e > 1))
    if any(any(e > 1.1))
        e(e > 1.1)
        error('Excitation much greater than 1!');
    else
        warning('Excitation slightly greater than 1...clipping at 1.');
        e(e > 1) = 1;
    end
end
aT      = control(:,numMuscles+1:numMuscles+numDOFs);
vMtilde = control(:,numMuscles+numDOFs+1:end);

% Get states
a       = state(:,1:numMuscles);
lMtilde = state(:,numMuscles+1:end);

% Metabolic cost
musclesApoorva = model.getMuscles();

probeSet = model.getProbeSet();
probe = probeSet.get('metabolic_power');
probeUmberger = Umberger2010MuscleMetabolicsProbe.safeDownCast(probe);

rho = 1059.7; % Muscle density [kg/m^3]

lMT = NaN(numColPoints, numMuscles);
for m = 1:numMuscles
    lMT(:, m) = ppval(auxdata.LMTSpline(m), time);
end
[~, ~, F, Fiso] = DeGroote2016Muscle_lMtildeState(a,lMtilde,vMtilde, ...
    lMT, auxdata.params, auxdata.Fvparam, auxdata.Fpparam, ...
    auxdata.Faparam);

musc_energy_rate = NaN(numColPoints,numMuscles);
for m = 1:numMuscles
    
    musc = musclesApoorva.get(MuscleNames{m});
    Fmax = musc.getMaxIsometricForce;   % Max isometric force [N]
    Lceopt = musc.getOptimalFiberLength;         % Optimal fiber length [m]
    maxFiberVel = musc.getMaxContractionVelocity();
    rST = probeUmberger.getRatioSlowTwitchFibers(MuscleNames{m});
    param_rFT = 1 - rST;        % Proportion of fast-twitch muscle fibers
    
    sigma = probeUmberger.getSpecificTension(MuscleNames{m}); % Specific tension [N/m^2]
    PCSA = Fmax/sigma;      % Physiological cross sectional area [m^2]
    mass = PCSA*rho*Lceopt; % Muscle mass [kg]
    
    paramsUmb = struct('Lceopt',Lceopt, 'rFT',param_rFT, ...
                'VceMax_LceoptsPerSecond',maxFiberVel, ...
                'muscleMass',mass, 'scalingFactorS',1.5, ... % 1.5: aerobic.
                'versionNumber',2010);
    VCEmax_mps = paramsUmb.VceMax_LceoptsPerSecond * Lceopt; % [m/s]
    
    heatRates = NaN(numColPoints,5);
    for i = 1:numColPoints
        Lce = lMtilde(i,m)*Lceopt;
        Vce = vMtilde(i,m)*VCEmax_mps;
        heatRates(i,:) = calcUmbergerProbe(Lce,Vce,F(i,m),Fiso(i,m),e(i,m),a(i,m),paramsUmb);
    end
    musc_energy_rates(:,:,m) = heatRates(:,:) * mass;
      
end

% heatRates = [Activation, Maintenance, Shortening/Lengthening (shortening
% is negative), Mechanical work rate (positive when shortening), sum]

state = model.initSystem();
bodyMass = model.getTotalMass(state);
duration = time(end) - time(1);

for i=1:5
    for m=1:numMuscles
        norm_average_muscle_energy_rate(i,m) = trapz(mat.Time,musc_energy_rates(:,i,m))/bodyMass/duration;
    end
end

end
