function [avg_total_rate] = calcWholeBodyMetabolicRate(model, mat)
import org.opensim.modeling.*

Time = mat.Time;
numColPoints = length(Time);

DatStore = mat.DatStore;
OptInfo = mat.OptInfo;
MuscleNames = mat.MuscleNames;
auxdata = mat.OptInfo.result.setup.auxdata;

numDOFs = DatStore.nDOF;
numMuscles = DatStore.nMuscles;

control = OptInfo.result.solution.phase.control;
state = OptInfo.result.solution.phase.state;

% Get controls
if strcmp(DatStore.formulation,'Ftilde')
    vA = 100*control(:,1:numMuscles);
    aT = control(:,numMuscles+1:numMuscles+numDOFs);
    dFtilde = 10*control(:,numMuscles+numDOFs+1:numMuscles+numDOFs+numMuscles);
   
else
    e = control(:,1:numMuscles);
    aT = control(:,numMuscles+1:numMuscles+numDOFs);
    vMtilde = control(:,numMuscles+numDOFs+1:numMuscles+numDOFs+numMuscles);
end

% Get states
if strcmp(DatStore.formulation,'Ftilde')
    a = state(:,1:numMuscles);
    e = computeExcitationRaasch(a, vA, auxdata.tauDeact, auxdata.tauAct);

    Ftilde = state(:,numMuscles+1:numMuscles+numMuscles);
    
    lMTinterp = interp1(DatStore.time,DatStore.LMT,Time);
    for m = 1:numMuscles
        LMTSpline(m) = spline(Time,lMTinterp(:,m));
        [LMT(:,m),VMT(:,m),~] = SplineEval_ppuval(LMTSpline(m),Time,1);
    end
        
   [lM,lMtilde,vM,vMtilde] = FiberLengthVelocity_Ftilde(Ftilde,dFtilde, auxdata.params, LMT, VMT, auxdata.Fpparam);
else
    a = state(:,1:numMuscles);
    lMtilde = state(:,numMuscles+1:numMuscles+numMuscles);
end

% Check excitation
if any(any(e < 0))
    if any(any(e < -0.1))
        e(e < -0.1)
        warning('VERY negative excitation! ')
        if length(e(e < -0.1)) < 5
            fprintf('But very few time points...still clipping at 0.')
            e(e < 0) = 0;
        else
            fprint('Exiting...')
            avg_total_rate = NaN;
            return
        end
    else
        warning('Slightly negative excitation...clipping at 0.');
        e(e < 0) = 0;
    end
end
if any(any(e > 1))
    if any(any(e > 1.1))
        e(e > 1.1)
        warning('Excitation much greater than 1!');
        if length(e(e > 1.1)) < 5
            fprintf('But very few time points...still clipping at 1.')
            e(e > 1) = 1;
        else
            fprint('Exiting...')
            avg_total_rate = NaN;
            return
        end
    else
        warning('Excitation slightly greater than 1...clipping at 1.');
        e(e > 1) = 1;
    end
end

% Metabolic cost
musclesApoorva = model.getMuscles();

probeSet = model.getProbeSet();
probe = probeSet.get('metabolic_power');
probeUmberger = Umberger2010MuscleMetabolicsProbe.safeDownCast(probe);

rho = 1059.7; % Muscle density [kg/m^3]

lMT = NaN(numColPoints, numMuscles);
for m = 1:numMuscles
    lMT(:, m) = ppval(auxdata.LMTSpline(m), Time);
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
    
    musc_energy_rate(:,m) = heatRates(:,5) * mass;
    
end

state = model.initSystem();
bodyMass = model.getTotalMass(state);
wholebody_energy_rate = sum(musc_energy_rate,2);
duration = Time(end) - Time(1);
norm_average_wholebody_energy_rate = trapz(mat.Time, wholebody_energy_rate) / bodyMass / duration;
avg_total_rate = norm_average_wholebody_energy_rate;
end
