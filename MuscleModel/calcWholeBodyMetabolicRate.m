function [UmbergerKoelewijn2018, MinettiAlexander1997] = calcWholeBodyMetabolicRate(model, mat)
import org.opensim.modeling.*

Time = mat.Time;
numColPoints = length(Time);

DatStore = mat.DatStore;
OptInfo = mat.OptInfo;
MuscleNames = mat.MuscleNames;
auxdata = mat.OptInfo.result.setup.auxdata;
walking_speed = auxdata.speed;

numDOFs = DatStore.nDOF;
numMuscles = DatStore.nMuscles;

control = OptInfo.result.solution.phase.control;
state = OptInfo.result.solution.phase.state;

% Get controls
if strcmp(DatStore.formulation,'Ftilde')
    if strcmp(auxdata.actdyn, 'implicit')
        vA = 100*control(:,1:numMuscles);
    elseif strcmp(auxdata.actdyn, 'explicit')
        e = control(:,1:numMuscles);
    end
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
    if strcmp(auxdata.actdyn, 'implicit')
        e = computeExcitationRaasch(a, vA, auxdata.tauDeact, auxdata.tauAct);
    end

    Ftilde = state(:,numMuscles+1:numMuscles+numMuscles);
    
    lMTinterp = interp1(DatStore.time,DatStore.LMT,Time);
    for m = 1:numMuscles
        LMTSpline(m) = spline(Time,lMTinterp(:,m));
        [LMT(:,m),VMT(:,m),~] = SplineEval_ppuval(LMTSpline(m),Time,1);
    end
        
   [lM,lMtilde,vM,vMtilde] = FiberLengthVelocity_Ftilde(Ftilde, dFtilde, ...
       auxdata.params, LMT, VMT, auxdata.Fpparam);
else
    a = state(:,1:numMuscles);
    lMtilde = state(:,numMuscles+1:numMuscles+numMuscles);
end

% Check excitation. We'll tolerate 
if any(any(e < 0))
    if any(any(e < -0.1))
        e(e < -0.1)
        warning('VERY negative excitation! ')
        if length(e(e < -0.1)) < 5
            fprintf('But very few time points...still clipping at 0.')
            e(e < 0) = 0;
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
muscleData = DeGroote2016Muscle_lMtildeState(a,lMtilde,vMtilde, ...
    lMT, auxdata.params, auxdata.Fvparam, auxdata.Fpparam, ...
    auxdata.Faparam);

musc_metabolic_rate_UmbKoel = NaN(numColPoints,numMuscles);
musc_metabolic_rate_MinAlex = NaN(numColPoints,numMuscles);
for m = 1:numMuscles
    
    musc = musclesApoorva.get(MuscleNames{m});
    Fmax = musc.getMaxIsometricForce; % Max isometric force [N]
    Lceopt = musc.getOptimalFiberLength; % Optimal fiber length [m]
    maxFiberVel = musc.getMaxContractionVelocity();
    
    rST = probeUmberger.getRatioSlowTwitchFibers(MuscleNames{m});
    rFT = 1 - rST; % Proportion of fast-twitch muscle fibers
    
    % Specific tension [N/m^2]
    sigma = probeUmberger.getSpecificTension(MuscleNames{m}); 
    PCSA = Fmax/sigma;      % Physiological cross sectional area [m^2]
    mass = PCSA*rho*Lceopt; % Muscle mass [kg]
    
    S = 1.5; % 1.5: aerobic    
    metabolicParams = [rFT, Lceopt, maxFiberVel, mass, S];
    VCEmax_mps = maxFiberVel * Lceopt; % [m/s]
    
    % Store UmbergerKoelewijn2018 muscle energy rate. 
    % calcUmbergerKoelewijn2018Cost() returns individual muscle rates in W/kg 
    % where the mass normalization is muscle mass. Multiply by muscle mass to 
    % return muscle energy rates in W.
    Lce = lMtilde(:,m)*Lceopt;
    Vce = vMtilde(:,m)*VCEmax_mps;
    musc_metabolic_rate_UmbKoel(:,m) = mass * calcUmbergerKoelewijn2018Cost(...
            Lce, Vce, muscleData.Fce(:,m), muscleData.FMltilde(:,m), e(:,m), ...  
            a(:,m), metabolicParams);
        
    % Store MinettiAlexander1997 muscle energy rate.
    musc_metabolic_rate_MinAlex(:,m) = calcMinettiAlexanderProbe(...
            vMtilde(:,m)*VCEmax_mps, VCEmax_mps, Fmax, a(:,m));
end

% Compute and store metabolic info.
state = model.initSystem();
bodyMass = model.getTotalMass(state);
duration = Time(end) - Time(1);

wholebody_metabolic_rate_UmbKoel = sum(musc_metabolic_rate_UmbKoel, 2); % [W]
UmbergerKoelewijn2018.average_wholebody_metabolic_rate = ... 
    trapz(mat.Time, wholebody_metabolic_rate_UmbKoel) / bodyMass / duration; % [W/kg]
UmbergerKoelewijn2018.wholebody_cost_of_transport = ...
    UmbergerKoelewijn2018.average_wholebody_metabolic_rate / walking_speed; % [W/kg/(m/s)]
UmbergerKoelewijn2018.muscle_average_metabolic_rates = ...
    trapz(mat.Time, musc_metabolic_rate_UmbKoel) / bodyMass / duration; % [W/kg]

wholebody_metabolic_rate_MinAlex = sum(musc_metabolic_rate_MinAlex, 2); % [W]
MinettiAlexander1997.average_wholebody_metabolic_rate = ...
    trapz(mat.Time, wholebody_metabolic_rate_MinAlex) / bodyMass / duration; % [W/kg]
MinettiAlexander1997.cost_of_transport = ...
    MinettiAlexander1997.average_wholebody_metabolic_rate / walking_speed; % [W/kg/(m/s)]
MinettiAlexander1997.muscle_average_metabolic_rates = ...
    trapz(mat.Time, musc_metabolic_rate_MinAlex) / bodyMass / duration; % [W/kg]

end
