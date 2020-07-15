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
        fprintf('\nWhole Body calculation:\n')
        if any(any(e < -0.1))
            e(e < -0.1)
            warning('VERY negative excitation! ')
            if length(e(e < -0.1)) < 5
                fprintf('But very few time points...still clipping at 0.')
                culprits = e(e<0)
                e(e < 0) = 0;
            else
                fprint('Exiting...')
                avg_total_rate = NaN;
                return
            end
        else
            warning('Slightly negative excitation...clipping at 0.');
            culprits = e(e<0)
            e(e < 0) = 0;
        end
    end
    if any(any(e > 1))
        fprintf('\nWhole Body calculation:\n')
        if any(any(e > 1.1))
            e(e > 1.1)
            warning('Excitation much greater than 1!');
            if length(e(e > 1.1)) < 5
                fprintf('But very few time points...still clipping at 1.')
                culprits = e(e>1)
                e(e > 1) = 1;
            else
                fprint('Exiting...')
                avg_total_rate = NaN;
                return
            end
        else
            warning('Excitation slightly greater than 1...clipping at 1.');
            culprits = e(e>1)
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
        
    %     sigma = probeUmberger.getSpecificTension(MuscleNames{m}); % Specific tension [N/m^2]
        sigma = 300000;
        PCSA = Fmax/sigma;      % Physiological cross sectional area [m^2]
        mass = PCSA*rho*Lceopt; % Muscle mass [kg]
        
        paramsUmb = struct('Lceopt',Lceopt, 'rFT',param_rFT, ...
                    'VceMax_LceoptsPerSecond',maxFiberVel, ...
                    'muscleMass',mass, 'scalingFactorS',1.5, ... % 1.5: aerobic.
                    'versionNumber',2003); % 2010
        VCEmax_mps = paramsUmb.VceMax_LceoptsPerSecond * Lceopt; % [m/s]
        
        heatRates = NaN(numColPoints,5);
        for i = 1:numColPoints
            % Lce = lMtilde(i,m)*Lceopt;
            % Vce = vMtilde(i,m)*VCEmax_mps;
            % heatRates(i,:) = calcUmbergerProbe(Lce,Vce,F(i,m),Fiso(i,m),e(i,m),a(i,m),paramsUmb);
            % going to try the tom recruitment thing
            u_slow = sin((pi/2).*e(i,m));
            u_fast = 1 - cos((pi/2).*e(i,m));
            
            if e(i,m) == 0
                f_rec_slow = 1;
            else
                f_rec_slow = (rST*u_slow)/((rST*u_slow) + ((1-rST)*u_fast));
            end
            
            param_rFT_new = 1 - f_rec_slow;
            paramsUmb.rFT = param_rFT_new;
            
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


%{
  Need
  - Lce:            Length of the contractile element, this is the length of the 
                    muscle at each given timestep
  - Vce             velocity of the muscle fiber also at each timestep
  - F               force in the muscle tendon unit
  - Fiso            the maximum isometric force the muscle is capable of
  - u               excitation in the muscle
  - a               activation in the muscle
  - params          optimal muscle length, ratio of fast and slow twitch 
                    (tweaked based on excitation), max fiber velocity,
                    muscle mass, scaling factor for aerobic or anerobic
                    version of the umberger equation

Need to find all of these values for each timestep and for each muscle
then:

loop through muscles
    loop through time
        calcUmbergerProbe(Lce, Vce, F, Fiso, u, a, params);
    end
end

%} 

% get time vector and the number of time points evaluated
Time = solution.getTimeMat();
numColPoints = solution.getNumTimes();


% get excitations
controlData = solution.getControlsTrajectoryMat();
controlNames_os = solution.getControlNames();
controlNames = [];
for i = 0:controlNames_os.size()-1
    controlNames= [controlNames, controlNames_os.get(i)];
end

% get activations
stateData = solution.getStatesTrajectoryMat();
stateNames_os = solution.getStateNames();
stateNames = [];
for i=0:stateNames_os.size()-1
    stateNames = [stateNames, stateNames_os.get(i)];
end

% now need to get the muscle length, velocity, and force in the muscle... etc. 




