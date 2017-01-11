import org.opensim.modeling.*

numDOFs = DatStore.nDOF;
numMuscles = DatStore.nMuscles;

time = OptInfo.result.solution.phase.time;
NPTS = length(time);

% Extract experimental data.
expTime = DatStore.time;
qExp = DatStore.q_exp;
momArmsExp = DatStore.dM;
momArms = interp1(expTime, momArmsExp, time);
jointAngles = pi / 180. * interp1(expTime, qExp, time);

% Extract parts of the solution related to the device.
control = OptInfo.result.solution.phase.control;
state = OptInfo.result.solution.phase.state;

% Get controls
e       = control(:,1:numMuscles); e(e<0)=0;
aT      = control(:,numMuscles+1:numMuscles+numDOFs);
vMtilde = control(:,numMuscles+numDOFs+1:end);

% Get states
a       = state(:,1:numMuscles);
lMtilde = state(:,numMuscles+1:end);

% Joint moment breakdown.
deviceIndices = strmatch('ankle_angle', DatStore.DOFNames);
assert(length(deviceIndices) == 1);

for idof = 1:numDOFs
    subplot(numDOFs, 1, idof);
    hold on;
    plot(expTime, DatStore.T_exp(:, idof), 'k', 'LineWidth', 2);
    legendEntries = {'net'};
    sumMoment = zeros(length(TForce(:, 1)), 1);
    for imusc = 1:numMuscles
        if any(momArms(:, idof, imusc)) > 0.00001
            thisMoment = TForce(:, imusc) .* momArms(:, idof, imusc);
            plot(time(1:end-1), thisMoment(1:end-1));
            legendEntries = [legendEntries {MuscleNames{imusc}}];
            sumMoment = sumMoment + thisMoment;
        end
    end
     deviceIndex = find(deviceIndices == idof);
     if ~isempty(deviceIndex)
         normSpringStiff = OptInfo.result.solution.parameter(1);
         maxSpringStiff = 400; % N-m/rad.
         rest_angle = OptInfo.result.solution.parameter(2);
         ankleAngle = -(jointAngles(:, idof) - rest_angle);
         deviceMoment = maxSpringStiff * normSpringStiff .* ankleAngle;
         plot(time, deviceMoment);
         legendEntries = [legendEntries {'device'}];
         sumMoment = sumMoment + deviceMoment;
     end
    plot(time(1:end-1), sumMoment(1:end-1), 'r', 'LineWidth', 2);
    legendEntries = [legendEntries {'sum'}];
    legend(legendEntries, 'Interpreter', 'none');
    title(DatStore.DOFNames{idof}, 'Interpreter', 'none');
    if idof == numDOFs
        xlabel('time (s)');
    end
    ylabel('moment (N-m)');
end

% Metabolic cost
modelApoorva = Model('subject05.osim');
musclesApoorva = modelApoorva.getMuscles();
MuscleNamesApoorva = {'bifemlh_r','bifemsh_r','glut_max1_r', ...
                      'psoas_r','rect_fem_r','vas_med_r', ...
                      'med_gas_r','soleus_r','tib_ant_r'};                  
muscleMap = containers.Map(MuscleNames,MuscleNamesApoorva);

probeSet = modelApoorva.getProbeSet();
probe = probeSet.get('metabolic_power');
probeUmberger = Umberger2010MuscleMetabolicsProbe.safeDownCast(probe);

rho = 1059.7; % Muscle density [kg/m^3]
maxFiberVel = 12;  % Fiber-lengths per second

for m = 1:numMuscles
    m
    
    muscleNameApoorva = muscleMap(MuscleNames{m});
    musc = musclesApoorva.get(muscleNameApoorva);
    Fmax = musc.getMaxIsometricForce;   % Max isometric force [N]
    Lceopt = 0.055;         % Optimal fiber length [m]
    
    rST = probeUmberger.getRatioSlowTwitchFibers(muscleNameApoorva);
    param_rFT = 1 - rST;        % Proportion of fast-twitch muscle fibers
    param_Arel = 0.1 + 0.4*param_rFT;
    param_Brel = param_Arel*maxFiberVel;
    
    sigma = probeUmberger.getSpecificTension(muscleNameApoorva); % Specific tension [N/m^2]
    PCSA = Fmax/sigma;      % Physiological cross sectional area [m^2]
    mass = PCSA*rho*Lceopt; % Muscle mass [kg]
    
    params = struct('width',0.80,'Lceopt',Lceopt, 'Arel',param_Arel, ...
                'Brel',param_Brel, 'Fmax',Fmax, 'rFT',param_rFT, ...
                'VceMax_LceoptsPerSecond',param_Brel/param_Arel, ...
                'muscleMass',mass, 'scalingFactorS',1.0, ...
                'versionNumber',2010);
    VCEmax_mps = params.VceMax_LceoptsPerSecond * Lceopt; % [m/s]
    
    heatRates = NaN(NPTS,5);
    for i = 1:NPTS
        Lce = lMtilde(i,m)*Lceopt;
        Vce = vMtilde(i,m)*VCEmax_mps;
        forces = calcMuscleForces(a(i,m),Lce,Vce,params);
        F = forces(1);
        Fiso = forces(2);
        heatRates(i,:) = calcUmbergerProbe(Lce,Vce,F,Fiso,e(i,m),a(i,m),params);
    end
    
    total_energy_rate(:,m) = heatRates(:,5) * mass;
    
end


% TODO left and right limbs together.





















