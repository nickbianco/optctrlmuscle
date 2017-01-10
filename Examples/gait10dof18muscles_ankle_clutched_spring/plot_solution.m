import org.opensim.modeling.*

numDOFs = DatStore.nDOF;
numMuscles = DatStore.nMuscles;

time = OptInfo.result.solution.phase.time;

% Extract experimental data.
expTime = DatStore.time;
qExp = DatStore.q_exp;
momArmsExp = DatStore.dM;
momArms = interp1(expTime, momArmsExp, time);

%for idof = 1:numDOFs
%    qSpline(idof) = spline(expTime, pi / 180. * qExp(:, idof));
%end

% Extract parts of the solution related to the device.
control = OptInfo.result.solution.phase.control;
state = OptInfo.result.solution.phase.state;

% Get controls
e       = control(:,1:numMuscles);
aT      = control(:,numMuscles+1:numMuscles+numDOFs);
vMtilde = control(:,numMuscles+numDOFs+1:end);

% Get states
a       = state(:,1:numMuscles);
lMtilde = state(:,numMuscles+1:end);

% Joint moment breakdown.
deviceIndices = strmatch('ankle_angle', DatStore.DOFNames);

Topt = 2*150.0;
for idof = 1:numDOFs
    %subplot(numDOFs, 1, idof);
    figure;
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
%     deviceIndex = find(deviceIndices == idof);
%     if ~isempty(deviceIndex)
%         deviceMoment = TODO;
%         plot(time, deviceMoment);
%         legendEntries = [legendEntries {'device'}];
%         sumMoment = sumMoment + deviceMoment;
%     end
    plot(time(1:end-1), sumMoment(1:end-1), 'r', 'LineWidth', 2);
    legendEntries = [legendEntries {'sum'}];
    legend(legendEntries, 'Interpreter', 'none');
    title(DatStore.DOFNames{idof}, 'Interpreter', 'none');
end

% TODO percent reduction in metabolic cost / sum squared activation and
% squared excitation.
model = Model('subject05.osim');
muscles = model.getMuscles();
muscles.getName()
keyboard

for m = 1:numMuscles
    
    for i = 1:length(time)
        
        
    end
    
end


rho = 1059.7;           % Muscle density [kg/m^3]
Fmax = 3127;            % Max isometric force [N]
sigma = 600000;         % Specific tension [N/m^2]
Lceopt = 0.055;         % Optimal fiber length [m]
PCSA = Fmax/sigma;      % Physiological cross sectional area [m^2]
mass = PCSA*rho*Lceopt; % Muscle mass [kg]

maxFiberVel = 12;       % Fiber-lengths per second
param_rFT = 0.2;        % Proportion of fast-twitch muscle fibers
param_Arel = 0.1 + 0.4*param_rFT;
param_Brel = param_Arel*maxFiberVel;
params = struct('Lceopt',Lceopt, 'Arel',param_Arel, ...
                'Brel',param_Brel, 'Fmax',Fmax, 'rFT',param_rFT, ...
                'VceMax_LceoptsPerSecond',param_Brel/param_Arel, ...
                'muscleMass',mass, 'scalingFactorS',1.0, ...
                'versionNumber',2010);
heatRates = calcUmbergerProbe(Lce,Vce,F,Fiso,u,a,params)


% TODO left and right limbs together.





















