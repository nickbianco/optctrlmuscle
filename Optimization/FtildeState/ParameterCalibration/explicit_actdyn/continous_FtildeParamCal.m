function phaseout = continous_FtildeParamCal(input)

% Get input data
NMuscles        = input.auxdata.NMuscles;
Ndof            = input.auxdata.Ndof;
tauAct          = input.auxdata.tauAct;
tauDeact        = input.auxdata.tauDeact;
params          = input.auxdata.params;
splinestruct    = input.auxdata.splinestruct;
numColPoints    = size(input.phase.state,1);

% Get controls
e        = input.phase.control(:,1:NMuscles);
aT       = input.phase.control(:,NMuscles+1:NMuscles+Ndof);
dFtilde  = 10*input.phase.control(:,NMuscles+Ndof+1:end);

% Get states
a       = input.phase.state(:,1:NMuscles);
Ftilde  = input.phase.state(:,NMuscles+1:end);

% Get parameters
parameters = input.phase.parameter;

% Modify muscle-tendon properties based on parameters
terms = input.auxdata.parameterCalibrationTerms;
paramIndices = input.auxdata.parameterCalibrationIndices;
musclesToCalibrate = fieldnames(paramIndices);
MuscleNames = input.auxdata.MuscleNames;

for m = 1:length(musclesToCalibrate)
    muscIdx = find(contains(MuscleNames, musclesToCalibrate{m}));
    paramsToCalibrate = fieldnames(paramIndices.(musclesToCalibrate{m}));
    for p = 1:length(paramsToCalibrate)
        idx = paramIndices.(musclesToCalibrate{m}).(paramsToCalibrate{p});
        paramVal = parameters(1,idx);
        switch paramsToCalibrate{p}
            case 'optimal_fiber_length'
                params(9,muscIdx) = paramVal;
            case 'tendon_slack_length'
                params(10,muscIdx) = paramVal;
            case 'pennation_angle'
                params(11,muscIdx) = paramVal;
            case 'muscle_strain'
                error('Muscle strain calbiration not currently supported.');
        end
    end
end

% PATH CONSTRAINTS
% Hill-equilibrium constraint
muscleData = DeGroote2016Muscle_FtildeState(a,Ftilde,dFtilde,splinestruct.LMT,splinestruct.VMT,params,input.auxdata.Fvparam,input.auxdata.Fpparam,input.auxdata.Faparam);

% Moments constraint
Topt = 150;
Tdiff = zeros(numColPoints,Ndof);
for dof = 1:Ndof
    T_exp=splinestruct.ID(:,dof);
    index_sel=(dof-1)*(NMuscles)+1:(dof-1)*(NMuscles)+NMuscles;
    T_sim=sum(muscleData.FT.*splinestruct.MA(:,index_sel),2) + Topt*aT(:,dof);
    Tdiff(:,dof) =  (T_exp-T_sim);
end

% DYNAMIC CONSTRAINTS
% Activation dynamics
dadt = ones(numColPoints,NMuscles);
for m = 1:NMuscles
    dadt(:,m) = ActivationDynamics(e(:,m),a(:,m),tauAct(m),tauDeact(m),input.auxdata.b);
end


% OBJECTIVE FUNCTION
costMuscles = fieldnames(terms);
EMGdiff = zeros(numColPoints, length(costMuscles));
emgScaleCount = 0;
for m = 1:length(costMuscles)
    if isfield(terms.(costMuscles{m}), 'costs')
        muscIdx = find(contains(MuscleNames, costMuscles{m}));
        
        if isfield(paramIndices.(costMuscles{m}), 'emgScale')
            idx = paramIndices.(costMuscles{m}).emgScale;
            emgScale = parameters(1,idx);
            emgScaleCount = emgScaleCount + 1;
        else
            emgScale = 1;
        end

        calibrationCosts = terms.(costMuscles{m}).costs;
        for c = 1:length(calibrationCosts)
            switch calibrationCosts{c}
                case 'emg'
                    EMGdiff(:,m) = (e(:,muscIdx) - emgScale*splinestruct.EMG(:,muscIdx));
                case 'fiber_length'
                    error('Fiber length cost not currently supported.');
                case 'fiber_velocity'
                    error('Fiber velocity cost not currently supported.');  
            end
        end        
    end
end

% Penalty on parameters deviating from nominal (except scale factors)
param_dev = 4*(parameters(1,1:(end-emgScaleCount))-1).^2;


% Outputs
phaseout.path = [Tdiff muscleData.err];
% Contraction dynamics is implicit
phaseout.dynamics = [dadt dFtilde];
w1 = 1000;
wAct = 0.1;
wParam = 0.1;
phaseout.integrand = w1.*sum(aT.^2,2) + sum(EMGdiff.^2,2) + wAct*sum(a.^2,2); %+ wParam*sum(param_dev,2)






