function phaseout = continous_FtildeParamCalMultiPhase(input)

% Get input data
NMuscles        = input.auxdata.NMuscles;
Ndof            = input.auxdata.Ndof;
tauAct          = input.auxdata.tauAct;
tauDeact        = input.auxdata.tauDeact;
params          = input.auxdata.params;

% Get parameters
parameters = input.phase(1).parameter;

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
                error('cannot optimize muscle strain');
%                 params(7,muscIdx) = paramVal;
        end
    end
end

for p = 1:input.auxdata.numPhases
    
    numColPoints = size(input.phase(p).state,1);
    LMT = input.auxdata.splinestruct(p).p.LMT;
    VMT = input.auxdata.splinestruct(p).p.VMT;
    ID = input.auxdata.splinestruct(p).p.ID;
    MA = input.auxdata.splinestruct(p).p.MA;
    EMG = input.auxdata.splinestruct(p).p.EMG;
    
    % Get controls
    e        = input.phase(p).control(:,1:NMuscles);
    aT       = input.phase(p).control(:,NMuscles+1:NMuscles+Ndof);
    dFtilde  = 10*input.phase(p).control(:,NMuscles+Ndof+1:end);
    
    % Get states
    a       = input.phase(p).state(:,1:NMuscles);
    Ftilde  = input.phase(p).state(:,NMuscles+1:end);
    
    % PATH CONSTRAINTS
    % Hill-equilibrium constraint
    muscleData = DeGroote2016Muscle_FtildeState(a,Ftilde,dFtilde,LMT,VMT,params,input.auxdata.Fvparam,input.auxdata.Fpparam,input.auxdata.Faparam);
    
    % Moments constraint
    Topt = 150;
    Tdiff = zeros(numColPoints,Ndof);
    for dof = 1:Ndof
        T_exp=ID(:,dof);
        index_sel=(dof-1)*(NMuscles)+1:(dof-1)*(NMuscles)+NMuscles;
        T_sim=sum(muscleData.FT.*MA(:,index_sel),2) + Topt*aT(:,dof);
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
    for m = 1:length(costMuscles)
        if isfield(terms.(costMuscles{m}), 'costs')
            muscIdx = find(contains(MuscleNames, costMuscles{m}));
           
            calibrationCosts = terms.(costMuscles{m}).costs;
            for c = 1:length(calibrationCosts)
                switch calibrationCosts{c}
                    case 'emg'
                        EMGdiff(:,m) = (e(:,muscIdx) - EMG(:,muscIdx));
                    case 'fiber_length'
                        error('Fiber length cost not currently supported.');
                    case 'fiber_velocity'
                        error('Fiber velocity cost not currently supported.');
                end
            end
        end
    end
        
    % Outputs
    phaseout(p).path = [Tdiff muscleData.err];
    % Contraction dynamics is implicit
    phaseout(p).dynamics = [dadt dFtilde];
    w1 = 1000;
    wAct = 0.01;
    wPass = 1;
    phaseout(p).integrand = w1.*sum(aT.^2,2) + sum(EMGdiff.^2,2) + wAct*sum(a.^2,2) + wPass*sum(muscleData.fpe.^2, 2);
    
end






