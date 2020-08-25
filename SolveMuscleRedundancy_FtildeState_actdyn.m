% SolveMuscleRedundancy_lMtildeState, version 0.1 (August 2016)
%
% This function solves the muscle redundancy problem in the leg using the
% direct collocation optimal control software GPOPS-II as described in De
% Groote F, Kinney AL, Rao AV, Fregly BJ. Evaluation of direct
% collocation optimal control problem formulations for solving the muscle
% redundancy problem. Annals of Biomedical Engineering (2016).
%
% Authors:  F. De Groote, M. Afschrift, A. Falisse
% Emails:   friedl.degroote@kuleuven.be
%           maarten.afschrift@kuleuven.be
%           antoine.falisse@kuleuven.be
%
% ----------------------------------------------------------------------- %
% This function uses the tendon force Ft as a state (see
% aforementioned publication for more details)
%
% INPUTS:
%           model_path: path to the .osim model
%           IK_path: path to the inverse kinematics results
%           ID_path: path to the inverse dynamics results
%           time: time window
%           OutPath: path to folder where results will be saved
%           Misc: structure of input data (see manual for more details)
%
% OUTPUTS:
%           Time: time window (as used when solving the optimal control
%           problem)
%           MExcitation: muscle excitation
%           MActivation: muscle activation
%           RActivation: activation of the reserve actuators
%           TForce_tilde: normalized tendon force
%           TForce: tendon force
%           lMtilde: normalized muscle fiber length
%           lM: muscle fiber length
%           MuscleNames: names of muscles
%           OptInfo: output of GPOPS-II
%           DatStore: structure with data used for solving the optimal
%           control problem
%
% ----------------------------------------------------------------------- %
%%

function [Time,MExcitation,MActivation,RActivation,TForcetilde,TForce,MuscleNames,MuscleData,OptInfo,DatStore]=SolveMuscleRedundancy_FtildeState_actdyn(model_path,IK_path,ID_path,time,OutPath,Misc)

%% ---------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% PART I: INPUTS FOR OPTIMAL CONTROL PROBLEM ---------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
% Based on study and cost function, decide which continuous and endpoint  % 
% functions to use ------------------------------------------------------ %
if ~isfield(Misc, 'study') || isempty(Misc.study)
    Misc.study = 'DeGroote2016/';
end
study = strsplit(Misc.study,'/');
switch study{1}
    case 'DeGroote2016'
        tag = '';
    case 'SoftExosuitDesign'
        tag = ['Exo' study{2}];
    otherwise
        error('Study not recognized')
end
% Cost Function
if ~isfield(Misc,'costfun') || isempty(Misc.costfun)
   Misc.costfun='Exc_Act';
end
if ~strcmp(Misc.costfun,'Default')
    tag = [tag '_' Misc.costfun];
end

% ----------------------------------------------------------------------- %
% Check for optional input arguments, see manual for details------------- %

% Default low-pass filter:
%   Butterworth order: 6
%   Cutoff frequency: 6Hz
% Inverse Dynamics
if ~isfield(Misc,'f_cutoff_ID') || isempty(Misc.f_cutoff_ID)
    Misc.f_cutoff_ID=6;
end
if ~isfield(Misc,'f_order_ID') || isempty(Misc.f_order_ID)
    Misc.f_order_ID=6;
end
% Muscle-tendon lengths
if ~isfield(Misc,'f_cutoff_lMT') || isempty(Misc.f_cutoff_lMT)
    Misc.f_cutoff_lMT=6;
end
if ~isfield(Misc,'f_order_lMT') || isempty(Misc.f_order_lMT)
    Misc.f_order_lMT=6;
end
% Moment arms
if ~isfield(Misc,'f_cutoff_dM') || isempty(Misc.f_cutoff_dM)
    Misc.f_cutoff_dM=6;
end
if ~isfield(Misc,'f_order_dM') || isempty(Misc.f_order_dM)
    Misc.f_order_dM=6;
end
% Inverse Kinematics
if ~isfield(Misc,'f_cutoff_IK') || isempty(Misc.f_cutoff_IK)
    Misc.f_cutoff_IK=6;
end
if ~isfield(Misc,'f_order_IK') || isempty(Misc.f_order_IK)
    Misc.f_order_IK=6;
end
% Mesh Frequency
if ~isfield(Misc,'Mesh_Frequency') || isempty(Misc.Mesh_Frequency)
   Misc.Mesh_Frequency=100;
end
% Activation dynamics ('explicit' or 'implicit')
if ~isfield(Misc,'actdyn') || isempty(Misc.actdyn)
   Misc.actdyn = 'explicit';          
end
% Walking speed
if ~isfield(Misc,'speed') || isempty(Misc.speed)
   Misc.speed = 1.0;          
end
% Device force level for Quinlivan study
if ~isfield(Misc,'exo_force_level') || isempty(Misc.exo_force_level)
   Misc.exo_force_level = -1;          
end
% Variable tendon stiffness
if ~isfield(Misc, 'tendonStiffnessCoeff') || isempty(Misc.tendonStiffnessCoeff)
    Misc.tendonStiffnessCoeff = 35;
end
% Modify individual tendon stiffnesses
if ~isfield(Misc, 'tendonStiffnessModifiers') || isempty(Misc.tendonStiffnessModifiers)
    Misc.tendonStiffnessModifiers = [];
end
% Modify individual passive muscle strain due to maximum isometric force (e0)
if ~isfield(Misc, 'muscleStrainModifiers') || isempty(Misc.muscleStrainModifiers)
    Misc.muscleStrainModifiers = [];
end
% Modify individual passive muscle force exponential shape factor (kpe)
if ~isfield(Misc, 'muscleShapeFactModifiers') || isempty(Misc.muscleShapeFactModifiers)
    Misc.muscleShapeFactModifiers = [];
end
% Modify individual muscle optimal fiber lengths (lMo)
if ~isfield(Misc, 'optimalFiberLengthModifiers') || isempty(Misc.optimalFiberLengthModifiers)
    Misc.optimalFiberLengthModifiers = [];
end
% Modify individual tendon slack lengths (lTs)
if ~isfield(Misc, 'tendonSlackLengthModifiers') || isempty(Misc.tendonSlackLengthModifiers)
    Misc.tendonSlackLengthModifiers = [];
end
% Modify individual pennation angles @ optimal fiber length (alphao)
if ~isfield(Misc, 'pennationAngleModifiers') || isempty(Misc.pennationAngleModifiers)
    Misc.pennationAngleModifiers = [];
end
if strcmp(study{2},'Topology')
    % ExoTopology: DOF's assisted by exoskeleton device
    if ~isfield(Misc, 'deviceDOFs') || isempty(Misc.deviceDOFs)
       Misc.deviceDOFs = []; 
    end
    % ExoTopology: fix moment arms to a constant value
    if ~isfield(Misc, 'fixMomentArms') || isempty(Misc.fixMomentArms)
       Misc.fixMomentArms = []; 
    end
    % ExoTopology: fix a torque parameter to a constant value
    if ~isfield(Misc, 'fixParams') || isempty(Misc.fixParams)
       Misc.fixParams = []; 
    end
    % ExoTopology: option to couple control signals between DOFs
    if ~isfield(Misc, 'coupledControl') || isempty(Misc.coupledControl)
       Misc.coupledControl = false; 
    end
end

% ----------------------------------------------------------------------- %
% Check that options are being specified correctly -----------------------%
if ~strcmp(study{2},'Topology')
    exoFlags = {'deviceDOFs','fixMomentArms', 'fixParams', ...
        'coupledControl', 'paramGuess'};
    
    for i = 1:length(exoFlags) 
        errmsg = [study{2} ': flag ' exoFlags{i} ' unused'];
        assert(~isfield(Misc, exoFlags{i}), errmsg)
    end
end

% ------------------------------------------------------------------------%
% Compute ID -------------------------------------------------------------%
if isempty(ID_path) || ~exist(ID_path,'file')
    disp('ID path was not specified or the file does not exist, computation ID started');
    if ~isfield(Misc,'Loads_path') || isempty(Misc.Loads_path) || ~exist(Misc.Loads_path,'file')
        error('External loads file was not specified or does not exist, please add the path to the external loads file: Misc.Loads_path');
    else
        %check the output path for the ID results
        if isfield(Misc,'ID_ResultsPath')
            [idpath,~]=fileparts(Misc.ID_ResultsPath);
            if ~isdir(idpath); mkdir(idpath); end
        else 
            % save results in the directory of the external loads
            [Lpath,name,~]=fileparts(Misc.Loads_path);
            Misc.ID_ResultsPath=fullfile(Lpath,name);
        end
        [ID_outPath,ID_outName,ext]=fileparts(Misc.ID_ResultsPath);
        output_settings=fullfile(ID_outPath,[ID_outName '_settings.xml']);
        Opensim_ID(model_path,[time(1)-0.1 time(2)+0.1],Misc.Loads_path,IK_path,ID_outPath,[ID_outName ext],output_settings);
        ID_path=[Misc.ID_ResultsPath '.sto'];
    end    
end

% ----------------------------------------------------------------------- %
% Muscle analysis ------------------------------------------------------- %

Misc.time=time;
MuscleAnalysisPath=fullfile(OutPath,'MuscleAnalysis'); if ~exist(MuscleAnalysisPath,'dir'); mkdir(MuscleAnalysisPath); end
disp('MuscleAnalysis Running .....');
OpenSim_Muscle_Analysis(IK_path,model_path,MuscleAnalysisPath,[time(1) time(end)])
disp('MuscleAnalysis Finished');
Misc.MuscleAnalysisPath=MuscleAnalysisPath;

% ----------------------------------------------------------------------- %
% Extract muscle information -------------------------------------------- %

% Get number of degrees of freedom (dofs), muscle-tendon lengths and moment
% arms for the selected muscles.
[~,Misc.trialName,~]=fileparts(IK_path);
if ~isfield(Misc,'MuscleNames_Input') || isempty(Misc.MuscleNames_Input)    
    Misc=getMuscles4DOFS(Misc);
end
[DatStore] = getMuscleInfo(IK_path,ID_path,Misc);

% ----------------------------------------------------------------------- %
% Solve the muscle redundancy problem using static optimization --------- %

% The solution of the static optimization is used as initial guess for the
% dynamic optimization
% Extract the muscle-tendon properties
[DatStore.params, DatStore.lOpt, DatStore.L_TendonSlack, DatStore.Fiso, ... 
    DatStore.PennationAngle, DatStore.metabolicParams] = ...
        ReadMuscleParameters(model_path,DatStore.MuscleNames);

% Modify tendon stiffnesses
fprintf('Muscles with modified tendon stiffness: \n')
for m = 1:DatStore.nMuscles
    muscle_name = DatStore.MuscleNames{m};
    if isfield(Misc.tendonStiffnessModifiers, muscle_name) 
        DatStore.params(6,m) = Misc.tendonStiffnessModifiers.(muscle_name);
        fprintf('--> %s tendon coefficient set to %f \n',muscle_name,DatStore.params(6,m))
    else
        DatStore.params(6,m) = 1;
    end
    if isfield(Misc.muscleStrainModifiers, muscle_name)
        DatStore.params(7,m) = Misc.muscleStrainModifiers.(muscle_name);
        fprintf('--> %s muscle strain set to %f \n',muscle_name,DatStore.params(7,m))
    else
        DatStore.params(7,m) = 1;
    end
    if isfield(Misc.muscleShapeFactModifiers, muscle_name)
        DatStore.params(8,m) = Misc.muscleShapeFactModifiers.(muscle_name);
        fprintf('--> %s muscle shape factor set to %f \n',muscle_name,DatStore.params(8,m))
    else
        DatStore.params(8,m) = 1;
    end
    if isfield(Misc.optimalFiberLengthModifiers, muscle_name)
        DatStore.params(9,m) = Misc.optimalFiberLengthModifiers.(muscle_name);
        fprintf('--> %s muscle optimal fiber length set to %f \n',muscle_name,DatStore.params(9,m))
    else
        DatStore.params(9,m) = 1;
    end
    if isfield(Misc.tendonSlackLengthModifiers, muscle_name)
        DatStore.params(10,m) = Misc.tendonSlackLengthModifiers.(muscle_name);
        fprintf('--> %s muscle tendon slack length set to %f \n',muscle_name,DatStore.params(10,m))
    else
        DatStore.params(10,m) = 1;
    end
    if isfield(Misc.pennationAngleModifiers, muscle_name)
        DatStore.params(11,m) = Misc.pennationAngleModifiers.(muscle_name);
        fprintf('--> %s muscle pennation angle set to %f \n',muscle_name,DatStore.params(11,m))
    else
        DatStore.params(11,m) = 1;
    end
end
fprintf('\n')

% Static optimization using IPOPT solver
% DatStore = SolveStaticOptimization_IPOPT(DatStore);


%% ---------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% PART II: OPTIMAL CONTROL PROBLEM FORMULATION -------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %

% Input arguments
DatStore.formulation = 'Ftilde';
auxdata.NMuscles = DatStore.nMuscles;   % number of muscles
auxdata.Ndof = DatStore.nDOF;           % humber of dofs
% DatStore.time = DatStore.time;          % time window
auxdata.ID = DatStore.T_exp;            % inverse dynamics
auxdata.params = DatStore.params;       % Muscle-tendon parameters
 % Parameters for calculating metabolic cost
auxdata.metabolicParams = DatStore.metabolicParams;
auxdata.speed = Misc.speed;             % walking speed

% ExoTopology study: handling possible device cases
% Create indicies for parameter array
model = org.opensim.modeling.Model(model_path);
state = model.initSystem();
model_mass = model.getTotalMass(state);
auxdata.model_mass = model_mass;
auxdata.numDeviceDOFs = 1;
config = ReadYaml('C:\Users\Nick\Projects\ExoTopology\exotopology\config.yaml');
norm_hip_max_torque = config.norm_hip_max_torque; % N-m/kg
norm_knee_max_torque = config.norm_knee_max_torque; % N-m/kg
norm_ankle_max_torque = config.norm_ankle_max_torque; % N-m/kg
if strcmp(study{2}, 'Topology')  
    minMomentArm = 0.05;
    maxMomentArm = 1.00;
    
    numExoParams = 0;
    paramsLower = [];
    paramsUpper = [];
    paramsGuess = [];
    isBidirectional = [];
    auxdata.Tmax_act_hip = norm_hip_max_torque * model_mass;
    auxdata.Tmax_act_knee = norm_knee_max_torque * model_mass;
    auxdata.Tmax_act_ankle = norm_ankle_max_torque * model_mass;
    auxdata.device.hip = 0;
    auxdata.device.knee = 0;
    auxdata.device.ankle = 0;
    for i = 1:length(Misc.deviceDOFs)
        dofInfo = split(Misc.deviceDOFs{i},'/');
        switch dofInfo{1}
            case 'hip'
                numExoParams = numExoParams + 1;
                auxdata.device.hip = numExoParams;
                paramsLower(numExoParams) = -1.0; %#ok<*AGROW>
                paramsUpper(numExoParams) = 1.0;
                isBidirectional(numExoParams) = 0;
                if ~isempty(Misc.fixMomentArms)
                    paramsLower(numExoParams) = Misc.fixMomentArms;
                    paramsUpper(numExoParams) = Misc.fixMomentArms;
                    paramsGuess(numExoParams) = Misc.fixMomentArms;
                end
                if length(dofInfo) > 1
                    switch dofInfo{2}
                        case 'flex'
                            if ~isempty(Misc.fixMomentArms)
                                paramsLower(numExoParams) = Misc.fixMomentArms;
                                paramsUpper(numExoParams) = Misc.fixMomentArms;
                                paramsGuess(numExoParams) = Misc.fixMomentArms;
                            else
                                paramsLower(numExoParams) = minMomentArm;
                                paramsUpper(numExoParams) = maxMomentArm;
                            end
                        case 'ext'
                            if ~isempty(Misc.fixMomentArms)
                                paramsLower(numExoParams) = -Misc.fixMomentArms;
                                paramsUpper(numExoParams) = -Misc.fixMomentArms;
                                paramsGuess(numExoParams) = -Misc.fixMomentArms;
                            else
                                paramsLower(numExoParams) = -maxMomentArm;
                                paramsUpper(numExoParams) = -minMomentArm;
                            end
                    end
                else
                    assert(~Misc.coupledControl, ...
                       'Cannot use bidirectional actuators in coupled control.');
                    isBidirectional(numExoParams) = 1;
                end
            case 'knee'
                numExoParams = numExoParams + 1;
                auxdata.device.knee = numExoParams;
                paramsLower(numExoParams) = -1.0;
                paramsUpper(numExoParams) = 1.0;
                isBidirectional(numExoParams) = 0;
                if ~isempty(Misc.fixMomentArms)
                    paramsLower(numExoParams) = Misc.fixMomentArms;
                    paramsUpper(numExoParams) = Misc.fixMomentArms;
                    paramsGuess(numExoParams) = Misc.fixMomentArms;
                end
                if length(dofInfo) > 1
                    switch dofInfo{2}
                        case 'ext'
                            if ~isempty(Misc.fixMomentArms)
                                paramsLower(numExoParams) = Misc.fixMomentArms;
                                paramsUpper(numExoParams) = Misc.fixMomentArms;
                                paramsGuess(numExoParams) = Misc.fixMomentArms;
                            else
                                paramsLower(numExoParams) = minMomentArm;
                                paramsUpper(numExoParams) = maxMomentArm;
                            end
                        case 'flex'
                            if ~isempty(Misc.fixMomentArms)
                                paramsLower(numExoParams) = -Misc.fixMomentArms;
                                paramsUpper(numExoParams) = -Misc.fixMomentArms;
                                paramsGuess(numExoParams) = -Misc.fixMomentArms;
                            else
                                paramsLower(numExoParams) = -maxMomentArm;
                                paramsUpper(numExoParams) = -minMomentArm;
                            end
                    end
                else
                    assert(~Misc.coupledControl, ...
                        'Cannot use bidirectional actuators in coupled control.');
                    isBidirectional(numExoParams) = 1;
                end
            case 'ankle'
                numExoParams = numExoParams + 1;
                auxdata.device.ankle = numExoParams;
                paramsLower(numExoParams) = -1.0;
                paramsUpper(numExoParams) = 1.0;
                isBidirectional(numExoParams) = 0;
                if ~isempty(Misc.fixMomentArms)
                    paramsLower(numExoParams) = Misc.fixMomentArms;
                    paramsUpper(numExoParams) = Misc.fixMomentArms;
                    paramsGuess(numExoParams) = Misc.fixMomentArms;
                end
                if length(dofInfo) > 1
                    switch dofInfo{2}
                        case 'dorsi'
                            if ~isempty(Misc.fixMomentArms)
                                paramsLower(numExoParams) = Misc.fixMomentArms;
                                paramsUpper(numExoParams) = Misc.fixMomentArms;
                                paramsGuess(numExoParams) = Misc.fixMomentArms;
                            else
                                paramsLower(numExoParams) = minMomentArm;
                                paramsUpper(numExoParams) = maxMomentArm;
                            end
                        case 'plantar'
                            if ~isempty(Misc.fixMomentArms)
                                paramsLower(numExoParams) = -Misc.fixMomentArms;
                                paramsUpper(numExoParams) = -Misc.fixMomentArms;
                                paramsGuess(numExoParams) = -Misc.fixMomentArms;
                            else
                                paramsLower(numExoParams) = -maxMomentArm;
                                paramsUpper(numExoParams) = -minMomentArm;
                            end
                    end
                else
                    assert(~Misc.coupledControl, ...
                        'Cannot use bidirectional actuators in coupled control.');
                    isBidirectional(numExoParams) = 1;
                end
        end
    end
    if ~Misc.coupledControl
        auxdata.numDeviceDOFs = numExoParams;
    end
    
    
    % Pass parameter index info to auxdata so it can be used in static
    % optimization initial guess
    auxdata.numExoParams = numExoParams;
    auxdata.paramsLower = paramsLower;
    auxdata.paramsUpper = paramsUpper;
    auxdata.isBidirectional = isBidirectional;
    
    auxdata.hip_DOF = strmatch('hip_flexion',DatStore.DOFNames);
    auxdata.knee_DOF = strmatch('knee_angle',DatStore.DOFNames);
    auxdata.ankle_DOF = strmatch('ankle_angle',DatStore.DOFNames);
    
    % Check knee coordinate convention 
    model = org.opensim.modeling.Model(model_path);
    coord_set = model.getCoordinateSet();
    knee_coord = coord_set.get('knee_angle_r');
    knee_range_min = (180/pi) * knee_coord.getRangeMin;
    knee_range_max = (180/pi) * knee_coord.getRangeMax;
    
    % Knee should increase joint angle during anterior swing
    auxdata.kneeAngleSign = 1;
    if (knee_range_min < -100) && (knee_range_max >= 0)
        auxdata.kneeAngleSign = 1;
    elseif (-10 <= knee_range_min && knee_range_min <= 0) && (knee_range_max > 100)
        auxdata.kneeAngleSign = -1;
    end
    
    auxdata.mod_name = Misc.mod_name;
end

% ADiGator works with 2D: convert 3D arrays to 2D structure (moment arms)
for i = 1:auxdata.Ndof
    auxdata.MA(i).Joint(:,:) = DatStore.dM(:,i,:);  % moment arms
end
auxdata.DOFNames = DatStore.DOFNames;   % names of dofs

tau_act = 0.015; auxdata.tauAct = tau_act * ones(1,auxdata.NMuscles);       % activation time constant (activation dynamics)
tau_deact = 0.06; auxdata.tauDeact = tau_deact * ones(1,auxdata.NMuscles);  % deactivation time constant (activation dynamics)
auxdata.b = 0.1;   

% Parameters of active muscle force-velocity characteristic
load ActiveFVParameters.mat
Fvparam(1) = 1.475*ActiveFVParameters(1); Fvparam(2) = 0.25*ActiveFVParameters(2);
Fvparam(3) = ActiveFVParameters(3) + 0.75; Fvparam(4) = ActiveFVParameters(4) - 0.027;
auxdata.Fvparam = Fvparam;

% Parameters of active muscle force-length characteristic
load Faparam.mat                            
auxdata.Faparam = Faparam;

% Parameters of passive muscle force-length characteristic, and tendon
% characteristic
e0 = 0.6*DatStore.params(7,:); 
kpe = 4*DatStore.params(8,:); 
t50 = exp(kpe .* (0.2 - 0.10e1) ./ e0);
pp1 = (t50 - 0.10e1); t7 = exp(kpe); pp2 = (t7 - 0.10e1);
auxdata.Fpparam = [pp1;pp2;ones(1,length(pp1))*Misc.tendonStiffnessCoeff];

% Problem bounds 
% --------------
% Time bounds
t0 = DatStore.time(1); 
tf = DatStore.time(end);
bounds.phase.initialtime.lower = t0; 
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tf; 
bounds.phase.finaltime.upper = tf;
auxdata.initialtime = t0;
auxdata.finaltime = tf;

% Controls bounds
N = length(DatStore.time);
vA_min = -1/100; vA_max = 1/100;  % bounds on derivative of muscle activation (if implicit)
e_min = 0; e_max = 1;             % bounds on excitation (if explicit)
dF_min = -100; dF_max = 100;      % bounds on derivative of normalzied tendon force
actdyn = Misc.actdyn;
auxdata.actdyn = actdyn;
if strcmp(actdyn, 'implicit')
    umin = vA_min./auxdata.tauDeact;
    umax = vA_max./auxdata.tauAct;
    uguess = zeros(N,auxdata.NMuscles);
elseif strcmp(actdyn, 'explicit')
    umin = e_min*ones(1,auxdata.NMuscles); 
    umax = e_max*ones(1,auxdata.NMuscles);
    uguess = 0.1*ones(N,auxdata.NMuscles);
end
dFMin = dF_min*ones(1,auxdata.NMuscles);
dFMax = dF_max*ones(1,auxdata.NMuscles);
aTmin = -1*ones(1,auxdata.Ndof); 
aTmax = 1*ones(1,auxdata.Ndof);
if strcmp(study{2},'Topology')
    aDmin = zeros(1, auxdata.numDeviceDOFs); 
    aDmax = ones(1, auxdata.numDeviceDOFs);
    for i = 1:auxdata.numDeviceDOFs
        if auxdata.isBidirectional(i)
            aDmin(i) = -1;
        end
    end
    control_bounds_lower = [umin aTmin dFMin aDmin];
    control_bounds_upper = [umax aTmax dFMax aDmax];
else
    control_bounds_lower = [umin aTmin dFMin];
    control_bounds_upper = [umax aTmax dFMax];
end
bounds.phase.control.lower = control_bounds_lower; 
bounds.phase.control.upper = control_bounds_upper;

% States bounds
a_min = 0; a_max = 1;             % bounds on muscle activation
F_min = 0; F_max = 5;             % bounds on normalized tendon force
actMin = a_min*ones(1,auxdata.NMuscles); 
actMax = a_max*ones(1,auxdata.NMuscles);
F0min = F_min*ones(1,auxdata.NMuscles);
F0max = F_max*ones(1,auxdata.NMuscles);
Ffmin = F_min*ones(1,auxdata.NMuscles);
Ffmax = F_max*ones(1,auxdata.NMuscles);
FMin = F_min*ones(1,auxdata.NMuscles);
FMax = F_max*ones(1,auxdata.NMuscles);
if strcmp(study{2}, 'Topology')
    bounds.phase.initialstate.lower = [actMin, F0min,];
    bounds.phase.initialstate.upper = [actMax, F0max];
    bounds.phase.state.lower = [actMin, FMin];
    bounds.phase.state.upper = [actMax, FMax];
    bounds.phase.finalstate.lower = [actMin, Ffmin];
    bounds.phase.finalstate.upper = [actMax, Ffmax];
else
    bounds.phase.initialstate.lower = [actMin, F0min];
    bounds.phase.initialstate.upper = [actMax, F0max];
    bounds.phase.state.lower = [actMin, FMin];
    bounds.phase.state.upper = [actMax, FMax];
    bounds.phase.finalstate.lower = [actMin, Ffmin];
    bounds.phase.finalstate.upper = [actMax, Ffmax];
end

% Integral bounds
bounds.phase.integral.lower = 0;
bounds.phase.integral.upper = 10000*(tf-t0);

% Parameter bounds
if strcmp(study{2},'Topology')
    % Parameter variable in the optimization problem always lie in the range 
    % [-0.5 0.5], and should be converted as necessary in the continuous function
    % for to make the correct computations.
    bounds.parameter.lower = -0.5*ones(size(auxdata.paramsLower));
    bounds.parameter.upper = 0.5*ones(size(auxdata.paramsUpper));
end

% Path constraints
HillEquil = zeros(1, auxdata.NMuscles);
ID_bounds = zeros(1, auxdata.Ndof);
act1_lower = zeros(1, auxdata.NMuscles);
act1_upper = inf*ones(1, auxdata.NMuscles);
act2_lower = -inf*ones(1, auxdata.NMuscles);
act2_upper = 1*ones(1, auxdata.NMuscles)./auxdata.tauAct;
if strcmp(actdyn, 'implicit')
    bounds.phase.path.lower = [ID_bounds,HillEquil,act1_lower,act2_lower];
    bounds.phase.path.upper = [ID_bounds,HillEquil,act1_upper,act2_upper];
elseif strcmp(actdyn, 'explicit')
    bounds.phase.path.lower = [ID_bounds,HillEquil];
    bounds.phase.path.upper = [ID_bounds,HillEquil];
end

% Eventgroup
% Impose mild periodicity
periodcity = 1;
pera_lower = -periodcity * ones(1, auxdata.NMuscles);
pera_upper = periodcity * ones(1, auxdata.NMuscles);
perFtilde_lower = -periodcity * ones(1,auxdata.NMuscles);
perFtilde_upper = periodcity * ones(1,auxdata.NMuscles);
if strcmp(study{2},'Topology')
    bounds.eventgroup.lower = [pera_lower perFtilde_lower];
    bounds.eventgroup.upper = [pera_upper perFtilde_upper];
else
    bounds.eventgroup.lower = [pera_lower perFtilde_lower];
    bounds.eventgroup.upper = [pera_upper perFtilde_upper];
end

% Initial guesses

% Time guess
if strcmp(study{2}, 'Topology')
    % Load unassisted solution to use as initial guess below.
    unassisted = load(Misc.mrs_path);
    solution = unassisted.OptInfo.result.solution;
    guess.phase.time = unassisted.Time;
    N = size(unassisted.Time, 1);
else
    guess.phase.time = DatStore.time;
end


% Control guess
dF_guess = zeros(N,auxdata.NMuscles);
if strcmp(study{2},'Topology')
    % Use controls from unassisted solution in guess.
    if strcmp(actdyn, 'implicit')
        uguess = (1/100)*unassisted.vA;
    elseif strcmp(actdyn, 'explicit')
        uguess = unassisted.MExcitation;
    end
    dF_guess = solution.phase.control(:,auxdata.NMuscles+auxdata.Ndof+1:auxdata.NMuscles+auxdata.Ndof+auxdata.NMuscles);
    control_guess = [uguess unassisted.RActivation dF_guess 0.1*ones(N,auxdata.numDeviceDOFs)];
else
    control_guess = [uguess zeros(N,auxdata.Ndof) dF_guess];

end
guess.phase.control = control_guess;

% State guess
if strcmp(study{2}, 'Topology') 
    a_guess = unassisted.MActivation;
    F_guess = unassisted.TForcetilde;
    aD_guess = 0.1*ones(size(a_guess,1), auxdata.numDeviceDOFs);
    state_guess = [a_guess F_guess];
else
    a_guess = 0.1*ones(N,auxdata.NMuscles);
    F_guess = 0.1*ones(N,auxdata.NMuscles); 
    state_guess = [a_guess F_guess];
end

guess.phase.state = state_guess;

% Integral guess
if strcmp(study{2}, 'Topology') 
    guess.phase.integral = solution.phase.integral;
else 
    guess.phase.integral = 0;
end

% Parameter guess
if strcmp(study{2},'Topology')
    guess.parameter = (bounds.parameter.upper - bounds.parameter.lower) / 2;
end

DatStore.T_exo = zeros(length(DatStore.time),auxdata.Ndof);

% Spline structures
for dof = 1:auxdata.Ndof
    for m = 1:auxdata.NMuscles       
        auxdata.JointMASpline(dof).Muscle(m) = spline(DatStore.time,auxdata.MA(dof).Joint(:,m));       
    end
    auxdata.JointIDSpline(dof) = spline(DatStore.time,DatStore.T_exp(:,dof));
    auxdata.JointEXOSpline(dof) = spline(DatStore.time,DatStore.T_exo(:,dof));
    auxdata.JointIKSpline(dof) = spline(DatStore.time,DatStore.q_exp(:,dof));
end

for m = 1:auxdata.NMuscles
    auxdata.LMTSpline(m) = spline(DatStore.time,DatStore.LMT(:,m));
end

% GPOPS setup
setup.name = 'DynamicOptimization_Ftildestate_vA';
setup.auxdata = auxdata;
setup.bounds = bounds;
setup.guess = guess;
setup.nlp.solver = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.nlp.ipoptoptions.tolerance = 1e-3;
setup.nlp.ipoptoptions.maxiterations = 10000;
setup.derivatives.supplier = 'sparseCD';
setup.derivatives.derivativelevel = 'first';
setup.derivatives.dependencies = 'sparse';
setup.scales.method = 'none';
setup.mesh.method = 'hp-PattersonRao';
setup.mesh.tolerance = 1e-3;
setup.mesh.maxiterations = 0;
setup.mesh.colpointsmin = 5;
setup.mesh.colpointsmax = 10;
setup.method = 'RPM-integration';
setup.displaylevel = 2;
NMeshIntervals = round((tf-t0)*Misc.Mesh_Frequency);
setup.mesh.phase.colpoints = 5*ones(1,NMeshIntervals);
setup.mesh.phase.fraction = (1/(NMeshIntervals))*ones(1,NMeshIntervals);
setup.functions.continuous = str2func(['continous_Ftilde_vA' tag]);
setup.functions.endpoint = str2func(['endpoint_Ftilde' tag]);
    
% ADiGator setup
persistent splinestruct
input.auxdata = auxdata;

% Path locking for Linux platforms
if isunix
    pathLock='/tmp/adigator3.lock'
    % Try to create and lock this file.
    if ismac
        lockfilecommand = 'dotlockfile'; % Get from homebrew.
    else
        lockfilecommand = 'lockfile';
    end
    if ~system(sprintf('%s %s',lockfilecommand, pathLock))
        % We succeeded, so perform some task which needs to be serialized.
        tdummy = guess.phase.time;
        splinestruct = SplineInputData(tdummy,input);
        splinenames = fieldnames(splinestruct);
        for Scount = 1:length(splinenames)
            secdim = size(splinestruct.(splinenames{Scount}),2);
            splinestructad.(splinenames{Scount}) = adigatorCreateAuxInput([Inf,secdim]);
            splinestruct.(splinenames{Scount}) = zeros(0,secdim);
        end
        setup.auxdata.splinestruct = splinestructad;
        if strcmp(setup.derivatives.supplier, 'adigator')
        	adigatorGenFiles4gpops2(setup)
        end
        
        % Now remove the lockfile
        system(sprintf('rm -f %s',pathLock));
    end
    
elseif ispc
    
    lockDir = 'C:\Users\Nick\tmp\adigatorLock\';
    pathLock=fullfile(lockDir, 'lockFile.mat');
    
    % If lock file exists, wait until it is deleted by a parallel process
    while true
        pause(randi(5,1)) % wait 5 seconds
        if ~(exist(pathLock,'file')==2)
           break 
        end
    end
    
    % Create a new lock file for this process
    emptyVar = [];
    save(pathLock, 'emptyVar')
    
    % Perform serialzied task
    tdummy = guess.phase.time;
    splinestruct = SplineInputData(tdummy,input);
    splinenames = fieldnames(splinestruct);
    for Scount = 1:length(splinenames)
        secdim = size(splinestruct.(splinenames{Scount}),2);
        splinestructad.(splinenames{Scount}) = adigatorCreateAuxInput([Inf,secdim]);
        splinestruct.(splinenames{Scount}) = zeros(0,secdim);
    end
    setup.auxdata.splinestruct = splinestructad;
    if strcmp(setup.derivatives.supplier, 'adigator')
        adigatorGenFiles4gpops2(setup)
    end

    % Remove the lockfile
    system(sprintf('del %s',pathLock))

else
    error('Platform unknown.')
end

if strcmp(actdyn, 'implicit')
    wrap_continuous = 'Wrap4continous_Ftilde_vA';
    continuous = 'continous_Ftilde_vA';
elseif strcmp(actdyn, 'explicit')
    wrap_continuous = 'Wrap4continous_Ftilde';
    continuous = 'continous_Ftilde';
end

setup.functions.continuous = str2func([wrap_continuous tag]);
setup.adigatorgrd.continuous = str2func([continuous tag 'GrdWrap']);
setup.adigatorgrd.endpoint   = str2func(['endpoint_Ftilde' tag 'ADiGatorGrd']);
setup.adigatorhes.continuous = str2func(['continous_Ftilde_vA' tag 'HesWrap']);
setup.adigatorhes.endpoint   = str2func([continuous tag 'ADiGatorHes']);

%% ---------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% PART III: SOLVE OPTIMAL CONTROL PROBLEM ------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %

setup.auxdata.regularizationWeight = 1;
regOutput = gpops2(setup);

setup.auxdata.regularizationWeight = 1e-2;
setup.guess = regOutput.result.solution;
output = gpops2(setup);

MuscleNames = DatStore.MuscleNames;
res = output.result.solution.phase(1);
Time = res.time;
MActivation = res.state(:,1:auxdata.NMuscles);
TForcetilde = res.state(:,auxdata.NMuscles+1:auxdata.NMuscles*2);
TForce = TForcetilde.*(ones(size(Time))*DatStore.Fiso);
dTForcetilde = 10*res.control(:,auxdata.NMuscles+auxdata.Ndof+1:auxdata.NMuscles+auxdata.Ndof+auxdata.NMuscles); 
if strcmp(actdyn, 'implicit')
    vA = 100*res.control(:,1:auxdata.NMuscles);
    MExcitation = computeExcitationRaasch(MActivation, vA, auxdata.tauDeact, auxdata.tauAct);
elseif strcmp(actdyn, 'explicit')
    MExcitation = res.control(:,1:auxdata.NMuscles);
end

RActivation = res.control(:,auxdata.NMuscles+1:auxdata.NMuscles+auxdata.Ndof);
OptInfo = output;

% Calculate muscle metabolic rates
mat.Time = Time;
mat.DatStore = DatStore;
mat.OptInfo = OptInfo;
mat.MuscleNames = MuscleNames;


[UmbergerKoelewijn2018, MinettiAlexander1997] = ...
        calcWholeBodyMetabolicRate(model, mat);
MetabolicInfo.UmbergerKoelewijn2018 = UmbergerKoelewijn2018;
MetabolicInfo.MinettiAlexander1997 = MinettiAlexander1997;
DatStore.MetabolicInfo = MetabolicInfo;

% Tendon force from lMtilde
% Interpolation lMT
lMTinterp = interp1(DatStore.time,DatStore.LMT,Time);
for m = 1:auxdata.NMuscles
    LMTSpline(m) = spline(Time,lMTinterp(:,m));
    [LMT(:,m),VMT(:,m),~] = SplineEval_ppuval(LMTSpline(m),Time,1);
end

MuscleData = DeGroote2016Muscle_FtildeState(MActivation, TForcetilde, ...
    dTForcetilde, LMT, VMT, auxdata.params, auxdata.Fvparam, auxdata.Fpparam, ... 
    auxdata.Faparam);

if strcmp(study{2},'Topology')
    [DatStore.ExoTorques, DatStore.MomentArms] = ...
        calcExoTorques_Ftilde_vAExoTopology(OptInfo, DatStore);
end
    
end

