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
% Subcase
if ~isfield(Misc,'subcase') || isempty(Misc.subcase)
   Misc.subcase = ''; 
else
   tag = [tag '_' Misc.subcase];
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
% ExoTopology: DOF's assisted by passive device
if ~isfield(Misc, 'passiveDOFs') || isempty(Misc.passiveDOFs)
   Misc.passiveDOFs = []; 
end
% ExoTopology: DOF's assisted by active device
if ~isfield(Misc, 'activeDOFs') || isempty(Misc.activeDOFs)
   Misc.activeDOFs = []; 
end
% ExoTopology: fix moment arms to a constant value
if ~isfield(Misc, 'fixMomentArms') || isempty(Misc.fixMomentArms)
   Misc.fixMomentArms = []; 
end
% ExoTopology: option to set individual control signals for each DOF
if ~isfield(Misc, 'mult_controls') || isempty(Misc.mult_controls)
   Misc.mult_controls = false; 
end
% ExoTopology: option to fix gains on experimental torque controls to be
% the same across DOFs
if ~isfield(Misc, 'same_torque_gain') || isempty(Misc.same_torque_gain)
   Misc.same_torque_gain = false; 
end
% Shift prescribed exoskeleton torque peaks to match ID peaks
if ~isfield(Misc, 'shift_exo_peaks') || isempty(Misc.shift_exo_peaks)
   Misc.shift_exo_peaks = false; 
end
% % Match device power to specified value
% if ~isfield(Misc, 'powerMatchType') || isempty(Misc.powerMatchType)
%    Misc.powerMatchType = [];
% end
% if ~isfield(Misc, 'powerMatchValue') || isempty(Misc.powerMatchValue)
%    Misc.powerMatchValue = []; 
% end
% ExoTopology (ActParam): guess for Zhang2017 parameterization
if ~isfield(Misc, 'paramGuess') || isempty(Misc.paramGuess)
   Misc.paramGuess = []; 
end

% ----------------------------------------------------------------------- %
% Check that options are being specified correctly -----------------------%
if ~strcmp(study{2},'Topology')
    errmsg = [study{2} ': active device DOFs unused'];
    assert(isempty(Misc.activeDOFs), errmsg)
    
    errmsg = [study{2} ': passive device DOFs unused'];
    assert(isempty(Misc.passiveDOFs), errmsg)
end
% if ~isempty(Misc.powerMatchValue) && isempty(Misc.powerMatchType)
%    error(['Power type (average positive power, average net power, etc) not ' ...
%           'specified.'])
% end
% if isempty(Misc.powerMatchValue) && ~isempty(Misc.powerMatchType)
%    error('Power value not specified for match type %s.', Misc.powerMatchType)
% end

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
[DatStore.params,DatStore.lOpt,DatStore.L_TendonSlack,DatStore.Fiso,DatStore.PennationAngle,DatStore.metabolicParams]=ReadMuscleParameters(model_path,DatStore.MuscleNames);

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
DatStore = SolveStaticOptimization_IPOPT(DatStore);


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
auxdata.metabolicParams = DatStore.metabolicParams; % Parameters for calculating metabolic cost

% ExoTopology study: handling possible active + passive device cases
% Create indicies for parameter array
model = org.opensim.modeling.Model(model_path);
state = model.initSystem();
model_mass = model.getTotalMass(state);
auxdata.model_mass = model_mass;
auxdata.numActiveDOFs = 1;

if strcmp(study{2}, 'Topology')
    auxdata.hasActiveDevice = false;
    auxdata.hasPassiveDevice = false;
    maxMomentArm = 1.00;
    
    numExoParams = 0;
    paramsLower = [];
    paramsUpper = [];
    % Active device indicies
    if ~isempty(Misc.activeDOFs)
        auxdata.hasActiveDevice = true;
        auxdata.Tmax_act = 1.2*max(max(abs(DatStore.T_exp)));
        auxdata.active.hip = 0;
        auxdata.active.knee = 0;
        auxdata.active.ankle = 0;
        for i = 1:length(Misc.activeDOFs)
            dofInfo = split(Misc.activeDOFs{i},'/');
            switch dofInfo{1}
                case 'hip'
                    numExoParams = numExoParams + 1;
                    auxdata.active.hip = numExoParams;
                    paramsLower(numExoParams) = -1.0; %#ok<*AGROW>
                    paramsUpper(numExoParams) = 1.0;
                    if length(dofInfo) > 1
                        switch dofInfo{2}
                            case 'flex'
                                if ~isempty(Misc.fixMomentArms)
                                    paramsLower(numExoParams) = Misc.fixMomentArms;
                                    paramsUpper(numExoParams) = Misc.fixMomentArms;
                                else
                                    paramsLower(numExoParams) = 0.1;
                                    paramsUpper(numExoParams) = maxMomentArm;
                                end
                            case 'ext'
                                if ~isempty(Misc.fixMomentArms)
                                    paramsLower(numExoParams) = -Misc.fixMomentArms;
                                    paramsUpper(numExoParams) = -Misc.fixMomentArms;
                                else
                                    paramsLower(numExoParams) = -maxMomentArm;
                                    paramsUpper(numExoParams) = -0.1;
                                end
                        end
                    end
                case 'knee'
                    numExoParams = numExoParams + 1;
                    auxdata.active.knee = numExoParams;
                    paramsLower(numExoParams) = -1.0;
                    paramsUpper(numExoParams) = 1.0;
                    if length(dofInfo) > 1
                        switch dofInfo{2}
                            case 'ext'
                                if ~isempty(Misc.fixMomentArms)
                                    paramsLower(numExoParams) = Misc.fixMomentArms;
                                    paramsUpper(numExoParams) = Misc.fixMomentArms;
                                else
                                    paramsLower(numExoParams) = 0.1;
                                    paramsUpper(numExoParams) = maxMomentArm;
                                end
                            case 'flex'
                                if ~isempty(Misc.fixMomentArms)
                                    paramsLower(numExoParams) = -Misc.fixMomentArms;
                                    paramsUpper(numExoParams) = -Misc.fixMomentArms;
                                else
                                    paramsLower(numExoParams) = -maxMomentArm;
                                    paramsUpper(numExoParams) = -0.1;
                                end
                        end
                    end
                case 'ankle'
                    numExoParams = numExoParams + 1;
                    auxdata.active.ankle = numExoParams;
                    paramsLower(numExoParams) = -1.0;
                    paramsUpper(numExoParams) = 1.0;
                    if length(dofInfo) > 1
                        switch dofInfo{2}
                            case 'dorsi'
                                if ~isempty(Misc.fixMomentArms)
                                    paramsLower(numExoParams) = Misc.fixMomentArms;
                                    paramsUpper(numExoParams) = Misc.fixMomentArms;
                                else
                                    paramsLower(numExoParams) = 0.1;
                                    paramsUpper(numExoParams) = maxMomentArm;
                                end
                            case 'plantar'
                                if ~isempty(Misc.fixMomentArms)
                                    paramsLower(numExoParams) = -Misc.fixMomentArms;
                                    paramsUpper(numExoParams) = -Misc.fixMomentArms;
                                else
                                    paramsLower(numExoParams) = -maxMomentArm;
                                    paramsUpper(numExoParams) = -0.1;
                                end
                        end
                    end
            end
        end
        if Misc.mult_controls
            auxdata.numActiveDOFs = numExoParams;
        end
    end
    % Passive device indicies
%     if ~isempty(Misc.passiveDOFs)
%         auxdata.hasPassiveDevice = true;
%         auxdata.passiveStiffness = 1250*auxdata.model_mass; % k ~= 100 kN/m for 80kg subject (van den Bogert 2003) 
%         auxdata.passive.hip = 0;
%         auxdata.passive.knee = 0;
%         auxdata.passive.ankle = 0;
%         for i = 1:length(Misc.passiveDOFs)
%             dofInfo = split(Misc.passiveDOFs{i},'/');
%             switch dofInfo{1}
%                 case 'hip'
%                     numExoParams = numExoParams + 1;
%                     auxdata.passive.hip = numExoParams;
%                     paramsLower(numExoParams) = -0.10;
%                     paramsUpper(numExoParams) = 0.10;
%                     if length(dofInfo) > 1
%                         switch dofInfo{2}
%                             case 'flex'
%                                 if ~isempty(Misc.fixMomentArms)
%                                     paramsLower(numExoParams) = Misc.fixMomentArms;
%                                     paramsUpper(numExoParams) = Misc.fixMomentArms;
%                                 else
%                                     paramsLower(numExoParams) = 0;
%                                     paramsUpper(numExoParams) = maxMomentArm;
%                                 end
%                             case 'ext'
%                                 if ~isempty(Misc.fixMomentArms)
%                                     paramsLower(numExoParams) = -Misc.fixMomentArms;
%                                     paramsUpper(numExoParams) = -Misc.fixMomentArms;
%                                 else
%                                     paramsLower(numExoParams) = -maxMomentArm;
%                                     paramsUpper(numExoParams) = 0;
%                                 end
%                         end
%                     end
%                 case 'knee'
%                     numExoParams = numExoParams + 1;
%                     auxdata.passive.knee = numExoParams;
%                     paramsLower(numExoParams) = -0.10;
%                     paramsUpper(numExoParams) = 0.10;
%                     if length(dofInfo) > 1
%                         switch dofInfo{2}
%                             case 'ext'
%                                 if ~isempty(Misc.fixMomentArms)
%                                     paramsLower(numExoParams) = Misc.fixMomentArms;
%                                     paramsUpper(numExoParams) = Misc.fixMomentArms;
%                                 else
%                                     paramsLower(numExoParams) = 0;
%                                     paramsUpper(numExoParams) = maxMomentArm;
%                                 end
%                             case 'flex'
%                                 if ~isempty(Misc.fixMomentArms)
%                                     paramsLower(numExoParams) = -Misc.fixMomentArms;
%                                     paramsUpper(numExoParams) = -Misc.fixMomentArms;
%                                 else
%                                     paramsLower(numExoParams) = -maxMomentArm;
%                                     paramsUpper(numExoParams) = 0;
%                                 end
%                         end
%                     end
%                     
%                 case 'ankle'
%                     numExoParams = numExoParams + 1;
%                     auxdata.passive.ankle = numExoParams;
%                     paramsLower(numExoParams) = -0.10;
%                     paramsUpper(numExoParams) = 0.10;
%                     if length(dofInfo) > 1
%                         switch dofInfo{2}
%                             case 'dorsi'
%                                 if ~isempty(Misc.fixMomentArms)
%                                     paramsLower(numExoParams) = Misc.fixMomentArms;
%                                     paramsUpper(numExoParams) = Misc.fixMomentArms;
%                                 else
%                                     paramsLower(numExoParams) = 0;
%                                     paramsUpper(numExoParams) = maxMomentArm;
%                                 end
%                             case 'plantar'
%                                 if ~isempty(Misc.fixMomentArms)
%                                     paramsLower(numExoParams) = -Misc.fixMomentArms;
%                                     paramsUpper(numExoParams) = -Misc.fixMomentArms;
%                                 else
%                                     paramsLower(numExoParams) = -maxMomentArm;
%                                     paramsUpper(numExoParams) = 0;
%                                 end
%                         end
%                     end
%             end
%         end
%         % Extra index for exotendon slack length
%         numExoParams = numExoParams + 1;
%         paramsLower(numExoParams) = 0.75;
%         paramsUpper(numExoParams) = 1.25;
%         auxdata.passive.slack_length = numExoParams;
%     end
    
    % Pass parameter index info to auxdata so it can be used in static
    % optimization initial guess
    auxdata.numExoParams = numExoParams;
    auxdata.paramsLower = paramsLower;
    auxdata.paramsUpper = paramsUpper;
    auxdata.subcase = Misc.subcase;
    
    % Pass whether or not exoskeleton torques should be shifted based on 
    % net joint moment timings
    auxdata.shift_exo_peaks = Misc.shift_exo_peaks;
    
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
    
end

% For the ActParam subcase, set indices for peak torque, peak time,
% rise time, and fall time (ref. Zhang et al. 2017)
if strcmp(Misc.subcase, 'ActParam')
    
    for i = 1:length(Misc.activeDOFs)
        dofInfo = split(Misc.activeDOFs{i},'/');
        for dof = 1:auxdata.Ndof
            if contains(DatStore.DOFNames{dof}, dofInfo{1})
                % Max possible torque is smallest ID torque peak
                if max(abs(DatStore.T_exp(:,dof))) < auxdata.Tmax_act
                    auxdata.Tmax_act = max(abs(DatStore.T_exp(:,dof)));
                end               
            end
        end
    end
       
    % peak torque
    peakTorqueGuess = Misc.paramGuess(1) / auxdata.Tmax_act;
    numExoParams = numExoParams + 1;
    auxdata.active.params.peak_torque.idx = numExoParams;
    paramsLower(numExoParams) = peakTorqueGuess*0.75;
    paramsUpper(numExoParams) = peakTorqueGuess*1.25;

    % peak time
    peakTimeGuess = (Misc.paramGuess(2)-time(1))/(time(2)-time(1));
    numExoParams = numExoParams + 1;
    auxdata.active.params.peak_time.idx = numExoParams;
    paramsLower(numExoParams) = peakTimeGuess*0.75;
    paramsUpper(numExoParams) = peakTimeGuess*1.25;
    
    % rise time
    riseTimeGuess = Misc.paramGuess(3) / (time(2)-time(1));
    numExoParams = numExoParams + 1;
    auxdata.active.params.rise_time.idx = numExoParams;
    paramsLower(numExoParams) = riseTimeGuess*0.75;
    paramsUpper(numExoParams) = riseTimeGuess*1.25;
    
    % fall time
    fallTimeGuess = Misc.paramGuess(4) / (time(2)-time(1));
    numExoParams = numExoParams + 1;
    auxdata.active.params.fall_time.idx = numExoParams;
    paramsLower(numExoParams) = fallTimeGuess*0.75;
    paramsUpper(numExoParams) = fallTimeGuess*1.25;
    
    % Overwrite auxdata fields from non-parameterized cases
    auxdata.numExoParams = numExoParams;
    auxdata.paramsLower = paramsLower;
    auxdata.paramsUpper = paramsUpper;
end

% TODO: figure out what to do with this case
% if strcmp(Misc.subcase, 'Exp')
%     
%     if isempty(Misc.fixMomentArms)
%         error('Optimizing moment arms not supported for Exp subcase');
%     end
%     
%     numExoParams = 0;
%     paramsLower = [];
%     paramsUpper = [];  
%     
%     for i = 1:length(Misc.activeDOFs)
%         dofInfo = split(Misc.activeDOFs{i},'/');
%         for dof = 1:auxdata.Ndof
%             if contains(DatStore.DOFNames{dof}, dofInfo{1})
%                 % Max possible torque is smallest ID torque peak
%                 if max(abs(DatStore.T_exp(:,dof))) < auxdata.Tmax_act
%                     auxdata.Tmax_act = max(abs(DatStore.T_exp(:,dof)));
%                 end
%                 
%                 switch dofInfo{1}
%                     case 'hip'
%                         signMoment = sign(auxdata.paramsUpper(auxdata.active.hip));
%                         auxdata.signMoment.hip = signMoment;
%                     case 'knee'
%                         signMoment = auxdata.kneeAngleSign*sign(auxdata.paramsUpper(auxdata.active.knee));
%                         auxdata.signMoment.knee = signMoment;
%                     case 'ankle'
%                         signMoment = sign(auxdata.paramsUpper(auxdata.active.ankle));
%                         auxdata.signMoment.ankle = signMoment;
%                 end
%                 
%             end
%         end
%     end
%        
%     % peak torque
%     peakTorqueGuess = Misc.paramGuess(1) / auxdata.Tmax_act;
%     numExoParams = numExoParams + 1;
%     auxdata.active.params.peak_torque.idx = numExoParams;
%     paramsLower(numExoParams) = peakTorqueGuess*0.75;
%     paramsUpper(numExoParams) = peakTorqueGuess*1.25;
% 
%     % peak time
%     peakTimeGuess = Misc.paramGuess(2);
%     auxdata.T_exo_peak_time = peakTimeGuess;
%     numExoParams = numExoParams + 1;
%     auxdata.active.params.shift_time.idx = numExoParams;
%     paramsLower(numExoParams) = -0.25;
%     paramsUpper(numExoParams) = 0.25;
%     
%     auxdata.same_torque_gain = true;
%     auxdata.numExoParams = numExoParams;
%     auxdata.paramsLower = paramsLower;
%     auxdata.paramsUpper = paramsUpper;
% end

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

% Solve static optimization problem with devices for ExoTopology study to
% improve quality of initial guesses
if strcmp(study{2},'Topology') && ~strcmp(Misc.subcase, 'ActParam') && ~strcmp(Misc.subcase, 'Exp') && ~strcmp(Misc.subcase, 'FitOpt')
    DatStore = SolveStaticOptimization_ExoTopology(DatStore, auxdata, Misc);
end

% Problem bounds 
a_min = 0; a_max = 1;             % bounds on muscle activation
vA_min = -1/100; vA_max = 1/100;  % bounds on derivative of muscle activation
F_min = 0; F_max = 5;             % bounds on normalized tendon force
dF_min = -100; dF_max = 100;      % bounds on derivative of normalzied tendon force

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
vAmin = vA_min./auxdata.tauDeact;
vAmax = vA_max./auxdata.tauAct;
dFMin = dF_min*ones(1,auxdata.NMuscles);
dFMax = dF_max*ones(1,auxdata.NMuscles);
aTmin = -1*ones(1,auxdata.Ndof); 
aTmax = 1*ones(1,auxdata.Ndof);
aDmin = zeros(1, auxdata.numActiveDOFs); 
aDmax = ones(1, auxdata.numActiveDOFs);
if strcmp(study{2},'Topology')
    if auxdata.hasActiveDevice && ~strcmp(Misc.subcase, 'ActParam') && ~strcmp(Misc.subcase, 'Exp') && ~strcmp(Misc.subcase, 'FitOpt')
        control_bounds_lower = [vAmin aTmin dFMin aDmin];
        control_bounds_upper = [vAmax aTmax dFMax aDmax];
    else
        control_bounds_lower = [vAmin aTmin dFMin];
        control_bounds_upper = [vAmax aTmax dFMax];
    end
else
    control_bounds_lower = [vAmin aTmin dFMin];
    control_bounds_upper = [vAmax aTmax dFMax];
end
bounds.phase.control.lower = control_bounds_lower; 
bounds.phase.control.upper = control_bounds_upper;

% States bounds
actMin = a_min*ones(1,auxdata.NMuscles); 
actMax = a_max*ones(1,auxdata.NMuscles);
F0min = F_min*ones(1,auxdata.NMuscles);
F0max = F_max*ones(1,auxdata.NMuscles);
Ffmin = F_min*ones(1,auxdata.NMuscles);
Ffmax = F_max*ones(1,auxdata.NMuscles);
FMin = F_min*ones(1,auxdata.NMuscles);
FMax = F_max*ones(1,auxdata.NMuscles);
bounds.phase.initialstate.lower = [actMin, F0min]; 
bounds.phase.initialstate.upper = [actMax, F0max];
bounds.phase.state.lower = [actMin, FMin]; 
bounds.phase.state.upper = [actMax, FMax];
bounds.phase.finalstate.lower = [actMin, Ffmin]; 
bounds.phase.finalstate.upper = [actMax, Ffmax];
if auxdata.shift_exo_peaks
    bounds.phase.initialstate.lower = [bounds.phase.initialstate.lower a_min]; 
    bounds.phase.initialstate.upper = [bounds.phase.initialstate.upper a_max];
    bounds.phase.state.lower = [bounds.phase.state.lower a_min]; 
    bounds.phase.state.upper = [bounds.phase.state.upper a_max];
    bounds.phase.finalstate.lower = [bounds.phase.finalstate.lower a_min]; 
    bounds.phase.finalstate.upper = [bounds.phase.finalstate.upper a_max];
end

% Integral bounds
bounds.phase.integral.lower = 0;
bounds.phase.integral.upper = 10000*(tf-t0);

% Parameter bounds
if auxdata.shift_exo_peaks
    auxdata.paramsLower = [auxdata.paramsLower 0.1];
    auxdata.paramsUpper = [auxdata.paramsUpper 0.5];
end
if strcmp(study{2},'Topology') && ~strcmp(Misc.subcase, 'FitOpt')
    % Parameter variable in the optimization problem always lie in the range 
    % [-1 1], and should be converted as necessary in the continuous function
    % for to make the correct computations.
    bounds.parameter.lower = -1*ones(size(auxdata.paramsLower));
    bounds.parameter.upper = 1*ones(size(auxdata.paramsUpper));
end


% Path constraints
HillEquil = zeros(1, auxdata.NMuscles);
ID_bounds = zeros(1, auxdata.Ndof);
act1_lower = zeros(1, auxdata.NMuscles);
act1_upper = inf*ones(1, auxdata.NMuscles);
act2_lower = -inf*ones(1, auxdata.NMuscles);
act2_upper = 1*ones(1, auxdata.NMuscles)./auxdata.tauAct;
bounds.phase.path.lower = [ID_bounds,HillEquil,act1_lower,act2_lower];
bounds.phase.path.upper = [ID_bounds,HillEquil,act1_upper,act2_upper];
% if ~isempty(Misc.powerMatchType)
%    bounds.phase.path.lower = [bounds.phase.path.lower, Misc.powerMatchValue*0.8];
%    bounds.phase.path.upper = [bounds.phase.path.upper, Misc.powerMatchValue*1.20];
%    auxdata.powerMatchType = Misc.powerMatchType;
% end

% Eventgroup
% Impose mild periodicity
pera_lower = -1 * ones(1, auxdata.NMuscles);
pera_upper = 1 * ones(1, auxdata.NMuscles);
perFtilde_lower = -1*ones(1,auxdata.NMuscles);
perFtilde_upper = 1*ones(1,auxdata.NMuscles);
bounds.eventgroup.lower = [pera_lower perFtilde_lower]; 
bounds.eventgroup.upper = [pera_upper perFtilde_upper];

% Initial guesses

% Time guess
N = length(DatStore.time);
guess.phase.time = DatStore.time;

% Control guess 
if strcmp(study{2},'Topology') && auxdata.hasActiveDevice && ~strcmp(Misc.subcase, 'ActParam') && ~strcmp(Misc.subcase, 'Exp') && ~strcmp(Misc.subcase, 'FitOpt')
    control_guess = [DatStore.SO_MAct DatStore.SO_RAct 0.01*ones(N,auxdata.NMuscles) DatStore.SO_ExoAct];
else
    control_guess = [DatStore.SoAct DatStore.SoRAct 0.01*ones(N,auxdata.NMuscles)];
end
guess.phase.control = control_guess;

% State guess
guess.phase.state =  [DatStore.SoAct 0.2*ones(N,auxdata.NMuscles)];
if auxdata.shift_exo_peaks
    guess.phase.state = [guess.phase.state 0.2*ones(N,1)];
end

% Integral guess
guess.phase.integral = 0;

% Parameter guess
if strcmp(study{2},'Topology') && ~strcmp(Misc.subcase, 'FitOpt')
    guess.parameter = DatStore.SO_parameter;
end
if auxdata.shift_exo_peaks
    guess.parameter = [guess.parameter 0];
end

% Empty exosuit force and torque data structures
DatStore.T_exo = zeros(length(DatStore.time),auxdata.Ndof);

% TODO: figure out what to do with this case
% if strcmp(study{2},'Topology') && strcmp(Misc.subcase, 'Exp')
%     % Exosuit moment curves
%     currentFile = mfilename('fullpath');
%     [currentDir,~] = fileparts(currentFile);
%     ExoCurves = load(fullfile(currentDir,'Data','Quinlivan2017','ExoCurves.mat'));
%     exoTime = ExoCurves.time;
%     exoNormalizedMoment = ExoCurves.hm_norm;
%     exoNormalizedMomentInterp = interp1(linspace(0,100,length(exoTime)), ...
%                                 exoNormalizedMoment, ...
%                                 linspace(0,100,length(DatStore.time)));
%     
%     for i = 1:length(Misc.activeDOFs)
%         dofInfo = split(Misc.activeDOFs{i},'/');
%         for dof = 1:auxdata.Ndof
%                         
%             if contains(DatStore.DOFNames{dof}, dofInfo{1})
%                 
%                 switch dofInfo{1}
%                     case 'hip'
%                         signMoment = sign(auxdata.paramsUpper(auxdata.active.hip));
%                     case 'knee'
%                         signMoment = auxdata.kneeAngleSign*sign(auxdata.paramsUpper(auxdata.active.knee));
%                     case 'ankle'
%                         signMoment = sign(auxdata.paramsUpper(auxdata.active.ankle));
%                 end
%                 
%                 DatStore.T_exo(:,dof) = exoNormalizedMomentInterp;
%   
%             end              
%         end
%     end
%     
% end

if strcmp(study{2},'Topology') && strcmp(Misc.subcase, 'FitOpt')
    exoTime = Misc.exoTime;
    exoApproximatedMoment = Misc.exoApproximatedMoment;
    
    for i = 1:length(Misc.activeDOFs)
        dofInfo = split(Misc.activeDOFs{i},'/');
        for dof = 1:auxdata.Ndof
            if contains(DatStore.DOFNames{dof}, dofInfo{1})
                switch dofInfo{1}
                    case 'hip'
                        signMoment = sign(auxdata.paramsUpper(auxdata.active.hip));
                    case 'knee'
                        signMoment = auxdata.kneeAngleSign*sign(auxdata.paramsUpper(auxdata.active.knee));
                    case 'ankle'
                        signMoment = sign(auxdata.paramsUpper(auxdata.active.ankle));
                end
                
                DatStore.T_exo(:,dof) = signMoment*interp1(linspace(0,100,length(exoTime)), ...
                    exoApproximatedMoment, ...
                    linspace(0,100,length(DatStore.time)));
            end
        end
    end
end

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
setup.nlp.ipoptoptions.maxiterations = 100000;
setup.derivatives.supplier = 'sparseCD';
setup.derivatives.derivativelevel = 'first';
setup.derivatives.dependencies = 'sparse';
setup.scales.method = 'none';
setup.mesh.method = 'hp-PattersonRao';
setup.mesh.tolerance = 1e-4;
setup.mesh.maxiterations = 20;
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

setup.functions.continuous = str2func(['Wrap4continous_Ftilde_vA' tag]);
setup.adigatorgrd.continuous = str2func(['continous_Ftilde_vA' tag 'GrdWrap']);
setup.adigatorgrd.endpoint   = str2func(['endpoint_Ftilde' tag 'ADiGatorGrd']);
setup.adigatorhes.continuous = str2func(['continous_Ftilde_vA' tag 'HesWrap']);
setup.adigatorhes.endpoint   = str2func(['endpoint_Ftilde' tag 'ADiGatorHes']);


%% ---------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% PART III: SOLVE OPTIMAL CONTROL PROBLEM ------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %

output = gpops2(setup);

MuscleNames = DatStore.MuscleNames;
res = output.result.solution.phase(1);
Time = res.time;
MActivation = res.state(:,1:auxdata.NMuscles);
TForcetilde = res.state(:,auxdata.NMuscles+1:auxdata.NMuscles*2);
TForce = TForcetilde.*(ones(size(Time))*DatStore.Fiso);
vA=100*res.control(:,1:auxdata.NMuscles);
dTForcetilde = 10*res.control(:,auxdata.NMuscles+auxdata.Ndof+1:auxdata.NMuscles+auxdata.Ndof+auxdata.NMuscles); 
MExcitation = computeExcitationRaasch(MActivation, vA, auxdata.tauDeact, auxdata.tauAct);
RActivation = res.control(:,auxdata.NMuscles+1:auxdata.NMuscles+auxdata.Ndof);
OptInfo = output;

% Calculate muscle metabolic rates
mat.Time = Time;
mat.DatStore = DatStore;
mat.OptInfo = OptInfo;
mat.MuscleNames = MuscleNames;

MetabolicRate.whole_body = calcWholeBodyMetabolicRate(model, mat);
muscle_energy_rates = calcIndividualMetabolicRate(model, mat);
if isnan(muscle_energy_rates)
    warning(['ERROR: NaN values returned for muscle metabolic rates, ' ...
             'setting to zeros'])
    muscle_energy_rates = zeros(1,length(MuscleNames));
end
MetabolicRate.individual_muscles = muscle_energy_rates(end,:);
DatStore.MetabolicRate = MetabolicRate;

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
    if strcmp(Misc.subcase, 'Act')
        [DatStore.ExoTorques_Act, DatStore.MomentArms_Act] = ...
            calcExoTorques_Ftilde_vAExoTopology_Act(OptInfo, DatStore);
    elseif strcmp(Misc.subcase, 'Pass')
        [DatStore.ExoTorques_Pass, DatStore.MomentArms_Pass, DatStore.passiveForce, ...
            DatStore.pathLength, DatStore.jointAngles, DatStore.slackLength] = ...
            calcExoTorques_Ftilde_vAExoTopology_Pass(OptInfo, DatStore);
    elseif strcmp(Misc.subcase, 'ActPass')
        [DatStore.ExoTorques_Act, DatStore.MomentArms_Act, DatStore.ExoTorques_Pass, DatStore.MomentArms_Pass, ...
            DatStore.passiveForce, DatStore.pathLength, DatStore.jointAngles, DatStore.slackLength] = ...
            calcExoTorques_Ftilde_vAExoTopology_ActPass(OptInfo, DatStore);
    elseif strcmp(Misc.subcase, 'ActParam')
        [DatStore.ExoTorques_Act] = ...
            calcExoTorques_Ftilde_vAExoTopology_ActParam(OptInfo, DatStore);
    elseif strcmp(Misc.subcase, 'Exp')
        [DatStore.ExoTorques_Act] = ...
            calcExoTorques_Ftilde_vAExoTopology_Exp(OptInfo, DatStore);
    elseif strcmp(Misc.subcase, 'FitOpt')
        [DatStore.ExoTorques_Act] = ...
            calcExoTorques_Ftilde_vAExoTopology_FitOpt(OptInfo, DatStore);
    end
end

end

