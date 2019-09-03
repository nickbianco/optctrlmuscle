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
if length(study) == 1
   study{2} = ''; 
end
switch study{1}
    case 'DeGroote2016'
        tag = '';
    case 'SoftExosuitDesign'
        tag = ['Exo' study{2}];
    case 'ISB2017'
        tag = ['ISB' study{2}];
    case 'AnkleHipExosuit'
        tag = 'AHE';
    otherwise
        tag = '';
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
   Misc.Mesh_Frequency=100; % 100
end
% Which study?
if ~isfield(Misc,'study') || isempty(Misc.study)
   Misc.study='DeGroote2016';
end
% Device force level for Quinlivan study
if ~isfield(Misc,'exo_force_level') || isempty(Misc.exo_force_level)
   Misc.exo_force_level = -1;          
end
% Clutched spring stiffness for Collins study
if ~isfield(Misc, 'ankle_clutched_spring_stiffness') || isempty(Misc.ankle_clutched_spring_stiffness)
    Misc.ankle_clutched_spring_stiffness = -1;
end
% Fixed spring rest length for Collins study
if ~isfield(Misc, 'fixed_rest_length') || isempty(Misc.fixed_rest_length)
    Misc.fixed_rest_length = true; 
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
% Modify individual muscle optimal fiber length (lMo)
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
% Quinlivan: shift exoskeleton peaks to match ID peaks
if ~isfield(Misc, 'shift_exo_peaks') || isempty(Misc.shift_exo_peaks)
   Misc.shift_exo_peaks = false; 
end
if ~isfield(Misc, 'data_path') || isempty(Misc.data_path) 
   Misc.data_path = ''; 
end
if ~isfield(Misc, 'cycle') || isempty(Misc.cycle)
   Misc.cycle = -1; 
end
if ~isfield(Misc, 'subject') || isempty(Misc.subject) 
   Misc.data_path = ''; 
end
if ~isfield(Misc, 'condition') || isempty(Misc.condition) 
   Misc.data_path = ''; 
end
% ----------------------------------------------------------------------- %
% Check that options are being specified correctly -----------------------%
if (strcmp(study{1},'SoftExosuitDesign') && ~strcmp(study{2},'Topology')) ... 
    || strcmp(study{1}, 'AnkleHipExosuit')
    errmsg = [study{1} ': must specify force level'];
    assert(Misc.exo_force_level ~= -1,errmsg)
end

if (~strcmp(study{1},'SoftExosuitDesign') && strcmp(study{2},'Collins2015')) ... 
   || ~strcmp(study{1}, 'AnkleHipExosuit')
   errmsg = [study{1} '/' study{2} ': exosuit force level unused'];
   assert(Misc.exo_force_level == -1,errmsg)
end

if ~strcmp(study{2},'Collins2015')
    errmsg = [study{2} ': spring stiffness level unused'];
    assert(Misc.ankle_clutched_spring_stiffness == -1, errmsg)
end

if ~strcmp(study{2},'Quinlivan2017') && ~strcmp(study{1},'AnkleHipExosuit')
    errmsg = [study{2} ': shifting exoskeleton force peaks unused'];
    assert(Misc.shift_exo_peaks == false, errmsg)
end

if ~strcmp(study{2},'Topology')
    errmsg = [study{2} ': active device DOFs unused'];
    assert(isempty(Misc.activeDOFs), errmsg)
    
    errmsg = [study{2} ': passive device DOFs unused'];
    assert(isempty(Misc.passiveDOFs), errmsg)
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
[DatStore.params,DatStore.lOpt,DatStore.L_TendonSlack,DatStore.Fiso,DatStore.PennationAngle,DatStore.metabolicParams]=ReadMuscleParameters(model_path,DatStore.MuscleNames);

% Modify muscle properties
fprintf('Muscles with modified properties: \n')
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
if ~strcmp(study{2},'Topology')
    DatStore = SolveStaticOptimization_IPOPT(DatStore);
end

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

% Collins et al. 2015 study, optimizing for spring stiffness
if strcmp(study{2}, 'Collins2015')
    
    % Find clutched spring DOFs
    % TODO: support separately clutching left and right leg.
    auxdata.clutched_spring_dofs = strmatch('ankle_angle',DatStore.DOFNames);
      
    % Find fixed rest length based on first ankle angle peak after
    % heel strike
    if Misc.fixed_rest_length
        errmsg = 'Fixing spring stiffness only supported for one DOF assistance';
        assert(length(auxdata.clutched_spring_dofs)==1,errmsg)
        
        % Flip sign to match convention in Collins et al. 2015 study
        ankleAngle = -DatStore.q_exp(:,auxdata.clutched_spring_dofs);
        
        % To set fixed rest length, find first ankle angle peak 
        % within first ~40% of gait cycle
        percentGC = 0.4;
        [restLength,firstPeakIdx] = max(ankleAngle(1:floor(percentGC*length(ankleAngle))));
        auxdata.rest_length_first_peak = DatStore.time(firstPeakIdx);
        
        % Find where the rest length is reached after spring recoil
        [~,minIdx] = min(ankleAngle);
        [~,maxIdx] = max(ankleAngle);
        [~,idx] = min(abs(restLength-ankleAngle(minIdx:maxIdx)));
        restLengthAfterRecoilIdx = idx+minIdx-1;
        auxdata.rest_length_after_recoil = DatStore.time(restLengthAfterRecoilIdx);
        
        % Save rest length (in radians)
        auxdata.rest_length = -restLength * (pi/180);
        fprintf('Spring is active in [%f, %f], with rest length %f.\n', ...
            auxdata.rest_length_first_peak, ...
            auxdata.rest_length_after_recoil, ...
            -restLength);
    end
end

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
    
    smallMomentArm = 0.03;
    largeMomentArm = 0.10;
    
    numExoParams = 0;
    paramsLower = [];
    paramsUpper = [];
    % Active device indicies
    if ~isempty(Misc.activeDOFs)
        auxdata.hasActiveDevice = true;
        auxdata.Fmax_act = 25*auxdata.model_mass; % N/kg * kg
        auxdata.active.hip = 0;
        auxdata.active.knee = 0;
        auxdata.active.ankle = 0;
        for i = 1:length(Misc.activeDOFs)
            dofInfo = split(Misc.activeDOFs{i},'/');
            switch dofInfo{1}
                case 'hip'
                    numExoParams = numExoParams + 1;
                    auxdata.active.hip = numExoParams;
                    paramsLower(numExoParams) = -0.10; %#ok<*AGROW>
                    paramsUpper(numExoParams) = 0.10;
                    if length(dofInfo) > 1
                        switch dofInfo{2}
                            case 'flex'
                                if ~isempty(Misc.fixMomentArms)
                                    paramsLower(numExoParams) = Misc.fixMomentArms;
                                    paramsUpper(numExoParams) = Misc.fixMomentArms;
                                else
                                    paramsLower(numExoParams) = smallMomentArm;
                                    paramsUpper(numExoParams) = largeMomentArm;
                                end
                            case 'ext'
                                if ~isempty(Misc.fixMomentArms)
                                    paramsLower(numExoParams) = -Misc.fixMomentArms;
                                    paramsUpper(numExoParams) = -Misc.fixMomentArms;
                                else
                                    paramsLower(numExoParams) = -largeMomentArm;
                                    paramsUpper(numExoParams) = -smallMomentArm;
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
                                    paramsLower(numExoParams) = smallMomentArm;
                                    paramsUpper(numExoParams) = largeMomentArm;
                                end
                            case 'flex'
                                if ~isempty(Misc.fixMomentArms)
                                    paramsLower(numExoParams) = -Misc.fixMomentArms;
                                    paramsUpper(numExoParams) = -Misc.fixMomentArms;
                                else
                                    paramsLower(numExoParams) = -largeMomentArm;
                                    paramsUpper(numExoParams) = -smallMomentArm;
                                end
                        end
                    end
                case 'ankle'
                    numExoParams = numExoParams + 1;
                    auxdata.active.ankle = numExoParams;
                    paramsLower(numExoParams) = -0.10;
                    paramsUpper(numExoParams) = 0.10;
                    if length(dofInfo) > 1
                        switch dofInfo{2}
                            case 'dorsi'
                                if ~isempty(Misc.fixMomentArms)
                                    paramsLower(numExoParams) = Misc.fixMomentArms;
                                    paramsUpper(numExoParams) = Misc.fixMomentArms;
                                else
                                    paramsLower(numExoParams) = smallMomentArm;
                                    paramsUpper(numExoParams) = largeMomentArm;
                                end
                            case 'plantar'
                                if ~isempty(Misc.fixMomentArms)
                                    paramsLower(numExoParams) = -Misc.fixMomentArms;
                                    paramsUpper(numExoParams) = -Misc.fixMomentArms;
                                else
                                    paramsLower(numExoParams) = -largeMomentArm;
                                    paramsUpper(numExoParams) = -smallMomentArm;
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
    if ~isempty(Misc.passiveDOFs)
        auxdata.hasPassiveDevice = true;
        auxdata.passiveStiffness = 1250*auxdata.model_mass; % k ~= 100 kN/m for 80kg subject (van den Bogert 2003) 
        auxdata.passive.hip = 0;
        auxdata.passive.knee = 0;
        auxdata.passive.ankle = 0;
        for i = 1:length(Misc.passiveDOFs)
            dofInfo = split(Misc.passiveDOFs{i},'/');
            switch dofInfo{1}
                case 'hip'
                    numExoParams = numExoParams + 1;
                    auxdata.passive.hip = numExoParams;
                    paramsLower(numExoParams) = -0.10;
                    paramsUpper(numExoParams) = 0.10;
                    if length(dofInfo) > 1
                        switch dofInfo{2}
                            case 'flex'
                                if ~isempty(Misc.fixMomentArms)
                                    paramsLower(numExoParams) = Misc.fixMomentArms;
                                    paramsUpper(numExoParams) = Misc.fixMomentArms;
                                else
                                    paramsLower(numExoParams) = smallMomentArm;
                                    paramsUpper(numExoParams) = largeMomentArm;
                                end
                            case 'ext'
                                if ~isempty(Misc.fixMomentArms)
                                    paramsLower(numExoParams) = -Misc.fixMomentArms;
                                    paramsUpper(numExoParams) = -Misc.fixMomentArms;
                                else
                                    paramsLower(numExoParams) = -largeMomentArm;
                                    paramsUpper(numExoParams) = -smallMomentArm;
                                end
                        end
                    end
                case 'knee'
                    numExoParams = numExoParams + 1;
                    auxdata.passive.knee = numExoParams;
                    paramsLower(numExoParams) = -0.10;
                    paramsUpper(numExoParams) = 0.10;
                    if length(dofInfo) > 1
                        switch dofInfo{2}
                            case 'ext'
                                if ~isempty(Misc.fixMomentArms)
                                    paramsLower(numExoParams) = Misc.fixMomentArms;
                                    paramsUpper(numExoParams) = Misc.fixMomentArms;
                                else
                                    paramsLower(numExoParams) = smallMomentArm;
                                    paramsUpper(numExoParams) = largeMomentArm;
                                end
                            case 'flex'
                                if ~isempty(Misc.fixMomentArms)
                                    paramsLower(numExoParams) = -Misc.fixMomentArms;
                                    paramsUpper(numExoParams) = -Misc.fixMomentArms;
                                else
                                    paramsLower(numExoParams) = -largeMomentArm;
                                    paramsUpper(numExoParams) = -smallMomentArm;
                                end
                        end
                    end
                    
                case 'ankle'
                    numExoParams = numExoParams + 1;
                    auxdata.passive.ankle = numExoParams;
                    paramsLower(numExoParams) = -0.10;
                    paramsUpper(numExoParams) = 0.10;
                    if length(dofInfo) > 1
                        switch dofInfo{2}
                            case 'dorsi'
                                if ~isempty(Misc.fixMomentArms)
                                    paramsLower(numExoParams) = Misc.fixMomentArms;
                                    paramsUpper(numExoParams) = Misc.fixMomentArms;
                                else
                                    paramsLower(numExoParams) = smallMomentArm;
                                    paramsUpper(numExoParams) = largeMomentArm;
                                end
                            case 'plantar'
                                if ~isempty(Misc.fixMomentArms)
                                    paramsLower(numExoParams) = -Misc.fixMomentArms;
                                    paramsUpper(numExoParams) = -Misc.fixMomentArms;
                                else
                                    paramsLower(numExoParams) = -largeMomentArm;
                                    paramsUpper(numExoParams) = -smallMomentArm;
                                end
                        end
                    end
            end
        end
        % Extra index for exotendon slack length
        numExoParams = numExoParams + 1;
        paramsLower(numExoParams) = 0.75;
        paramsUpper(numExoParams) = 1.25;
        auxdata.passive.slack_length = numExoParams;
    end
    auxdata.numExoParams = numExoParams;
    auxdata.paramsLower = paramsLower;
    auxdata.paramsUpper = paramsUpper;
    
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

% Pre-calibration to passive moments from Silder et al. 2007
% need to handle if we want to calib or not, and if we have done it already
% conditional for if empty struct -> we don't want to calib
% conditional for if we have the file already -> don't recalibrate

% hardcode the path to the subject we are looking at
% TODO: figure out how to go to each subject

% fprintf(Misc.trail
calibrationparampath = ['C:\Users\JP\code\repos\Stanford\delplab\projects\AnkleHipExosuit\results\experiments\',Misc.trialName(17:26),'\calibratedModifiers.mat'];

if isempty(fieldnames(Misc.parameterCalibrationTerms))
    fprintf('\n***********************************************************')
    fprintf('\nSkipping passive moment calibration: no parameters chosen.\n')
    if isfile(calibrationparampath)
        fprintf('But we have previous calibration... loading it now.\n')
        load(calibrationparampath);
    else
        fprintf('No previous calibration, using default scaled model.\n')
    end
else
    % check here for the existing calibration
    if isfile(calibrationparampath)
        fprintf('\nWe have a previous passive calibration, loading it now.\n')
        load(calibrationparampath);                
    else
        fprintf('\nBeginning passive moment calibration.\n')
        calibratedModifiers = PassiveMomentCalibration(model_path, auxdata, DatStore(1), Misc.parameterCalibrationTerms);
        save(calibrationparampath,'calibratedModifiers');
        
    end
    % Update parameter structs
    indices = calibratedModifiers.indices;
    auxdata.params(9,indices.lMo) = calibratedModifiers.lMo;
    auxdata.params(10,indices.lTs) = calibratedModifiers.lTs;
    auxdata.params(7,indices.e0) = calibratedModifiers.e0;
end




% Parameters of passive muscle force-length characteristic, and tendon
% characteristic
e0 = 0.6*DatStore.params(7,:); 
kpe = 4*DatStore.params(8,:); 
t50 = exp(kpe .* (0.2 - 0.10e1) ./ e0);
pp1 = (t50 - 0.10e1); t7 = exp(kpe); pp2 = (t7 - 0.10e1);
auxdata.Fpparam = [pp1;pp2;ones(1,length(pp1))*Misc.tendonStiffnessCoeff];

% Solve static optimization problem with devices for ExoTopology study to
% improve quality of initial guesses
if strcmp(study{2},'Topology')
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
% !!! fix the control bounds for the excitations if doing explicit!!!
vAmin = vA_min./auxdata.tauDeact;
vAmax = vA_max./auxdata.tauAct;
dFMin = dF_min*ones(1,auxdata.NMuscles);
dFMax = dF_max*ones(1,auxdata.NMuscles);
aTmin = -1*ones(1,auxdata.Ndof); 
aTmax = 1*ones(1,auxdata.Ndof);
aDmin = zeros(1, auxdata.numActiveDOFs); 
aDmax = ones(1, auxdata.numActiveDOFs);
switch study{2}
    case {'HipAnkle','HipKneeAnkle','HipExtHipAbd'}
        control_bounds_lower = [vAmin aTmin dFMin aDmin];
        control_bounds_upper = [vAmax aTmax dFMax aDmax];
    case 'Topology'
        if auxdata.hasActiveDevice
            control_bounds_lower = [vAmin aTmin dFMin aDmin];
            control_bounds_upper = [vAmax aTmax dFMax aDmax];
        else
            control_bounds_lower = [vAmin aTmin dFMin];
            control_bounds_upper = [vAmax aTmax dFMax];
        end
    otherwise
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

% Integral bounds
bounds.phase.integral.lower = 0;
bounds.phase.integral.upper = 10000*(tf-t0);
% keyboard
% Parameter bounds
switch study{2}
    case 'HipAnkle'
        bounds.parameter.lower = -1.0;
        bounds.parameter.upper = 1.0;
    case 'HipKneeAnkle'
        bounds.parameter.lower = -1.0;
        bounds.parameter.upper = 1.0;
    case 'HipExtHipAbd'
        bounds.parameter.lower = -1.0;
        bounds.parameter.upper = 1.0;
    case 'HipAnkleMass'
        bounds.parameter.lower = -1.0;
        bounds.parameter.upper = 1.0;
    case 'Collins2015'
        stiffness_lower = 0;
        stiffness_upper = 1;
        if Misc.ankle_clutched_spring_stiffness ~= -1
            assert(Misc.ankle_clutched_spring_stiffness >= 0 && ...
                   Misc.ankle_clutched_spring_stiffness <= 1);
            stiffness_lower = Misc.ankle_clutched_spring_stiffness;
            stiffness_upper = Misc.ankle_clutched_spring_stiffness;
        end
        rest_length_lower = -0.5;
        rest_length_upper = 0.5;
        if Misc.fixed_rest_length
            rest_length_lower = auxdata.rest_length;
            rest_length_upper = auxdata.rest_length;
        end
        bounds.parameter.lower = [stiffness_lower, rest_length_lower];
        bounds.parameter.upper = [stiffness_upper, rest_length_upper];
    case 'Quinlivan2017'
        force_level_lower = 0;
        force_level_upper = 10;
        if Misc.exo_force_level ~= -1
            assert(Misc.exo_force_level >= 0 && ...
                   Misc.exo_force_level <= 10);
            force_level_lower = Misc.exo_force_level;
            force_level_upper = Misc.exo_force_level;
        end
        bounds.parameter.lower = force_level_lower;
        bounds.parameter.upper = force_level_upper;
    case 'Topology'
        bounds.parameter.lower = auxdata.paramsLower;
        bounds.parameter.upper = auxdata.paramsUpper;
end

% Path constraints
HillEquil = zeros(1, auxdata.NMuscles);
ID_bounds = zeros(1, auxdata.Ndof);
act1_lower = zeros(1, auxdata.NMuscles);
act1_upper = inf*ones(1, auxdata.NMuscles);
act2_lower = -inf*ones(1, auxdata.NMuscles);
act2_upper = ones(1, auxdata.NMuscles)./auxdata.tauAct;
bounds.phase.path.lower = [ID_bounds,HillEquil ,act1_lower,act2_lower];
bounds.phase.path.upper = [ID_bounds,HillEquil ,act1_upper,act2_upper];
% extras in the constraints for implicit dynamics formulation

% Eventgroup
% Impose mild periodicity
pera_lower = -1 * ones(1, auxdata.NMuscles);
pera_upper = 1 * ones(1, auxdata.NMuscles);
perFtilde_lower = -1*ones(1,auxdata.NMuscles);
perFtilde_upper = 1*ones(1,auxdata.NMuscles);
bounds.eventgroup.lower = [pera_lower perFtilde_lower]; 
bounds.eventgroup.upper = [pera_upper perFtilde_upper];

% Initial guess
N = length(DatStore.time);
guess.phase.time = DatStore.time;
switch study{2}
    case {'HipAnkle','HipKneeAnkle','HipExtHipAbd'}
        control_guess = [zeros(N,auxdata.NMuscles) zeros(N,auxdata.Ndof) 0.01*ones(N,auxdata.NMuscles) 0.5*ones(N,1)];
    case 'Topology'
        if auxdata.hasActiveDevice
            control_guess = [zeros(N,auxdata.NMuscles) DatStore.SO_RAct 0.01*ones(N,auxdata.NMuscles) DatStore.SO_ExoAct];
        else
            control_guess = [zeros(N,auxdata.NMuscles) DatStore.SO_RAct 0.01*ones(N,auxdata.NMuscles)];
        end
    otherwise
        control_guess = [zeros(N,auxdata.NMuscles) zeros(N,auxdata.Ndof) 0.01*ones(N,auxdata.NMuscles)];
end
guess.phase.control = control_guess;

switch study{2}
    case 'Topology'
        guess.phase.state =  [DatStore.SO_MAct 0.2*ones(N,auxdata.NMuscles)];
    otherwise
        guess.phase.state =  [0.2*ones(N,auxdata.NMuscles) 0.2*ones(N,auxdata.NMuscles)];
end
guess.phase.integral = 0;
switch study{2}
    case 'HipAnkle'
        guess.parameter = 0;
    case 'HipKneeAnkle'
        guess.parameter = 0;
    case 'HipExtHipAbd'
        guess.parameter = 0;
    case 'HipAnkleMass'
        guess.parameter = 0;
    case 'Collins2015'
        guess.parameter = [0.5, 0];
    case 'Quinlivan2017'
        guess.parameter = 4;
    case 'Topology'
       if auxdata.hasPassiveDevice
%            guess.parameter = [zeros(1,auxdata.numExoParams-1) 1];
           guess.parameter = DatStore.SO_parameter;
       else
%            guess.parameter = zeros(1,auxdata.numExoParams);
           guess.parameter = DatStore.SO_parameter;
       end
end

% Empty exosuit force and torque data structures
DatStore.T_exo = zeros(length(DatStore.time),auxdata.Ndof);
T_exo_shift = DatStore.T_exo; % Temp variable for shifting data
DatStore.p_linreg = zeros(2,auxdata.Ndof);

% Reproduce Quinlivan et al. 2017 study
if strcmp(study{1}, 'AnkleHipExosuit')
    
    fname = [Misc.subject '.mat'];
    mat_path = fullfile(Misc.data_path, 'segmented', Misc.subject, fname);
    data = load(mat_path);
    
    subject_mass = data.subject_data.bdy_mas.x.raw;
    cycle = Misc.cycle;
    switch Misc.exo_force_level
        case 1
            exoTime = data.stp_025.non_tme_r.x.raw(:,cycle);
            ankleTorque = data.stp_025.aca_tor_r.x.raw(:,cycle) * subject_mass;
            hipTorque = data.stp_025.ahf_tor_l.x.raw(:,cycle) * subject_mass;
        case 2
            exoTime = data.stp_050.non_tme_r.x.raw(:,cycle);
            ankleTorque = data.stp_050.aca_tor_r.x.raw(:,cycle) * subject_mass;
            hipTorque = data.stp_050.ahf_tor_l.x.raw(:,cycle) * subject_mass;
        case 3
            exoTime = data.stp_075.non_tme_r.x.raw(:,cycle);
            ankleTorque = data.stp_075.aca_tor_r.x.raw(:,cycle) * subject_mass;
            hipTorque = data.stp_075.ahf_tor_l.x.raw(:,cycle) * subject_mass;
            
        case 4
            exoTime = data.stp_100.non_tme_r.x.raw(:,cycle);
            ankleTorque = data.stp_100.aca_tor_r.x.raw(:,cycle) * subject_mass;
            hipTorque = data.stp_100.ahf_tor_l.x.raw(:,cycle) * subject_mass;
        otherwise
            error(['AnkleHipExosuit: exo_force_level must be in the range ' ...
                   '1-4, but the value ' + num2str(Misc.exo_force_level) + ...
                   'was provided.']);
    end
    
    if strcmp(Misc.condition, 'slack')
       exoTime = linspace(time(1), time(2), length(exoTime));
    end
    
    for dof = 1:auxdata.Ndof
        if strfind(DatStore.DOFNames{dof},'ankle_angle')
            % Negative to match ankle_angle_r coord convention
            DatStore.T_exo(:,dof) = -interp1(exoTime, ankleTorque, DatStore.time, 'pchip');
            if Misc.shift_exo_peaks
                [~,idxExo] = max(abs(DatStore.T_exo(:,dof)));
                if DatStore.T_exo(idxExo,dof) > 0
                    [~,idxID] = max(DatStore.T_exp(:,dof));
                else
                    [~,idxID] = min(DatStore.T_exp(:,dof));
                end
                
                shift = idxExo - idxID;
                if shift > 0
                    DatStore.T_exo(:,dof) = [DatStore.T_exo((shift+1):end, dof); ...
                                             DatStore.T_exo(1:shift, dof)];
                else
                    DatStore.T_exo(:,dof) = [DatStore.T_exo((end+shift+1):end, dof); ...
                                             DatStore.T_exo(1:(end+shift), dof)];
                end
            end
        elseif strfind(DatStore.DOFNames{dof},'hip_flexion')
            % Negative to match hip_flexion_r coord convention
            DatStore.T_exo(:,dof) = -interp1(exoTime, hipTorque, DatStore.time, 'pchip');
            if Misc.shift_exo_peaks
                [~,idxExo] = max(abs(DatStore.T_exo(:,dof)));
                if DatStore.T_exo(idxExo,dof) > 0
                    [~,idxID] = max(DatStore.T_exp(:,dof));
                else
                    [~,idxID] = min(DatStore.T_exp(:,dof));
                end
                
                shift = idxExo - idxID;
                if shift > 0
                    DatStore.T_exo(:,dof) = [DatStore.T_exo((shift+1):end,dof); ...
                                             DatStore.T_exo(1:shift, dof)];
                else
                    DatStore.T_exo(:,dof) = [DatStore.T_exo((end+shift+1):end,dof); ...
                                             DatStore.T_exo(1:(end+shift), dof)];
                end
                
            end
            
        end
    end
        
end

if strcmp(study{2},'Quinlivan2017') || strcmp(study{2},'Q2017')
    % Exosuit moment curves
    currentFile = mfilename('fullpath');
    [currentDir,~] = fileparts(currentFile);
    ExoCurves = load(fullfile(currentDir,'Data','Quinlivan2017','ExoCurves.mat'));
    exoTime = ExoCurves.time;
    % Peaks are body mass normalized so multiply by model mass
    exoAnkleMomentPeaks = ExoCurves.am_peak * model_mass;
    exoAnkleNormalizedMoment = ExoCurves.am_norm;
    exoHipMomentPeaks = ExoCurves.hm_peak * model_mass;
    exoHipNormalizedMoment = ExoCurves.hm_norm;
    
    % Interpolate exosuit moments to match data
    switch study{1}
        case 'SoftExosuitDesign'
            if Misc.exo_force_level
                exoAnkleMoment = exoAnkleMomentPeaks(Misc.exo_force_level) * exoAnkleNormalizedMoment;
                exoHipMoment = exoHipMomentPeaks(Misc.exo_force_level) * exoHipNormalizedMoment;
                for dof = 1:auxdata.Ndof
                    if strfind(DatStore.DOFNames{dof},'ankle_angle')
                        % Negative to match ankle_angle_r coord convention
                        DatStore.T_exo(:,dof) = -interp1(exoTime, exoAnkleMoment, DatStore.time);
                    elseif strfind(DatStore.DOFNames{dof},'hip_flexion')
                        % Positive to match hip_flexion_r coord convention
                        DatStore.T_exo(:,dof) = interp1(exoTime, exoHipMoment, DatStore.time);
                    end
                end
            end
        case 'ISB2017'
            % Normalized curves at hip and ankle
            % Linear regression on moment peaks for parameter optimization
            %
            % NOTE: interpolation scheme assumes that time range provided
            %       in problem represents a full gait cycle
            % TODO: find a way make this generic
            if Misc.exo_force_level
                for dof = 1:auxdata.Ndof
                    if strfind(DatStore.DOFNames{dof},'ankle_angle')
                        % Negative to match ankle_angle_r coord convention
                        DatStore.T_exo(:,dof) = -interp1(linspace(0,100,length(exoTime)), ...
                            exoAnkleNormalizedMoment, ...
                            linspace(0,100,length(DatStore.time)));
                        
                        if Misc.shift_exo_peaks
                            [~,exo_idx] = max(-DatStore.T_exo(:,dof));
                            [~,exp_idx] = max(-DatStore.T_exp(:,dof));
                            shift_idx = abs(exp_idx - exo_idx);
                            T_exo_shift(:,dof) = [zeros(shift_idx,1); DatStore.T_exo(1:end-shift_idx,dof)];
                            
                            DatStore.T_exo(:,dof) = T_exo_shift(:,dof);
                        end
                        
                        X = 1:4;
                        Y = exoAnkleMomentPeaks(1:4);
                        DatStore.p_linreg(:,dof) = polyfit(X,Y,1)';
                        % Output peak assistive force at ankle in %BW to use
                        % for abscissa in plots.
                        disp('ISB2017/Quinlivan2017 force levels (%BW):');
                        fprintf('slope: %f; intercept: %f\n', ...
                            DatStore.p_linreg(1,dof) / model_mass, ...
                            DatStore.p_linreg(2,dof) / model_mass);
                        for fl = 1:10
                            fprintf('force level %i: %f\n', fl, ...
                                (DatStore.p_linreg(1,dof)*fl + ...
                                DatStore.p_linreg(2,dof)) / model_mass);
                        end
                        for mp = 1:length(ExoCurves.am_peak)
                            fprintf('exo ankle moment peaks mass-norm %i: %f\n', ...
                                    mp, ExoCurves.am_peak(mp));
                        end
                    elseif strfind(DatStore.DOFNames{dof},'hip_flexion')
                        % Positive to match hip_flexion_r coord convention
                        DatStore.T_exo(:,dof) = interp1(linspace(0,100,length(exoTime)), ...
                            exoHipNormalizedMoment, ...
                            linspace(0,100,length(DatStore.time)));
                        
                        if Misc.shift_exo_peaks 
                            [~,exo_idx] = max(DatStore.T_exo(:,dof));
                            [~,exp_idx] = max(DatStore.T_exp(:,dof));
                            shift_idx = abs(exp_idx - exo_idx);
                            T_exo_shift(:,dof) = [zeros(shift_idx,1); DatStore.T_exo(1:end-shift_idx,dof)];
                            
                            DatStore.T_exo(:,dof) = T_exo_shift(:,dof);
                        end
                        
                        X = 1:4;
                        Y = exoHipMomentPeaks(1:4);
                        DatStore.p_linreg(:,dof) = polyfit(X,Y,1)';
                    end
                end
            end
            auxdata.p_linreg = DatStore.p_linreg;
    end   
end

% "Tradeoff" struct
%  +/- 1 for the two joints you are investigating
%  zeros for everything else
DatStore.tradeoff = zeros(auxdata.Ndof,1);

% Given one actuator, compare tradeoff between hip flexion and ankle 
% plantarflexion assistance
DatStore.Fopt_exo = zeros(auxdata.Ndof,1);
if strcmp(study{2},'HipAnkle') 
    % Exosuit moment curves
    currentFile = mfilename('fullpath');
    [currentDir,~] = fileparts(currentFile);
    ExoCurves = load(fullfile(currentDir,'Data','Quinlivan2017','ExoCurves.mat'));
    % Peaks are body mass normalized so multiply by model mass
    exoAnkleForcePeaks = ExoCurves.af_peak * model_mass;

    % Interpolate exosuit moments to match data
    if Misc.exo_force_level    
        exoForce = exoAnkleForcePeaks(Misc.exo_force_level);
        for dof = 1:auxdata.Ndof
            if strfind(DatStore.DOFNames{dof}, 'ankle_angle')
                % Negative to match ankle_angle_r coord convention
                DatStore.Fopt_exo(dof) = -exoForce;
                DatStore.tradeoff(dof) = -1;
            elseif strfind(DatStore.DOFNames{dof}, 'hip_flexion')
                % Positive to match hip_flexion_r coord convention
                DatStore.Fopt_exo(dof) = exoForce;
                DatStore.tradeoff(dof) = 1;
            end
        end
    end
    
    auxdata.Fopt_exo = DatStore.Fopt_exo;
    auxdata.tradeoff = DatStore.tradeoff;
end

% Given one actuator, compare tradeoff between hip flexion and ankle 
% plantarflexion assistance. Also allow moment arm at knee.
DatStore.Fopt_exo_knee = zeros(auxdata.Ndof,1);
if strcmp(study{2},'HipKneeAnkle') 
    % Exosuit moment curves
    currentFile = mfilename('fullpath');
    [currentDir,~] = fileparts(currentFile);
    ExoCurves = load(fullfile(currentDir,'Data','Quinlivan2017','ExoCurves.mat'));
    % Peaks are body mass normalized so multiply by model mass
    exoAnkleForcePeaks = ExoCurves.af_peak * model_mass;

    % Interpolate exosuit moments to match data
    if Misc.exo_force_level    
        exoForce = exoAnkleForcePeaks(Misc.exo_force_level);
        for dof = 1:auxdata.Ndof
            if strfind(DatStore.DOFNames{dof},'ankle_angle')
                % Negative to match ankle_angle_r coord convention
                DatStore.Fopt_exo(dof) = -exoForce;
                DatStore.tradeoff(dof) = -1;
            elseif strfind(DatStore.DOFNames{dof},'hip_flexion')
                % Positive to match hip_flexion_r coord convention
                DatStore.Fopt_exo(dof) = exoForce;
                DatStore.tradeoff(dof) = 1;
            elseif strfind(DatStore.DOFNames{dof},'knee_angle')
                DatStore.Fopt_exo_knee(dof) = exoForce;
            end
        end
    end
    
    auxdata.Fopt_exo = DatStore.Fopt_exo;
    auxdata.Fopt_exo_knee = DatStore.Fopt_exo_knee;
    auxdata.tradeoff = DatStore.tradeoff;
end

% Given one actuator, compare tradeoff between hip extension and hip 
% abduction assistance
if strcmp(study{2},'HipExtHipAbd') 
    % Exosuit moment curves
    currentFile = mfilename('fullpath');
    [currentDir,~] = fileparts(currentFile);
    ExoCurves = load(fullfile(currentDir,'Data','Quinlivan2017','ExoCurves.mat'));
    % Peaks are body mass normalized so multiply by model mass
    exoAnkleForcePeaks = ExoCurves.af_peak * model_mass;

    % Interpolate exosuit moments to match data
    if Misc.exo_force_level    
        exoForce = exoAnkleForcePeaks(Misc.exo_force_level);
        for dof = 1:auxdata.Ndof
            if strfind(DatStore.DOFNames{dof},'hip_flexion')
                % Negative to match hip_flexion_r coord convention
                DatStore.Fopt_exo(dof) = -exoForce;
                DatStore.tradeoff(dof) = 1;
            elseif strfind(DatStore.DOFNames{dof},'hip_adduction')
                % Negative to match hip_adduction_r coord convention
                DatStore.Fopt_exo(dof) = -exoForce;
                DatStore.tradeoff(dof) = -1;
            end
        end
    end
    
    auxdata.Fopt_exo = DatStore.Fopt_exo;
    auxdata.tradeoff = DatStore.tradeoff;
end

% Given one actuator, compare tradeoff between hip flexion and ankle 
%  plantarflexion assistance with addition of device mass
if strcmp(study{2},'HipAnkleMass') 
    
    % Exosuit moment curves
    currentFile = mfilename('fullpath');
    [currentDir,~] = fileparts(currentFile);
    ExoCurves = load(fullfile(currentDir,'Data','Quinlivan2017','ExoCurves.mat'));
    exoTime = ExoCurves.time;
    % Peaks are body mass normalized so multiply by model mass
    exoAnkleMomentPeaks = ExoCurves.am_peak * model_mass;
    exoAnkleNormalizedMoment = ExoCurves.am_norm;
    exoHipMomentPeaks = ExoCurves.hm_peak * model_mass;
    exoHipNormalizedMoment = ExoCurves.hm_norm;
    
    % Interpolate exosuit moments to match data
    if Misc.exo_force_level    
        exoAnkleMoment = exoAnkleMomentPeaks(Misc.exo_force_level) * exoAnkleNormalizedMoment;
        exoHipMoment = exoHipMomentPeaks(Misc.exo_force_level) * exoHipNormalizedMoment;
        for dof = 1:auxdata.Ndofdoit
            if strfind(DatStore.DOFNames{dof},'ankle_angle')
                % Negative to match ankle_angle_r coord convention
                DatStore.T_exo(:,dof) = -interp1(exoTime, exoAnkleMoment, DatStore.time, 'previous','extrap'); 
                DatStore.tradeoff(dof) = -1;
            elseif strfind(DatStore.DOFNames{dof},'hip_flexion')
                % Positive to match hip_flexion_r coord convention
                DatStore.T_exo(:,dof) = interp1(exoTime, exoHipMoment, DatStore.time, 'previous','extrap');
                DatStore.tradeoff(dof) = 1;
            end
        end
    end  
    auxdata.tradeoff = DatStore.tradeoff;
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
setup.derivatives.derivativelevel = 'first'; % first / second
setup.nlp.ipoptoptions.tolerance = 10^(-4);
setup.nlp.ipoptoptions.maxiterations = 10000;
setup.derivatives.supplier = 'sparseCD'; % sparseCD / adigator
setup.scales.method = 'none';
setup.mesh.method = 'hp-PattersonRao';
setup.mesh.tolerance = 1e-3;
setup.mesh.maxiterations = 2;
setup.mesh.colpointsmin = 3;  % match the setup.mesh.phase.colpoints leading coeff. 
setup.mesh.colpointsmax = 20;
setup.method = 'RPM-integration';
setup.displaylevel = 2;
NMeshIntervals = round((tf-t0)*Misc.Mesh_Frequency);
setup.mesh.phase.colpoints = 3*ones(1,NMeshIntervals);
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
    
    lockDir = 'C:\Users\JP\tmp\adigatorLock\';
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
    if strcmp(setup.derivatives.supplier, 'adigator')
        save(pathLock, 'emptyVar')
    end
    
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
    if strcmp(setup.derivatives.supplier, 'adigator')
        system(sprintf('del %s',pathLock))
    end

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

% Interpolation lMT
lMTinterp = interp1(DatStore.time,DatStore.LMT,Time);
for m = 1:auxdata.NMuscles
    LMTSpline(m) = spline(Time,lMTinterp(:,m));
    [LMT(:,m),VMT(:,m),~] = SplineEval_ppuval(LMTSpline(m),Time,1);
end

MuscleData = DeGroote2016Muscle_FtildeState(MActivation, TForcetilde, ...
    dTForcetilde, LMT, VMT, auxdata.params, auxdata.Fvparam, auxdata.Fpparam, ... 
    auxdata.Faparam);

if strcmp(study{1},'ISB2017')
    if strcmp(study{2},'Quinlivan2017')
        DatStore.ExoTorques_Act = calcExoTorques_FtildeISBQuinlivan2017(...
            OptInfo, DatStore);
    elseif strcmp(study{2},'Collins2015')
        DatStore.ExoTorques_Act = calcExoTorques_FtildeISBCollins2015(...
            OptInfo, DatStore);
    end
end

if strcmp(study{2},'Topology')
    if auxdata.hasActiveDevice && ~auxdata.hasPassiveDevice
        [DatStore.ExoTorques_Act, DatStore.MomentArms_Act] = ...
            calcExoTorques_Ftilde_vAExoTopology_Act(OptInfo, DatStore);
    elseif ~auxdata.hasActiveDevice && auxdata.hasPassiveDevice
        [DatStore.ExoTorques_Pass, DatStore.MomentArms_Pass, DatStore.passiveForce, ...
            DatStore.pathLength, DatStore.jointAngles, DatStore.slackLength] = ...
            calcExoTorques_Ftilde_vAExoTopology_Pass(OptInfo, DatStore);
    elseif auxdata.hasActiveDevice && auxdata.hasPassiveDevice
        [DatStore.ExoTorques_Act, DatStore.MomentArms_Act, DatStore.ExoTorques_Pass, DatStore.MomentArms_Pass, ...
            DatStore.passiveForce, DatStore.pathLength, DatStore.jointAngles, DatStore.slackLength] = ...
            calcExoTorques_Ftilde_vAExoTopology_ActPass(OptInfo, DatStore);
    end
end

if strcmp(study{1}, 'AnkleHipExosuit')
    DatStore.ExoTorques_Act = calcExoTorques_FtildeAHE(OptInfo, DatStore);
end

