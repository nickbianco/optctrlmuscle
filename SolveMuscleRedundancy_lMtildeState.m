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
% This function uses the normalized muscle fiber length as a state (see
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

function [Time,MExcitation,MActivation,RActivation,TForcetilde,TForce,lMtilde,lM,MuscleNames,OptInfo,DatStore]=SolveMuscleRedundancy_lMtildeState(model_path,IK_path,ID_path,time,OutPath,Misc)

%% ---------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% PART I: INPUTS FOR OPTIMAL CONTROL PROBLEM ---------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %

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
% Cost Function
if ~isfield(Misc,'costfun') || isempty(Misc.costfun)
   Misc.costfun='Exc';
end
% Which study?
if ~isfield(Misc,'study') || isempty(Misc.study)
   Misc.study='DeGroote2016';
end
% Device force level for Quinlivan study
if ~isfield(Misc,'exo_force_level') || isempty(Misc.exo_force_level)
   Misc.exo_force_level = 0;
end
% Model mass
if ~isfield(Misc,'model_mass') || isempty(Misc.model_mass)
   Misc.model_mass = 1;
end
% Subcase
if ~isfield(Misc,'subcase') || isempty(Misc.subcase)
    folder = Misc.costfun;
else
    folder = [Misc.costfun '_' Misc.subcase];
end
% variable tendon stiffness
if ~isfield(Misc, 'tendonStiffnessCoeff') || isempty(Misc.tendonStiffnessCoeff)
    Misc.tendonStiffnessCoeff = 35;
end
% Clutched spring for ankle plantarflexion (Collins 2015)
if ~isfield(Misc, 'ankle_clutched_spring') || isempty(Misc.ankle_clutched_spring)
    Misc.ankle_clutched_spring = false;
end
if ~isfield(Misc, 'ankle_clutched_spring_stiffness') || isempty(Misc.ankle_clutched_spring_stiffness)
    Misc.ankle_clutched_spring_stiffness = -1;
end
if ~isfield(Misc, 'phase_boundary') || isempty(Misc.phase_boundary)
    Misc.phase_boundary = NaN;
end

if Misc.ankle_clutched_spring 
    if ~strcmp(Misc.costfun, 'Exc_Act')
        error('ankle_clutched_spring == true requires costfun == ''Exc_Act''');
    end
    Misc.costfun = 'Exc_ActSpr';
    if Misc.ankle_clutched_spring_stiffness ~= -1
        assert(Misc.ankle_clutched_spring_stiffness >= 0 && ...
            Misc.ankle_clutched_spring_stiffness <= 1);
    end
else
    assert(Misc.ankle_clutched_spring_stiffness == -1);
end

% Based on study and cost function, decide which continuous and endpoint
% functions to use
study = strsplit(Misc.study,'/');
switch study{1}
    case 'DeGroote2016'
        tag = '';
    case 'SoftExosuitDesign'
        tag = ['Exo' study{2}];
end
tag = [tag '_' Misc.costfun];

% ------------------------------------------------------------------------%
% Compute ID -------------------------------------------------------------%
if isempty(ID_path) || ~exist(ID_path,'file')
    disp('ID path was not specified or the file does not exist, computation ID started');
    if ~isfield(Misc,'Loads_path') || isempty(Misc.Loads_path) || ~exist(Misc.Loads_path,'file');
        error('External loads file was not specified or does not exist, please add the path to the external loads file: Misc.Loads_path');
    else
        %check the output path for the ID results
        if isfield(Misc,'ID_ResultsPath');
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
[DatStore.params,DatStore.lOpt,DatStore.L_TendonSlack,DatStore.Fiso,DatStore.PennationAngle]=ReadMuscleParameters(model_path,DatStore.MuscleNames);
% Static optimization using IPOPT solver
DatStore = SolveStaticOptimization_IPOPT(DatStore);

%% ---------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% PART II: OPTIMAL CONTROL PROBLEM FORMULATION -------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %

% Input arguments
auxdata.NMuscles = DatStore.nMuscles;   % number of muscles
auxdata.Ndof = DatStore.nDOF;           % humber of dofs
% DatStore.time = DatStore.time;          % time window
auxdata.ID = DatStore.T_exp;            % inverse dynamics
auxdata.params = DatStore.params;       % Muscle-tendon parameters

auxdata.ankle_clutched_spring = Misc.ankle_clutched_spring;
if Misc.ankle_clutched_spring
    % TODO support separately clutching left and right leg.
    auxdata.clutched_spring_dofs = strmatch('ankle_angle', DatStore.DOFNames);
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
e0 = 0.6; kpe = 4; t50 = exp(kpe * (0.2 - 0.10e1) / e0);
pp1 = (t50 - 0.10e1); t7 = exp(kpe); pp2 = (t7 - 0.10e1);
auxdata.Fpparam = [pp1;pp2;Misc.tendonStiffnessCoeff];

% Problem bounds 
e_min = 0; e_max = 1;                   % bounds on muscle excitation
a_min = 0; a_max = 1;                   % bounds on muscle activation
vMtilde_min = -1; vMtilde_max = 1;      % bounds on normalized muscle fiber velocity
lMtilde_min = 0.2; lMtilde_max = 1.8;   % bounds on normalized muscle fiber length

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
umin = e_min*ones(1,auxdata.NMuscles); 
umax = e_max*ones(1,auxdata.NMuscles);
vMtildemin = vMtilde_min*ones(1,auxdata.NMuscles); 
vMtildemax = vMtilde_max*ones(1,auxdata.NMuscles);
aTmin = -1*ones(1,auxdata.Ndof); 
aTmax = 1*ones(1,auxdata.Ndof);
switch study{2}
    case 'HipAnkle'
        aDmin = 0; aDmax = 1;
        control_bounds_lower = [umin aTmin vMtildemin aDmin];
        control_bounds_upper = [umax aTmax vMtildemax aDmax];
    case 'HipKneeAnkle'
        aDmin = 0; aDmax = 1;
        control_bounds_lower = [umin aTmin vMtildemin aDmin];
        control_bounds_upper = [umax aTmax vMtildemax aDmax];
    case 'HipExtHipAbd'
        aDmin = 0; aDmax = 1;
        control_bounds_lower = [umin aTmin vMtildemin aDmin];
        control_bounds_upper = [umax aTmax vMtildemax aDmax];
    otherwise
        control_bounds_lower = [umin aTmin vMtildemin];
        control_bounds_upper = [umax aTmax vMtildemax];
end
bounds.phase.control.lower = control_bounds_lower; 
bounds.phase.control.upper = control_bounds_upper;

% States bounds
actMin = a_min*ones(1,auxdata.NMuscles); 
actMax = a_max*ones(1,auxdata.NMuscles);
lMtildemin = lMtilde_min*ones(1,auxdata.NMuscles); 
lMtildemax = lMtilde_max*ones(1,auxdata.NMuscles);
states_bounds_lower = [actMin, lMtildemin];
states_bounds_upper = [actMax, lMtildemax];
bounds.phase.initialstate.lower = states_bounds_lower; 
bounds.phase.initialstate.upper = states_bounds_upper;
bounds.phase.state.lower = states_bounds_lower; 
bounds.phase.state.upper = states_bounds_upper;
bounds.phase.finalstate.lower = states_bounds_lower; 
bounds.phase.finalstate.upper = states_bounds_upper;

% Integral bounds
bounds.phase.integral.lower = 0;
bounds.phase.integral.upper = 10000*(tf-t0);

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
    otherwise
        % No parameter case
end

% Path constraints
HillEquil = zeros(1, auxdata.NMuscles);
ID_bounds = zeros(1, auxdata.Ndof);
bounds.phase.path.lower = [ID_bounds,HillEquil];
bounds.phase.path.upper = [ID_bounds,HillEquil];

if Misc.ankle_clutched_spring
    stiffness_lower = 0;
    stiffness_upper = 1;
    if Misc.ankle_clutched_spring_stiffness ~= -1
        stiffness_lower = Misc.ankle_clutched_spring_stiffness;
        stiffness_upper = Misc.ankle_clutched_spring_stiffness;
    end
    bounds.parameter.lower = [stiffness_lower, -0.5];
    bounds.parameter.upper = [stiffness_upper, 0.5];
end

% Eventgroup
% Impose mild periodicity
pera_lower = -1 * ones(1, auxdata.NMuscles);
pera_upper = 1 * ones(1, auxdata.NMuscles);
perlMtilde_lower = -1*ones(1,auxdata.NMuscles);
perlMtilde_upper = 1*ones(1,auxdata.NMuscles);
bounds.eventgroup.lower = [pera_lower perlMtilde_lower]; 
bounds.eventgroup.upper = [pera_upper perlMtilde_upper];

% Initial guess
N = length(DatStore.time);
guess.phase.time = DatStore.time;
switch study{2}
    case 'HipAnkle'
        guess.phase.control = [DatStore.SoAct DatStore.SoRAct./150 0.01*ones(N,auxdata.NMuscles) 0.5*ones(N,1)];
    case 'HipKneeAnkle'
        guess.phase.control = [DatStore.SoAct DatStore.SoRAct./150 0.01*ones(N,auxdata.NMuscles) 0.5*ones(N,1)];
    case 'HipExtHipAbd'
        guess.phase.control = [DatStore.SoAct DatStore.SoRAct./150 0.01*ones(N,auxdata.NMuscles) 0.5*ones(N,1)];
    otherwise
        guess.phase.control = [DatStore.SoAct DatStore.SoRAct./150 0.01*ones(N,auxdata.NMuscles)];
end

guess.phase.state =  [DatStore.SoAct ones(N,auxdata.NMuscles)];
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
end

if Misc.ankle_clutched_spring
    guess.parameter = [0.5, 0];
end

% Empty exosuit force and torque data structures
DatStore.T_exo = zeros(length(DatStore.time),auxdata.Ndof);
DatStore.Fopt_exo = zeros(auxdata.Ndof,1);
DatStore.Fopt_exo_knee = zeros(auxdata.Ndof,1);

% "Tradeoff" struct
%  +/- 1 for the two joints you are investigating
%  zeros for everything else
DatStore.tradeoff = zeros(auxdata.Ndof,1);

% Reproduce Quinlivan et al. 2017 study
if strcmp(study{2},'Quinlivan2017') || strcmp(study{2},'Q2017')
    disp('debug')
    % Exosuit moment curves
    ExoCurves = load(['/Examples/SoftExosuitDesign/Quinlivan2017/' folder '/ExoCurves.mat']);
    exoTime = ExoCurves.time;
    % Peaks are body mass normalized so multiply by model mass
    exoAnkleMomentPeaks = ExoCurves.am_peak * Misc.model_mass;
    exoAnkleNormalizedMoment = ExoCurves.am_norm;
    exoHipMomentPeaks = ExoCurves.hm_peak * Misc.model_mass;
    exoHipNormalizedMoment = ExoCurves.hm_norm;

    % Interpolate exosuit moments to match data
    if Misc.exo_force_level    
        exoAnkleMoment = exoAnkleMomentPeaks(Misc.exo_force_level) * exoAnkleNormalizedMoment;
        exoHipMoment = exoHipMomentPeaks(Misc.exo_force_level) * exoHipNormalizedMoment;
        for dof = 1:auxdata.Ndof
            if strcmp('ankle_angle_r', DatStore.DOFNames{dof})
                % Negative to match ankle_angle_r coord convention
                DatStore.T_exo(:,dof) = -interp1(exoTime, exoAnkleMoment, DatStore.time);    
            elseif strcmp('hip_flexion_r', DatStore.DOFNames{dof})
                % Positive to match hip_flexion_r coord convention
                DatStore.T_exo(:,dof) = interp1(exoTime, exoHipMoment, DatStore.time);                
            end
        end
    end
end

% Given one actuator, compare tradeoff between hip flexion and ankle 
%  plantarflexion assistance
if strcmp(study{2},'HipAnkle') 
    % Exosuit moment curves
    ExoCurves = load('/Examples/SoftExosuitDesign/HipAnkle/ExoCurves.mat');
    % Peaks are body mass normalized so multiply by model mass
    exoAnkleForcePeaks = ExoCurves.af_peak * Misc.model_mass;

    % Interpolate exosuit moments to match data
    if Misc.exo_force_level    
        exoForce = exoAnkleForcePeaks(Misc.exo_force_level);
        for dof = 1:auxdata.Ndof
            if strcmp('ankle_angle_r', DatStore.DOFNames{dof})
                % Negative to match ankle_angle_r coord convention
                DatStore.Fopt_exo(dof) = -exoForce;
                DatStore.tradeoff(dof) = -1;
            elseif strcmp('hip_flexion_r', DatStore.DOFNames{dof})
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
if strcmp(study{2},'HipKneeAnkle') 
    % Exosuit moment curves
    ExoCurves = load('/Examples/SoftExosuitDesign/HipAnkle/ExoCurves.mat');
    % Peaks are body mass normalized so multiply by model mass
    exoAnkleForcePeaks = ExoCurves.af_peak * Misc.model_mass;

    % Interpolate exosuit moments to match data
    if Misc.exo_force_level    
        exoForce = exoAnkleForcePeaks(Misc.exo_force_level);
        for dof = 1:auxdata.Ndof
            if strcmp('ankle_angle_r', DatStore.DOFNames{dof})
                % Negative to match ankle_angle_r coord convention
                DatStore.Fopt_exo(dof) = -exoForce;
                DatStore.tradeoff(dof) = -1;
            elseif strcmp('hip_flexion_r', DatStore.DOFNames{dof})
                % Positive to match hip_flexion_r coord convention
                DatStore.Fopt_exo(dof) = exoForce;
                DatStore.tradeoff(dof) = 1;
            elseif strcmp('knee_angle_r', DatStore.DOFNames{dof})
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
    ExoCurves = load('/Examples/SoftExosuitDesign/HipAnkle/ExoCurves.mat');
    % Peaks are body mass normalized so multiply by model mass
    exoAnkleForcePeaks = ExoCurves.af_peak * Misc.model_mass;

    % Interpolate exosuit moments to match data
    if Misc.exo_force_level    
        exoForce = exoAnkleForcePeaks(Misc.exo_force_level);
        for dof = 1:auxdata.Ndof
            if strcmp('hip_flexion_r', DatStore.DOFNames{dof})
                % Negative to match hip_flexion_r coord convention
                DatStore.Fopt_exo(dof) = -exoForce;
                DatStore.tradeoff(dof) = 1;
            elseif strcmp('hip_adduction_r',DatStore.DOFNames{dof})
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
    ExoCurves = load('/Examples/SoftExosuitDesign/HipAnkleMass/ExoCurves.mat');
    exoTime = ExoCurves.time;
    % Peaks are body mass normalized so multiply by model mass
    exoAnkleMomentPeaks = ExoCurves.am_peak * Misc.model_mass;
    exoAnkleNormalizedMoment = ExoCurves.am_norm;
    exoHipMomentPeaks = ExoCurves.hm_peak * Misc.model_mass;
    exoHipNormalizedMoment = ExoCurves.hm_norm;
 
    % Interpolate exosuit moments to match data
    if Misc.exo_force_level    
        exoAnkleMoment = exoAnkleMomentPeaks(Misc.exo_force_level) * exoAnkleNormalizedMoment;
        exoHipMoment = exoHipMomentPeaks(Misc.exo_force_level) * exoHipNormalizedMoment;
        for dof = 1:auxdata.Ndof
            if strcmp('ankle_angle_r', DatStore.DOFNames{dof})
                % Negative to match ankle_angle_r coord convention
                DatStore.T_exo(:,dof) = -interp1(exoTime, exoAnkleMoment, DatStore.time); 
                DatStore.tradeoff(dof) = -1;
            elseif strcmp('hip_flexion_r', DatStore.DOFNames{dof})
                % Positive to match hip_flexion_r coord convention
                DatStore.T_exo(:,dof) = interp1(exoTime, exoHipMoment, DatStore.time);
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
setup.name = 'DynamicOptimization_lMtildestate';
setup.auxdata = auxdata;
setup.bounds = bounds;
setup.guess = guess;
setup.nlp.solver = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.derivatives.derivativelevel = 'second';
setup.nlp.ipoptoptions.tolerance = 1e-6;
setup.nlp.ipoptoptions.maxiterations = 2000;
setup.derivatives.supplier = 'adigator';
setup.scales.method = 'none';
setup.mesh.method = 'hp-PattersonRao';
setup.mesh.tolerance = 1e-3;
setup.mesh.maxiterations = 20;
setup.mesh.colpointsmin = 3;
setup.mesh.colpointsmax = 10;
setup.method = 'RPM-integration';
setup.displaylevel = 2;
NMeshIntervals = round((tf-t0)*Misc.Mesh_Frequency);
setup.mesh.phase.colpoints = 3*ones(1,NMeshIntervals);
setup.mesh.phase.fraction = (1/(NMeshIntervals))*ones(1,NMeshIntervals);
setup.functions.continuous = str2func(['continous_lMtilde' tag]);
setup.functions.endpoint = str2func(['endpoint_lMtilde' tag]);

    
% ADiGator setup
persistent splinestruct
input.auxdata = auxdata;

pathLock='/tmp/adigator3.lock'
% Try to create and lock this file.
if ~system(sprintf('lockfile %s',pathLock))
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
    adigatorGenFiles4gpops2(setup)

    % Now remove the lockfile
    system(sprintf('rm -f %s',pathLock));
end

setup.functions.continuous = str2func(['Wrap4continous_lMtilde' tag]);
setup.adigatorgrd.continuous = str2func(['continous_lMtilde' tag 'GrdWrap']);
setup.adigatorgrd.endpoint   = str2func(['endpoint_lMtilde' tag 'ADiGatorGrd']);
setup.adigatorhes.continuous = str2func(['continous_lMtilde' tag 'HesWrap']);
setup.adigatorhes.endpoint   = str2func(['endpoint_lMtilde' tag 'ADiGatorHes']);


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
lMtilde = res.state(:,auxdata.NMuscles+1:auxdata.NMuscles*2);
lM = lMtilde.*(ones(length(Time),1)*DatStore.lOpt);
MExcitation = res.control(:,1:auxdata.NMuscles);
RActivation = res.control(:,auxdata.NMuscles+1:auxdata.NMuscles+auxdata.Ndof);
OptInfo = output;

% Tendon force from lMtilde
% Interpolation lMT
lMTinterp = interp1(DatStore.time,DatStore.LMT,Time);
[TForcetilde,TForce] = TendonForce_lMtilde(...
    lMtilde,auxdata.params,lMTinterp);

end

