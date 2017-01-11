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
% Clutched spring for ankle plantarflexion (Collins 2015)
if ~isfield(Misc, 'ankle_clutched_spring') || isempty(Misc.ankle_clutched_spring)
    Misc.ankle_clutched_spring = false;
end
if ~isfield(Misc, 'phase_boundary') || isempty(Misc.phase_boundary)
    Misc.phase_boundary = NaN;
end

if Misc.ankle_clutched_spring 
    if ~strcmp(Misc.costfun, 'Exc_Act')
        error('ankle_clutched_spring == true requires costfun == ''Exc_Act''');
    end
    Misc.costfun = 'Exc_ActSpr';
end
if ~isnan(Misc.phase_boundary)
    assert(Misc.ankle_clutched_spring);
    boundary = Misc.phase_boundary;
    assert(time(1) < boundary);
    assert(boundary < time(2));
    numPhases = 2;
    Misc.costfun = 'Exc_ActPh';
else
    numPhases = 1;
end

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

% Parameters of passive muscle force-length characteristic
e0 = 0.6; kpe = 4; t50 = exp(kpe * (0.2 - 0.10e1) / e0);
pp1 = (t50 - 0.10e1); t7 = exp(kpe); pp2 = (t7 - 0.10e1);
auxdata.Fpparam = [pp1;pp2];

% Problem bounds 
e_min = 0; e_max = 1;                   % bounds on muscle excitation
a_min = 0; a_max = 1;                   % bounds on muscle activation
vMtilde_min = -1; vMtilde_max = 1;      % bounds on normalized muscle fiber velocity
lMtilde_min = 0.2; lMtilde_max = 1.8;   % bounds on normalized muscle fiber length

% Time bounds
t0 = DatStore.time(1); tf = DatStore.time(end);
if numPhases == 1
    bounds.phase.initialtime.lower = t0; bounds.phase.initialtime.upper = t0;
    bounds.phase.finaltime.lower = tf; bounds.phase.finaltime.upper = tf;
elseif numPhases == 2
    bounds.phase(1).initialtime.lower = t0;
    bounds.phase(1).initialtime.upper = t0;
    bounds.phase(1).finaltime.lower = Misc.phase_boundary;
    bounds.phase(1).finaltime.upper = Misc.phase_boundary;
    bounds.phase(2).initialtime.lower = Misc.phase_boundary;
    bounds.phase(2).initialtime.upper = Misc.phase_boundary;
    bounds.phase(2).finaltime.lower = tf;
    bounds.phase(2).finaltime.upper = tf;
else
    error('Invalid numPhases');
end
auxdata.initialtime = t0;
auxdata.finaltime = tf;
% Controls bounds
umin = e_min*ones(1,auxdata.NMuscles);
umax = e_max*ones(1,auxdata.NMuscles);
vMtildemin = vMtilde_min*ones(1,auxdata.NMuscles);
vMtildemax = vMtilde_max*ones(1,auxdata.NMuscles);
aTmin = -1*ones(1,auxdata.Ndof);
aTmax = 1*ones(1,auxdata.Ndof);
for ip = 1:numPhases
    bounds.phase(ip).control.lower = [umin aTmin vMtildemin];
    bounds.phase(ip).control.upper = [umax aTmax vMtildemax];
end
% States bunds
actMin = a_min*ones(1,auxdata.NMuscles);
actMax = a_max*ones(1,auxdata.NMuscles);
lMtildemin = lMtilde_min*ones(1,auxdata.NMuscles);
lMtildemax = lMtilde_max*ones(1,auxdata.NMuscles);
for ip = 1:numPhases
    bounds.phase(ip).initialstate.lower = [actMin, lMtildemin];
    bounds.phase(ip).initialstate.upper = [actMax, lMtildemax];
    bounds.phase(ip).state.lower = [actMin, lMtildemin];
    bounds.phase(ip).state.upper = [actMax, lMtildemax];
    bounds.phase(ip).finalstate.lower = [actMin, lMtildemin];
    bounds.phase(ip).finalstate.upper = [actMax, lMtildemax];
end
% Integral bounds
for ip = 1:numPhases
    bounds.phase(ip).integral.lower = 0;
    bounds.phase(ip).integral.upper = 10000*(tf-t0);
end

% Path constraints
HillEquil = zeros(1, auxdata.NMuscles);
ID_bounds = zeros(1, auxdata.Ndof);
for ip = 1:numPhases
    bounds.phase(ip).path.lower = [ID_bounds,HillEquil];
    bounds.phase(ip).path.upper = [ID_bounds,HillEquil];
end

if Misc.ankle_clutched_spring
    bounds.parameter.lower = [0, -0.5];
    bounds.parameter.upper = [1, 0.5];
end

% Eventgroup
% Impose mild periodicity
pera_lower = -1 * ones(1, auxdata.NMuscles);
pera_upper = 1 * ones(1, auxdata.NMuscles);
perlMtilde_lower = -1*ones(1,auxdata.NMuscles);
perlMtilde_upper = 1*ones(1,auxdata.NMuscles);
bounds.eventgroup(1).lower = [pera_lower perlMtilde_lower]; 
bounds.eventgroup(1).upper = [pera_upper perlMtilde_upper];
if numPhases == 2
    states_continuous = ...
            zeros(1, length(bounds.phase(1).state.lower));
            %+ ...
                     %length(bounds.phase(1).control.lower));
    bounds.eventgroup(2).lower = [states_continuous]; 
    bounds.eventgroup(2).upper = [states_continuous];
end

% Initial guess
N = length(DatStore.time);
if numPhases == 1
    guess.phase.time = DatStore.time;
    guess.phase.control = [DatStore.SoAct ...
                           DatStore.SoRAct./150 ...
                           0.01*ones(N,auxdata.NMuscles)];
    guess.phase.state =  [DatStore.SoAct ones(N,auxdata.NMuscles)];
    guess.phase.integral = 0;
elseif numPhases == 2
    phase1_indices = find(DatStore.time < Misc.phase_boundary);
    phase2_indices = find(DatStore.time >= Misc.phase_boundary);
    assert(length(phase1_indices) + length(phase2_indices) == N);
    phase_indices = {phase1_indices, phase2_indices};
    for ip = 1:numPhases
        idxs = phase_indices{ip};
        guess.phase(ip).time = DatStore.time(idxs);
        guess.phase(ip).control = [DatStore.SoAct(idxs, :) ...
                                   DatStore.SoRAct(idxs, :)./150 ...
                                   0.01*ones(length(idxs),auxdata.NMuscles)];
        guess.phase(ip).state =  [DatStore.SoAct(idxs, :) ones(length(idxs),auxdata.NMuscles)];
        guess.phase(ip).integral = 0;
    end
else
    error('Invalid numPhases');
end
if Misc.ankle_clutched_spring
    guess.parameter = [0.5, 0];
end

% Spline structures
for dof = 1:auxdata.Ndof
    for m = 1:auxdata.NMuscles       
        auxdata.JointMASpline(dof).Muscle(m) = spline(DatStore.time,auxdata.MA(dof).Joint(:,m));       
    end
    auxdata.JointIDSpline(dof) = spline(DatStore.time,DatStore.T_exp(:,dof));
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
setup.mesh.maxiterations = 5;
setup.mesh.colpointsmin = 3;
setup.mesh.colpointsmax = 10;
setup.method = 'RPM-integration';
setup.displaylevel = 2;
if numPhases == 1
    phaseDurations = {tf - t0};
elseif numPhases == 2
    phaseDurations = {Misc.phase_boundary - t0, tf - Misc.phase_boundary};
else
    error('Invalid numPhases.');
end
for ip = 1:numPhases
    NMeshIntervals = round(phaseDurations{ip} * Misc.Mesh_Frequency);
    setup.mesh.phase(ip).colpoints = 3*ones(1,NMeshIntervals);
    setup.mesh.phase(ip).fraction = (1/(NMeshIntervals))*ones(1,NMeshIntervals);
end
setup.functions.continuous = str2func(['musdynContinous_lMtildeState_' Misc.costfun]);
setup.functions.endpoint = str2func(['musdynEndpoint_lMtildeState_' Misc.costfun]);
    
% ADiGator setup
persistent splinestruct
input.auxdata = auxdata;
if numPhases == 1
    tdummy = guess.phase.time;
    splinestruct = SplineInputData(tdummy,input);
    splinenames = fieldnames(splinestruct);
    for Scount = 1:length(splinenames)
      secdim = size(splinestruct.(splinenames{Scount}),2);
      splinestructad.(splinenames{Scount}) = adigatorCreateAuxInput([Inf,secdim]);
      splinestruct.(splinenames{Scount}) = zeros(0,secdim);
    end
elseif numPhases == 2
    for ip = 1:numPhases
        tdummy = guess.phase(ip).time;
        splinestruct.phase(ip) = SplineInputData(tdummy,input);
        splinenames = fieldnames(splinestruct.phase(ip));
        for Scount = 1:length(splinenames)
          secdim = size(splinestruct.phase(ip).(splinenames{Scount}),2);
          splinestructad.phase(ip).(splinenames{Scount}) = adigatorCreateAuxInput([Inf,secdim]);
          splinestruct.phase(ip).(splinenames{Scount}) = zeros(0,secdim);
        end
    end
else
    error('Invalid numPhases.');
end
setup.auxdata.splinestruct = splinestructad;
adigatorGenFiles4gpops2(setup)

setup.functions.continuous = str2func(['Wrap4musdynContinous_lMtildeState_' Misc.costfun]);
setup.adigatorgrd.continuous = str2func(['musdynContinous_lMtildeState_' Misc.costfun 'GrdWrap']);
setup.adigatorgrd.endpoint   = str2func(['musdynEndpoint_lMtildeState_' Misc.costfun 'ADiGatorGrd']);
setup.adigatorhes.continuous = str2func(['musdynContinous_lMtildeState_' Misc.costfun 'HesWrap']);
setup.adigatorhes.endpoint   = str2func(['musdynEndpoint_lMtildeState_' Misc.costfun 'ADiGatorHes']);


%% ---------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% PART III: SOLVE OPTIMAL CONTROL PROBLEM ------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %

output = gpops2(setup);

MuscleNames = DatStore.MuscleNames;
Ntotal = 0;
for ip = 1:numPhases
    Ntotal = Ntotal + length(output.result.solution.phase(ip).time);
end
Time = zeros(Ntotal, 1);
MActivation = zeros(Ntotal, auxdata.NMuscles);
lMtilde = zeros(Ntotal, auxdata.NMuscles);
lM = zeros(Ntotal, auxdata.NMuscles);
MExcitation = zeros(Ntotal, auxdata.NMuscles);
RActivation = zeros(Ntotal, auxdata.Ndof);
TForcetilde = zeros(Ntotal, auxdata.NMuscles);
TForce = zeros(Ntotal, auxdata.NMuscles);
start_index = 1;
for ip = 1:numPhases
    idxs = start_index:(start_index-1+length(output.result.solution.phase(ip).time));

    res = output.result.solution.phase(ip);
    Time(idxs) = res.time;
    MActivation(idxs, :) = res.state(:,1:auxdata.NMuscles);
    lMtilde(idxs, :) = res.state(:,auxdata.NMuscles+1:auxdata.NMuscles*2);
    lM(idxs, :) = lMtilde(idxs, :).*(ones(length(idxs),1)*DatStore.lOpt);
    MExcitation(idxs, :) = res.control(:,1:auxdata.NMuscles);
    RActivation(idxs, :) = res.control(:,auxdata.NMuscles+1:auxdata.NMuscles+auxdata.Ndof);
    OptInfo = output;

    % Tendon force from lMtilde
    % Interpolation lMT
    lMTinterp = interp1(DatStore.time,DatStore.LMT,Time(idxs));
    [TForcetilde(idxs, :),TForce(idxs, :)] = TendonForce_lMtilde(...
            lMtilde(idxs, :),auxdata.params,lMTinterp);
    start_index = idxs(end) + 1;
end

end

