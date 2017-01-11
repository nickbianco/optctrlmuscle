%% Example
% add main folder and subfolder to matlab path (installation)
filepath=which('gait10dof18musc_ankle_clutched_spring.m');
[DirExample,~,~]=fileparts(filepath);
[DirExample2,~,~]=fileparts(DirExample);
[MainDir,~]=fileparts(DirExample2);
addpath(genpath(MainDir));

% Needed Input Arguments
if isempty(getenv('OPENSIM_HOME'))
    error('You must define the OPENSIM_HOME environment variable.');
end
Datapath = fullfile(getenv('OPENSIM_HOME'), 'Models', 'Gait10dof18musc', ...
    'OutputReference');
IK_path=fullfile(Datapath,'IK','subject01_walk_IK.mot');
ID_path=fullfile(Datapath,'ID','inversedynamics.sto');
%ID_path=[]; % compute ID from the external loads
model_path=fullfile(Datapath,'subject01.osim');
time=[0.7 1.4];     % Part of the right stance phase
OutPath=fullfile(DirExample,'Results');

Misc.DofNames_Input={'ankle_angle_r','knee_angle_r','hip_flexion_r'};
Misc.Loads_path=fullfile(Datapath,'ExperimentalData','subject01_walk_grf.xml');

Misc.costfun = 'Exc_Act';
Misc.ankle_clutched_spring = true;

%% Solve the problem
[Time,MExcitation,MActivation,RActivation,TForcetilde,TForce,lMtilde,lM,MuscleNames,OptInfo,DatStore]=SolveMuscleRedundancy_lMtildeState(model_path,IK_path,ID_path,time,OutPath,Misc);
