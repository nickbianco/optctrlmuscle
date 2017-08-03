function SynergyControl
%% Example
% add main folder and subfolder to matlab path (installation)
filepath=which('SynergyControl.m');
[DirCurrent,~,~]=fileparts(filepath);
[DirSynergyControl,~]=fileparts(DirCurrent);
[DirExamples,~]=fileparts(DirSynergyControl);
[MainDir,~]=fileparts(DirExamples);
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
OutPath=fullfile(DirCurrent,'Results');

Misc.DofNames_Input={'ankle_angle_r','knee_angle_r','hip_flexion_r'};
Misc.Loads_path=fullfile(Datapath,'ExperimentalData','subject01_walk_grf.xml');

% Optional Input Arguments
Misc.costfun = 'Exc_Act';
Misc.study = 'SynergyControl/HipAnkle';

%% Solve the problem
for numSyn = 2:6  
    Misc.synergy_control = numSyn;
    [Time,MExcitation,MActivation,RActivation,TForcetilde,TForce,lMtilde,lM,MuscleNames,OptInfo,DatStore] = SolveMuscleRedundancy_lMtildeState(model_path,IK_path,ID_path,time,OutPath,Misc);
    filename=['SynergyControl_MRS_solution_' num2str(numSyn) '_syns.mat'];
    savepath=fullfile(DirCurrent,filename);
    save(savepath,'Time','MExcitation','MActivation','RActivation','TForcetilde', ...
        'TForce','lMtilde','lM','MuscleNames','OptInfo','DatStore','ExoTorques');
end