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
IK_path='mrsdeviceopt_subject01_walk3_ik_solution.mot';
ID_path='mrsdeviceopt_subject01_walk3_id_solution.sto';
%ID_path=[]; % compute ID from the external loads
model_path='subject01.osim';
time=[0.7 1.4];     % Part of the right stance phase
OutPath=fullfile(DirCurrent,'Results');

Misc.DofNames_Input={'ankle_angle_r','knee_angle_r','hip_flexion_r'};
Misc.Loads_path='external_loads.xml';

% % Muscle excitation as control case
Misc.costfun = 'Exc_Act';
Misc.study = 'SoftExosuitDesign/HipAnkle';
Misc.exo_force_level = 2;

%% Solve the problem

[Time,MExcitation,MActivation,RActivation,TForcetilde,TForce,lMtilde,lM,MuscleNames,OptInfo,DatStore] = SolveMuscleRedundancy_lMtildeState(model_path,IK_path,ID_path,time,OutPath,Misc);
filename=['SynergyControl_MRS_solution_Inf_syns.mat'];
savepath=fullfile(DirCurrent,filename);
save(savepath,'Time','MExcitation','MActivation','RActivation','TForcetilde', ...
    'TForce','lMtilde','lM','MuscleNames','OptInfo','DatStore');

% Synergy contraint case
Misc.costfun = 'Exc_Act';
Misc.study = 'SynergyControl/HipAnkle';
Misc.exo_force_level = 2;

%% Solve the problem
for numSyn = 4:6  
    Misc.synergy_control = numSyn;
    [Time,MExcitation,RActivation,MuscleNames,OptInfo,DatStore] = SolveMuscleRedundancy_NoDynamics(model_path,IK_path,ID_path,time,OutPath,Misc);
    filename=['SynergyControl_MRS_solution_' num2str(numSyn) '_syns.mat'];
    savepath=fullfile(DirCurrent,filename);
    save(savepath,'Time','MExcitation','RActivation', ...
        'MuscleNames','OptInfo','DatStore');
end