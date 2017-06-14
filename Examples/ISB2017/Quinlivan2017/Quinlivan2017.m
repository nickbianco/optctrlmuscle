clear all; close all; clc

%% Example
% add main folder and subfolder to matlab path (installation)
filepath=which('Quinlivan2017.m');
[DirCurrent,~,~]=fileparts(filepath);
[DirISB2017,~]=fileparts(DirCurrent);
[DirExamples,~]=fileparts(DirISB2017);
[MainDir,~]=fileparts(DirExamples);
addpath(genpath(MainDir));

% Needed Input Arguments
if isempty(getenv('OPENSIM_HOME'))
    error('You must define the OPENSIM_HOME environment variable.');
end
Datapath = fullfile(getenv('OPENSIM_HOME'), 'Models', 'Gait2354_Simbody', ...
    'OutputReference');
IK_path=fullfile(Datapath,'subject01_walk1_ik.mot');
ID_path=fullfile(Datapath,'ResultsInverseDynamics','inverse_dynamics.sto');
model_path=fullfile(Datapath,'subject01_scaledOnly.osim');
time=[0.6 1.4];     % Right gait cycle
OutPath=fullfile(DirCurrent,'Results');

%Misc.MuscleNames_Input={};      % Selects all muscles for the Input DOFS when this is left empty.
Misc.DofNames_Input={'ankle_angle_r','knee_angle_r','hip_flexion_r','hip_adduction_r','hip_rotation_r'};
Misc.Loads_path = fullfile(getenv('OPENSIM_HOME'), 'Models', 'Gait2354_Simbody','subject01_walk1_grf.xml');

% Optional Input Arguments
Misc.costfun = 'Exc_Act';
Misc.study = 'ISB2017/Quinlivan2017';
Misc.model_mass = 75.1646; % kg (Gait2354 mass)

%% Solve the problem

% Device optimization condition
[Time,MExcitation,MActivation,RActivation,TForcetilde,TForce,lMtilde,lM,MuscleNames,OptInfo,DatStore]=SolveMuscleRedundancy_lMtildeState(model_path,IK_path,ID_path,time,OutPath,Misc);
ExoTorques = getExoTorques(OptInfo,DatStore,Misc);
filename = 'Quinlivan2017_MRS_solution_opt.mat';
savepath=fullfile(DirCurrent,filename);
save(savepath,'Time','MExcitation','MActivation','RActivation','TForcetilde', ...
        'TForce','lMtilde','lM','MuscleNames','OptInfo','DatStore','ExoTorques');

% No device condition
Misc.exo_force_level = 0;
[Time,MExcitation,MActivation,RActivation,TForcetilde,TForce,lMtilde,lM,MuscleNames,OptInfo,DatStore]=SolveMuscleRedundancy_lMtildeState(model_path,IK_path,ID_path,time,OutPath,Misc);
ExoTorques = getExoTorques(OptInfo,DatStore,Misc);
filename=strcat('Quinlivan2017_MRS_solution_force_level_',int2str(Misc.exo_force_level),'.mat');
savepath=fullfile(DirCurrent,filename);
save(savepath,'Time','MExcitation','MActivation','RActivation','TForcetilde', ...
        'TForce','lMtilde','lM','MuscleNames','OptInfo','DatStore','ExoTorques');
       
keyboard

for i = 1:4
    % Device force level
    % 1 --> MIN
    % 2 --> MED
    % 3 --> HIGH
    % 4 --> MAX
    Misc.exo_force_level = i;
    [Time,MExcitation,MActivation,RActivation,TForcetilde,TForce,lMtilde,lM,MuscleNames,OptInfo,DatStore]=SolveMuscleRedundancy_lMtildeState(model_path,IK_path,ID_path,time,OutPath,Misc);
    ExoTorques = getExoTorques(OptInfo,DatStore,Misc);
    filename=strcat('Quinlivan2017_MRS_solution_force_level_',int2str(Misc.exo_force_level),'.mat');
    savepath=fullfile(DirCurrent,filename);
    save(savepath,'Time','MExcitation','MActivation','RActivation','TForcetilde', ...
            'TForce','lMtilde','lM','MuscleNames','OptInfo','DatStore','ExoTorques');
end

function ExoTorques = getExoTorques(OptInfo,DatStore,Misc)

auxdata = OptInfo.result.setup.auxdata;
p = auxdata.p_linreg;
exo_force_level = OptInfo.result.solution.parameter;
for dof = 1:length(Misc.DofNames_Input)
    exoPeak = p(1,dof)*exo_force_level + p(2,dof);
    exoNormTorque = DatStore.T_exo(:,dof);
    if contains(Misc.DofNames_Input{dof},'ankle_angle')
        ExoTorques.ankle_angle = exoPeak*exoNormTorque;
    elseif contains(Misc.DofNames_Input{dof},'hip_flexion')
        ExoTorques.hip_flexion = exoPeak*exoNormTorque;
    end   
end

end
