function HOBL
%% Example
% add main folder and subfolder to matlab path (installation)
filepath=which('HOBL.m');
[DirCurrent,~,~]=fileparts(filepath);
[DirISB2017,~]=fileparts(DirCurrent);
[DirExamples,~]=fileparts(DirISB2017);
[MainDir,~]=fileparts(DirExamples);
addpath(genpath(MainDir));

% Needed Input Arguments
if isempty(getenv('OPENSIM_HOME'))
    error('You must define the OPENSIM_HOME environment variable.');
end

IK_path=fullfile(DirCurrent,'subject02_running_RRA_Kinematics_q.mot');
ID_path=fullfile(DirCurrent,'subject02_running_RRA_Actuation_force.sto');
model_path=fullfile(DirCurrent,'subject02_running_RRA_cycle02_07_v2_18muscle.osim');
time=[0.7 1.4];    % left heel strike to left heel strike
OutPath=fullfile(DirCurrent,'Results');

%Misc.DofNames_Input={'ankle_angle_r','knee_angle_r','hip_flexion_r','hip_adduction_r','hip_rotation_r','ankle_angle_l','knee_angle_l','hip_flexion_l','hip_adduction_l','hip_rotation_l'};
Misc.DofNames_Input={'ankle_angle_l','knee_angle_l','hip_flexion_l'};
%Misc.Loads_path=fullfile(DirCurrent,'subject02_run_grf.xml');

% Optional Input Arguments
Misc.costfun = 'Exc_Act';
Misc.study = 'ISB2017/HOBL';

%% Solve the problem
[Time,MExcitation,MActivation,RActivation,TForcetilde,TForce,lMtilde,lM,MuscleNames,OptInfo,DatStore] = SolveMuscleRedundancy_lMtildeState(model_path,IK_path,ID_path,time,OutPath,Misc);
filename='HOBL_MRS_solution_opt.mat';
savepath=fullfile(DirCurrent,filename);
save(savepath,'Time','MExcitation','MActivation','RActivation','TForcetilde', ...
        'TForce','lMtilde','lM','MuscleNames','OptInfo','DatStore');
    
end