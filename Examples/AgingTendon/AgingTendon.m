clear all; close all; clc

%% Choose formulation
formulation = 'lMtildeState';

%% Example
% add main folder and subfolder to matlab path (installation)
filepath=which('AgingTendon.m');
[DirCurrent,~,~]=fileparts(filepath); cd(DirCurrent);
[DirExamples,~]=fileparts(DirCurrent);
[MainDir,~]=fileparts(DirExamples);
addpath(genpath(MainDir));

% Needed Input Arguments
IK_path=fullfile(MainDir,'Data','WalkingData','inverse_kinematics.mot');
ID_path=fullfile(MainDir,'Data','WalkingData','inverse_dynamics.sto');
model_path=fullfile(MainDir,'Data','WalkingData','subject01.osim');
time=[0.516 1.95];     % Right stance phase (+50ms beginning and end of time interval, more details see manual and publication)
OutPath=fullfile(MainDir,'Examples','AgingTendon','Results');
Misc.MuscleNames_Input={};      % Selects all muscles for the Input DOFS when this is left empty.
Misc.DofNames_Input={'ankle_angle_r','knee_angle_r','hip_flexion_r'};

% Optional Input Arguments
Misc.Atendon = [];        % Tendon Stiffness for the selected muscles
Misc.f_cutoff_ID = 8;         % cutoff frequency filtering ID
Misc.f_order_ID = 5;             % order frequency filtering ID
Misc.f_cutoff_lMT = 8;         % cutoff frequency filtering lMT
Misc.f_order_lMT = 5;             % order frequency filtering lMT
Misc.f_cutoff_dM= 8;         % cutoff frequency filtering MA
Misc.f_order_dM = 5;             % order frequency filtering MA
Misc.f_cutoff_IK= 8;         % cutoff frequency filtering IK
Misc.f_order_IK = 5;             % order frequency filtering IK

%% Define problem
Misc.costfun = 'Exc_Act';
Misc.study = 'SoftExosuitDesign/HipAnkle';

%% Solve the problem

tendonMods = [1.0 0.9 0.8 0.7 0.6 0.5];
muscleMods = {'med_gas_r','lat_gas_r','soleus_r'};

for i = 1:length(tendonMods)
       
    for m = 1:length(muscleMods)
        Misc.tendonStiffnessModifiers.(muscleMods{m}) = tendonMods(i);
    end
    
    Misc.exo_force_level = 2;
    [Time,MExcitation,MActivation,RActivation,TForcetilde,TForce,lMtilde,lM,MuscleNames,OptInfo,DatStore] = ...
        SolveMuscleRedundancy_lMtildeState(model_path,IK_path,ID_path,time,OutPath,Misc);
    filename = ['AgingTendon_HipAnkle_stiffnessMod_' num2str(tendonMods(i)) '_MRS_solution.mat'];
    savepath=fullfile(DirCurrent,filename);
    save(savepath,'Time','MExcitation','MActivation','RActivation','TForcetilde', ...
        'TForce','lMtilde','lM','MuscleNames','OptInfo','DatStore','ExoTorques');
end

