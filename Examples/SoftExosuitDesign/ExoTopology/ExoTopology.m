% add main folder and subfolder to matlab path (installation)
filepath=which('ExoTopology.m');
[DirCurrent,~,~]=fileparts(filepath);
[DirSoftExosuitDesign,~]=fileparts(DirCurrent);
[DirExamples,~]=fileparts(DirSoftExosuitDesign);
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

Misc.DofNames_Input={'hip_flexion_r','knee_angle_r','ankle_angle_r'};
Misc.Loads_path=fullfile(Datapath,'ExperimentalData','subject01_walk_grf.xml');

% Optional Input Arguments
Misc.costfun = 'Default';
Misc.study = 'SoftExosuitDesign/Topology';
Misc.activeDOFs = {'hip','knee','ankle'};
Misc.passiveDOFs = {'hip','knee','ankle'};

tag = '';
if ~isempty(Misc.activeDOFs)
    Misc.costfun = [Misc.costfun 'Act'];
    active_tag = 'a';
    for i = 1:length(Misc.activeDOFs)
        switch Misc.activeDOFs{i}
            case 'hip'
                active_tag = [active_tag 'H'];
            case 'knee'
                active_tag = [active_tag 'K'];
            case 'ankle'
                active_tag = [active_tag 'A'];
        end
    end
    tag = [tag '_' active_tag];
end

if ~isempty(Misc.passiveDOFs)
    Misc.costfun = [Misc.costfun 'Pass'];
    passive_tag = 'p';
    for i = 1:length(Misc.passiveDOFs)
        switch Misc.passiveDOFs{i}
            case 'hip'
                passive_tag = [passive_tag 'H'];
            case 'knee'
                passive_tag = [passive_tag 'K'];
            case 'ankle'
                passive_tag = [passive_tag 'A'];
        end
    end
    tag = [tag '_' passive_tag];
end

%% Solve the problem
[Time,MExcitation,MActivation,RActivation,TForcetilde,TForce,lMtilde,lM,MuscleNames,OptInfo,DatStore] =  ... 
    SolveMuscleRedundancy_FtildeState_actdyn(model_path,IK_path,ID_path,time,OutPath,Misc);
filename='ExoTopology_MRS_solution_opt.mat';
savepath=fullfile(DirCurrent,filename);
save(savepath,'Time','MExcitation','MActivation','RActivation','TForcetilde', ...
        'TForce','lMtilde','lM','MuscleNames','OptInfo','DatStore');
    