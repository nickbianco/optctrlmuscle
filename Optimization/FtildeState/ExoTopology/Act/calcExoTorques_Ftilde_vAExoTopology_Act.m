function [ExoTorques_Act, MomentArms_Act] = calcExoTorques_Ftilde_vAExoTopology_Act(OptInfo,DatStore)

time = OptInfo.result.solution.phase.time;
numColPoints = length(time);
auxdata = OptInfo.result.setup.auxdata;
Ndof = auxdata.Ndof;

% Get device control
aD = OptInfo.result.solution.phase.control(:,end-(auxdata.numActiveDOFs-1):end);

% Get moment arms
parameter = OptInfo.result.solution.parameter;
exoMomentArms = zeros(numColPoints,3);
aD_hip = zeros(numColPoints,1);
aD_knee = zeros(numColPoints,1);
aD_ankle = zeros(numColPoints,1);
if auxdata.active.hip
    exoMomentArms(:,1) = parameter(:,auxdata.active.hip);
    if auxdata.numActiveDOFs > 1
        aD_hip = aD(:,auxdata.active.hip);
    else
        aD_hip = aD;
    end
end
if auxdata.active.knee
    exoMomentArms(:,2) = parameter(:,auxdata.active.knee);
    if auxdata.numActiveDOFs > 1
        aD_knee = aD(:,auxdata.active.knee);
    else
        aD_knee = aD;
    end
end
if auxdata.active.ankle
    exoMomentArms(:,3) = parameter(:,auxdata.active.ankle);
    if auxdata.numActiveDOFs > 1
        aD_ankle = aD(:,auxdata.active.ankle);
    else
        aD_ankle = aD;
    end
end
MomentArms_Act = exoMomentArms(1,:);

% Exosuit torques
Texo_act_hip = auxdata.Tmax_act*aD_hip.*exoMomentArms(:,1);
Texo_act_knee = auxdata.Tmax_act*aD_knee.*exoMomentArms(:,2)*auxdata.kneeAngleSign;
Texo_act_ankle = auxdata.Tmax_act*aD_ankle.*exoMomentArms(:,3);

ExoTorques_Act = zeros(length(time), Ndof);
for dof = 1:Ndof
    if dof==auxdata.hip_DOF
        ExoTorques_Act(:,dof) = Texo_act_hip;
    end
    if dof==auxdata.knee_DOF
        ExoTorques_Act(:,dof) = Texo_act_knee;
    end
    if dof==auxdata.ankle_DOF
        ExoTorques_Act(:,dof) = Texo_act_ankle;
    end
end

end