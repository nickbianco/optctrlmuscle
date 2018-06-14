function [ExoTorques_Exp, MomentArms_Exp] = calcExoTorques_Ftilde_vAExoTopology_Exp(OptInfo,DatStore)

time = OptInfo.result.solution.phase.time;
numColPoints = length(time);
auxdata = OptInfo.result.setup.auxdata;
Ndof = auxdata.Ndof;

% Get moment arms
parameter = OptInfo.result.solution.parameter;
exoMomentArms = zeros(numColPoints,3);
aD_hip = zeros(numColPoints,1);
aD_knee = zeros(numColPoints,1);
aD_ankle = zeros(numColPoints,1);
if auxdata.active.hip
    aD_hip = ppval(auxdata.JointEXOSpline(auxdata.hip_DOF),time);
    if auxdata.same_torque_gain
        exoMomentArms(:,1) = sign(auxdata.paramsUpper(auxdata.active.hip))*parameter;
    else
        exoMomentArms(:,1) = parameter(:,auxdata.active.hip);
    end
end
if auxdata.active.knee
    aD_knee = ppval(auxdata.JointEXOSpline(auxdata.knee_DOF),time);
    if auxdata.same_torque_gain
        exoMomentArms(:,2) = sign(auxdata.paramsUpper(auxdata.active.knee))*parameter;
    else
        exoMomentArms(:,2) = parameter(:,auxdata.active.knee);
    end
end
if auxdata.active.ankle
    aD_ankle = ppval(auxdata.JointEXOSpline(auxdata.ankle_DOF),time);
    if auxdata.same_torque_gain
        exoMomentArms(:,3) = sign(auxdata.paramsUpper(auxdata.active.ankle))*parameter;
    else
        exoMomentArms(:,3) = parameter(:,auxdata.active.ankle);
    end
end
MomentArms_Exp = exoMomentArms(1,:);

% Exosuit torques
Texo_act_hip = auxdata.Tmax_act*aD_hip.*exoMomentArms(:,1);
Texo_act_knee = auxdata.Tmax_act*aD_knee.*exoMomentArms(:,2)*auxdata.kneeAngleSign;
Texo_act_ankle = auxdata.Tmax_act*aD_ankle.*exoMomentArms(:,3);

ExoTorques_Exp = zeros(length(time), Ndof);
for dof = 1:Ndof
    if dof==auxdata.hip_DOF
        ExoTorques_Exp(:,dof) = Texo_act_hip;
    end
    if dof==auxdata.knee_DOF
        ExoTorques_Exp(:,dof) = Texo_act_knee;
    end
    if dof==auxdata.ankle_DOF
        ExoTorques_Exp(:,dof) = Texo_act_ankle;
    end
end

end