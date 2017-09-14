function [ExoTorques_Act] = calcExoTorques_Ftilde_vAExoTopology_Act(OptInfo,DatStore)

time = OptInfo.result.solution.phase.time;
auxdata = OptInfo.result.setup.auxdata;
Ndof = auxdata.Ndof;

% Get device control
aD = OptInfo.result.solution.phase.control(:,end);

% Get moment arms
parameter = OptInfo.result.solution.phase.parameter;
exoMomentArms = zeros(numColPoints,3);
if auxdata.active.hip
    exoMomentArms(:,1) = parameter(:,input.auxdata.active.hip);
end
if auxdata.active.knee
    exoMomentArms(:,2) = parameter(:,input.auxdata.active.knee);
end
if auxdata.active.ankle
    exoMomentArms(:,3) = parameter(:,input.auxdata.active.ankle);
end

% Exosuit torques
Texo_act_hip = auxdata.Fmax_act*aD.*exoMomentArms(:,1);
Texo_act_knee = auxdata.Fmax_act*aD.*exoMomentArms(:,2);
Texo_act_ankle = auxdata.Fmax_act*aD.*exoMomentArms(:,3);

ExoTorques_Act = zeros(length(time), Ndof);
for dof = 1:Ndof
    if dof==input.auxdata.hip_DOF
        ExoTorques_Act(:,dof) = Texo_act_hip;
    end
    if dof==input.auxdata.knee_DOF
        ExoTorques_Act(:,dof) = Texo_act_knee;
    end
    if dof==input.auxdata.ankle_DOF
        ExoTorques_Act(:,dof) = Texo_act_ankle;
    end
end

end