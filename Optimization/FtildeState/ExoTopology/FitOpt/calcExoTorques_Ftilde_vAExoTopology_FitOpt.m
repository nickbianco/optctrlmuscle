function [ExoTorques_FitOpt] = calcExoTorques_Ftilde_vAExoTopology_FitOpt(OptInfo,DatStore)

time = OptInfo.result.solution.phase.time;
numColPoints = length(time);
auxdata = OptInfo.result.setup.auxdata;
Ndof = auxdata.Ndof;

% Get moment arms
aD_hip = zeros(numColPoints,1);
aD_knee = zeros(numColPoints,1);
aD_ankle = zeros(numColPoints,1);
if auxdata.active.hip
    aD_hip = ppval(auxdata.JointEXOSpline(auxdata.hip_DOF),time);
end
if auxdata.active.knee
    aD_knee = ppval(auxdata.JointEXOSpline(auxdata.knee_DOF),time);
end
if auxdata.active.ankle
    aD_ankle = ppval(auxdata.JointEXOSpline(auxdata.ankle_DOF),time);
end

% Exosuit torques
Texo_act_hip = auxdata.Tmax_act*aD_hip;
Texo_act_knee = auxdata.Tmax_act*aD_knee; %.*auxdata.kneeAngleSign;
Texo_act_ankle = auxdata.Tmax_act*aD_ankle;

ExoTorques_FitOpt = zeros(length(time), Ndof);
for dof = 1:Ndof
    if dof==auxdata.hip_DOF
        ExoTorques_FitOpt(:,dof) = Texo_act_hip;
    end
    if dof==auxdata.knee_DOF
        ExoTorques_FitOpt(:,dof) = Texo_act_knee;
    end
    if dof==auxdata.ankle_DOF
        ExoTorques_FitOpt(:,dof) = Texo_act_ankle;
    end
end

end