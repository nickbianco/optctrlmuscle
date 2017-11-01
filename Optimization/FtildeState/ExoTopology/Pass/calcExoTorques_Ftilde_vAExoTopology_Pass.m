function [ExoTorques_Pass, MomentArms_Pass, Fexo_pass, Lexo, IK, slackLength] = calcExoTorques_Ftilde_vAExoTopology_Pass(OptInfo,DatStore)

time = OptInfo.result.solution.phase.time;
numColPoints = length(time);
auxdata = OptInfo.result.setup.auxdata;
Ndof = auxdata.Ndof;

% Get moment arms
parameter = OptInfo.result.solution.parameter;
exoMomentArms = zeros(numColPoints,3);
if auxdata.passive.hip
    exoMomentArms(:,1) = parameter(:,auxdata.passive.hip);
end
if auxdata.passive.knee
    exoMomentArms(:,2) = parameter(:,auxdata.passive.knee);
end
if auxdata.passive.ankle
    exoMomentArms(:,3) = parameter(:,auxdata.passive.ankle);
end
MomentArms_Pass = exoMomentArms(1,:);

% Slack length of passive elastic device
exoSlackLength = parameter(:,end);
slackLength = exoSlackLength(1);

% Exosuit path length
Lexo = ones(numColPoints,1);
IK = interp1(DatStore.time, (pi/180)*DatStore.q_exp, time);
for dof = 1:Ndof
    if auxdata.passive.hip && (dof==auxdata.hip_DOF)
        Lexo = Lexo + -exoMomentArms(:,1).*IK(:,auxdata.hip_DOF);
    end
    if auxdata.passive.knee && (dof==auxdata.knee_DOF)
        Lexo = Lexo + -exoMomentArms(:,2).*IK(:,auxdata.knee_DOF)*auxdata.kneeAngleSign;
    end
    if auxdata.passive.ankle && (dof==auxdata.ankle_DOF)
        Lexo = Lexo + -exoMomentArms(:,3).*IK(:,auxdata.ankle_DOF);
    end
end

% Calculate passive force based on normalized exo path length
k = auxdata.passiveStiffness;
positiveStiffnessAboveLslack = (1 ./ (1 + exp(100 * ((exoSlackLength+0.075) - Lexo))));
Fexo_pass = k*(Lexo - exoSlackLength) .* positiveStiffnessAboveLslack;
Texo_pass_hip = Fexo_pass.*exoMomentArms(:,1);
Texo_pass_knee = Fexo_pass.*exoMomentArms(:,2);
Texo_pass_ankle = Fexo_pass.*exoMomentArms(:,3);

ExoTorques_Pass = zeros(length(time), Ndof);
for dof = 1:Ndof
    if dof==auxdata.hip_DOF
        ExoTorques_Pass(:,dof) = Texo_pass_hip;
    end
    if dof==auxdata.knee_DOF
        ExoTorques_Pass(:,dof) = Texo_pass_knee;
    end
    if dof==auxdata.ankle_DOF
        ExoTorques_Pass(:,dof) = Texo_pass_ankle;
    end
end


end