function [ExoTorques_Pass, s] = calcExoTorques_Ftilde_vAExoTopology_Pass(OptInfo,DatStore)

time = OptInfo.result.solution.phase.time;
auxdata = OptInfo.result.setup.auxdata;
Ndof            = auxdata.Ndof;
splinestruct    = auxdata.splinestruct;

% Get slack variable
s = OptInfo.result.solution.phase.control(:,end);

% Get moment arms
exoMomentArms = zeros(numColPoints,3);
if input.auxdata.passive.hip
    exoMomentArms(:,1) = input.phase.parameter(:,auxdata.passive.hip);
end
if input.auxdata.passive.knee
    exoMomentArms(:,2) = input.phase.parameter(:,auxdata.passive.knee);
end
if input.auxdata.passive.ankle
    exoMomentArms(:,3) = input.phase.parameter(:,auxdata.passive.ankle);
end

% Slack length of passive elastic device
exoSlackLength = input.phase.parameter(:,end);

% Exosuit path length
Lexo = zeros(numColPoints,1);
for dof = 1:Ndof
    if input.auxdata.passive.hip && (dof==auxdata.hip_DOF)
        Lexo = Lexo + -exoMomentArms(:,1).*splinestruct.IK(:,input.auxdata.hip_DOF);
    end
    if input.auxdata.passive.knee && (dof==auxdata.knee_DOF)
        Lexo = Lexo + -exoMomentArms(:,2).*splinestruct.IK(:,input.auxdata.knee_DOF);
    end
    if input.auxdata.passive.ankle && (dof==auxdata.ankle_DOF)
        Lexo = Lexo + -exoMomentArms(:,3).*splinestruct.IK(:,input.auxdata.ankle_DOF);
    end
end

% Calculate passive force based on normalized exo path length
k = auxdata.passiveStiffness;
Fexo_pass = k*(Lexo - exoSlackLength) + s;
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