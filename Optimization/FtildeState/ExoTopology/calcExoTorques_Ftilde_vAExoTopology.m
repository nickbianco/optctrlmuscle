function [ExoTorques_Act, MomentArms_Act] = calcExoTorques_Ftilde_vAExoTopology(OptInfo,DatStore)

time = OptInfo.result.solution.phase.time;
numColPoints = length(time);
auxdata = OptInfo.result.setup.auxdata;
Ndof = auxdata.Ndof;

% Get device control
aD = OptInfo.result.solution.phase.control(:,end-(auxdata.numDeviceDOFs-1):end);

% Get net joint moments
T_exp = interp1(DatStore.time, DatStore.T_exp, time);

% Convert parameters to the correct range
paramsLower = auxdata.paramsLower;
paramsUpper = auxdata.paramsUpper;
parameter = OptInfo.result.solution.parameter;
parameter = (paramsUpper-paramsLower).*parameter + 0.5*(paramsLower+paramsUpper);

% Get moment arms
exoMomentArms = zeros(numColPoints,3);
aD_hip = zeros(numColPoints,1);
aD_knee = zeros(numColPoints,1);
aD_ankle = zeros(numColPoints,1);
if auxdata.device.hip
    exoMomentArms(:,1) = parameter(:,auxdata.device.hip);
    if auxdata.numDeviceDOFs > 1
        aD_hip = aD(:,auxdata.device.hip);
    else
        aD_hip = aD;
    end
end
if auxdata.device.knee
    exoMomentArms(:,2) = parameter(:,auxdata.device.knee);
    if auxdata.numDeviceDOFs > 1
        aD_knee = aD(:,auxdata.device.knee);
    else
        aD_knee = aD;
    end
end
if auxdata.device.ankle
    exoMomentArms(:,3) = parameter(:,auxdata.device.ankle);
    if auxdata.numDeviceDOFs > 1
        aD_ankle = aD(:,auxdata.device.ankle);
    else
        aD_ankle = aD;
    end
    if auxdata.shift_exo_peaks
        T_hip_ref = sign(exoMomentArms(1,1))*T_exp(:,auxdata.hip_DOF);
        T_ankle_ref = sign(exoMomentArms(1,3))*T_exp(:,auxdata.ankle_DOF);
        [~,hipIdx] = max(T_hip_ref);
        [~,ankleIdx] = max(T_ankle_ref);
        shift_factor = ankleIdx - hipIdx;
        aD_ankle = ShiftCurve(aD_ankle, shift_factor);
    end 
end
MomentArms_Act = exoMomentArms(1,:);

% Exosuit torques
Texo_act_hip = auxdata.Tmax_act_hip.*aD_hip.*exoMomentArms(:,1);
Texo_act_knee = auxdata.Tmax_act_knee.*aD_knee.*exoMomentArms(:,2).*auxdata.kneeAngleSign;
Texo_act_ankle = auxdata.Tmax_act_ankle.*aD_ankle.*exoMomentArms(:,3);

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