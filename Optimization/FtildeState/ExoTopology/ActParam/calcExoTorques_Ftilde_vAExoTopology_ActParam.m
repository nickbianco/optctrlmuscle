function [ExoTorques_Act] = calcExoTorques_Ftilde_vAExoTopology_ActParam(OptInfo,DatStore)

time = OptInfo.result.solution.phase.time;
numColPoints = length(time);
auxdata = OptInfo.result.setup.auxdata;
Ndof = auxdata.Ndof;

% Get moment arms
parameter = OptInfo.result.solution.parameter;
signMoment_hip = 1;
signMoment_knee = 1;
signMoment_ankle = 1;
aD_hip = zeros(numColPoints,1);
aD_knee = zeros(numColPoints,1);
aD_ankle = zeros(numColPoints,1);

torqueParams = auxdata.active.params;
peakTorque = parameter(1, torqueParams.peak_torque.idx);
peakTime = parameter(1, torqueParams.peak_time.idx);
riseTime = parameter(1, torqueParams.rise_time.idx);
fallTime = parameter(1, torqueParams.fall_time.idx);

peakTorque = 0.5*(torqueParams.peak_torque.upper-torqueParams.peak_torque.lower)*(peakTorque+1) + torqueParams.peak_torque.lower;
peakTime = 0.5*(torqueParams.peak_time.upper-torqueParams.peak_time.lower)*(peakTime+1) + torqueParams.peak_time.lower;
riseTime = 0.5*(torqueParams.rise_time.upper-torqueParams.rise_time.lower)*(riseTime+1) + torqueParams.rise_time.lower;
fallTime = 0.5*(torqueParams.fall_time.upper-torqueParams.fall_time.lower)*(fallTime+1) + torqueParams.fall_time.lower;
aD = getTorqueControlFromParameters(peakTorque, peakTime, riseTime, fallTime, numColPoints);

if auxdata.active.hip
    signMoment_hip = auxdata.signMoment.hip;
    if auxdata.numActiveDOFs > 1
        % TODO
    else
        aD_hip = aD;
    end
end
if auxdata.active.knee
    signMoment_knee = auxdata.signMoment.knee;
    if auxdata.numActiveDOFs > 1
        % TODO
    else
        aD_knee = aD;
    end
end
if auxdata.active.ankle
    signMoment_ankle = auxdata.signMoment.ankle;
    if auxdata.numActiveDOFs > 1
        % TODO
    else
        aD_ankle = aD;
    end
end

% Exosuit torques
Texo_act_hip = auxdata.Tmax_act.*aD_hip.*signMoment_hip;
Texo_act_knee = auxdata.Tmax_act.*aD_knee.*signMoment_knee;
Texo_act_ankle = auxdata.Tmax_act.*aD_ankle.*signMoment_ankle;

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