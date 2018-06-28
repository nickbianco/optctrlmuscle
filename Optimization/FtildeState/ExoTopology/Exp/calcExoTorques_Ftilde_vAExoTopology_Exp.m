function [ExoTorques_Exp] = calcExoTorques_Ftilde_vAExoTopology_Exp(OptInfo,DatStore)

time = OptInfo.result.solution.phase.time;
numColPoints = length(time);
auxdata = OptInfo.result.setup.auxdata;
Ndof = auxdata.Ndof;
parameter = OptInfo.result.solution.parameter;

% Shift exoskeleton torque curve
torqueParams = auxdata.active.params;
shiftTime = parameter(1, torqueParams.shift_time.idx);
shiftTime = 0.5*(torqueParams.shift_time.upper-torqueParams.shift_time.lower)*(shiftTime+1) + torqueParams.shift_time.lower;
new_peak_time = auxdata.T_exo_peak_time + shiftTime;

aD_exo = zeros(length(time), length(auxdata.JointEXOSpline));
for i = 1:length(auxdata.JointEXOSpline)
    aD_exo(:,i) = ppval(auxdata.JointEXOSpline(i),time);
end
[~,idxs] = max(aD_exo);
currIdx = max(idxs);
[~,newIdx] = min(abs(time-new_peak_time));

shift_dir = sign(newIdx - currIdx);
aD_exo_shift = aD_exo;
if nnz(shift_dir)
    shift_idx = abs(newIdx - currIdx);
    if shift_dir > 0
        aD_exo_shift = [zeros(shift_idx,3); aD_exo(1:end-shift_idx,:)];
    elseif shift_dir <= 0
        aD_exo_shift = [aD_exo(shift_idx:end,:); zeros(shift_idx-1,3)];
    end
end

% Get moment arms
peakTorque = parameter(1, torqueParams.peak_torque.idx);
peakTorque = 0.5*(torqueParams.peak_torque.upper-torqueParams.peak_torque.lower)*(peakTorque+1) + torqueParams.peak_torque.lower;
signMoment_hip = 1;
signMoment_knee = 1;
signMoment_ankle = 1;
aD_hip = zeros(numColPoints,1);
aD_knee = zeros(numColPoints,1);
aD_ankle = zeros(numColPoints,1);
if auxdata.active.hip
    signMoment_hip = auxdata.signMoment.hip;
    if auxdata.numActiveDOFs > 1
        % TODO
    else
        aD_hip = peakTorque*aD_exo_shift(:,auxdata.hip_DOF);
    end
end
if auxdata.active.knee
    signMoment_knee = auxdata.signMoment.knee;
    if auxdata.numActiveDOFs > 1
        % TODO
    else
        aD_knee = peakTorque*aD_exo_shift(:,auxdata.knee_DOF);
    end
end
if auxdata.active.ankle
    signMoment_ankle = auxdata.signMoment.ankle;
    if auxdata.numActiveDOFs > 1
        % TODO
    else
        aD_ankle = peakTorque*aD_exo_shift(:,auxdata.ankle_DOF);
    end
end

% Exosuit torques
Texo_act_hip = auxdata.Tmax_act.*aD_hip.*signMoment_hip;
Texo_act_knee = auxdata.Tmax_act.*aD_knee.*signMoment_knee;
Texo_act_ankle = auxdata.Tmax_act.*aD_ankle.*signMoment_ankle;

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