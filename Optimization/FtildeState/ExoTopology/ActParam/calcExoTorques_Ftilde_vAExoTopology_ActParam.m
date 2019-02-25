function [ExoTorques_Act] = calcExoTorques_Ftilde_vAExoTopology_ActParam(OptInfo,DatStore)

time = OptInfo.result.solution.phase.time;
numColPoints = length(time);
auxdata = OptInfo.result.setup.auxdata;
Ndof = auxdata.Ndof;

% Get parameters
parameter = OptInfo.result.solution.parameter;

% Initialize torque arrays
aD_hip = zeros(numColPoints,1);
aD_knee = zeros(numColPoints,1);
aD_ankle = zeros(numColPoints,1);

% Convert parameters to the correct range
paramsLower = auxdata.paramsLower;
paramsUpper = auxdata.paramsUpper;
parameter = 0.5*(paramsUpper-paramsLower).*(parameter+1) + paramsLower;
torqueParams = auxdata.active.params;
peakTorque = parameter(1, torqueParams.peak_torque.idx);
peakTime = parameter(1, torqueParams.peak_time.idx);
riseTime = parameter(1, torqueParams.rise_time.idx);
fallTime = parameter(1, torqueParams.fall_time.idx);
aD = getTorqueControlFromParameters(peakTorque, peakTime, riseTime, fallTime, numColPoints);

if contains(auxdata.mod_name,'actHf_') || contains(auxdata.mod_name,'actHe_') 
    aD_hip = aD;

elseif contains(auxdata.mod_name,'actKf_') || contains(auxdata.mod_name,'actKe_')
    aD_knee = aD;

elseif contains(auxdata.mod_name,'actAp_') || contains(auxdata.mod_name,'actAd_')
    aD_ankle = aD;

elseif strcmp(auxdata.mod_name,'fitreopt_zhang2017_actHfAp')
    aD_hip = aD;
    ankleTorqueScale = parameter(1, torqueParams.ankle_torque_scale.idx);
    aD_ankle = getTorqueControlFromParameters(peakTorque*ankleTorqueScale, ...
        peakTime, riseTime, fallTime, numColPoints);   
    
elseif strcmp(auxdata.mod_name,'fitreopt_zhang2017_actHfKf')
    aD_hip = aD;
    kneeTorqueScale = parameter(1, torqueParams.knee_torque_scale.idx);
    aD_knee = getTorqueControlFromParameters(peakTorque*kneeTorqueScale, ...
        peakTime, riseTime, fallTime, numColPoints);

elseif strcmp(auxdata.mod_name,'fitreopt_zhang2017_actKfAp')
    aD_knee = aD;
    ankleTorqueScale = parameter(1, torqueParams.ankle_torque_scale.idx);
    aD_ankle = getTorqueControlFromParameters(peakTorque*ankleTorqueScale, ...
        peakTime, riseTime, fallTime, numColPoints);
    
elseif strcmp(auxdata.mod_name,'fitreopt_zhang2017_actHfKfAp')
    aD_hip = aD;
    kneeTorqueScale = parameter(1, torqueParams.knee_torque_scale.idx);
    ankleTorqueScale = parameter(1, torqueParams.ankle_torque_scale.idx);
    aD_knee = getTorqueControlFromParameters(peakTorque*kneeTorqueScale, ...
        peakTime, riseTime, fallTime, numColPoints);
    aD_ankle = getTorqueControlFromParameters(peakTorque*ankleTorqueScale, ...
        peakTime, riseTime, fallTime, numColPoints);
    
elseif strcmp(auxdata.mod_name,'fitreopt_zhang2017_actHeKe')
    aD_hip = aD;
    kneeTorqueScale = parameter(1, torqueParams.knee_torque_scale.idx);
    aD_knee = getTorqueControlFromParameters(peakTorque*kneeTorqueScale, ...
        peakTime, riseTime, fallTime, numColPoints);   
end    

% Exosuit torques
Texo_act_hip = auxdata.Tmax_act.*aD_hip.*auxdata.signMoment.hip;
Texo_act_knee = auxdata.Tmax_act.*aD_knee.*auxdata.signMoment.knee;
Texo_act_ankle = auxdata.Tmax_act.*aD_ankle.*auxdata.signMoment.ankle;

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