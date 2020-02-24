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
parameter = (paramsUpper-paramsLower).*parameter + 0.5*(paramsLower+paramsUpper);
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
    anklePeakTorque = parameter(1, torqueParams.ankle_peak_torque.idx);
    aD_ankle = getTorqueControlFromParameters(anklePeakTorque, ...
        peakTime, riseTime, fallTime, numColPoints);   
    
elseif strcmp(auxdata.mod_name,'fitreopt_zhang2017_actHfKf')
    aD_hip = aD;
    kneePeakTorque = parameter(1, torqueParams.knee_peak_torque.idx);
    aD_knee = getTorqueControlFromParameters(kneePeakTorque, ...
        peakTime, riseTime, fallTime, numColPoints);

elseif strcmp(auxdata.mod_name,'fitreopt_zhang2017_actKfAp')
    aD_knee = aD;
    anklePeakTorque = parameter(1, torqueParams.ankle_peak_torque.idx);
    aD_ankle = getTorqueControlFromParameters(anklePeakTorque, ...
        peakTime, riseTime, fallTime, numColPoints);
    
elseif strcmp(auxdata.mod_name,'fitreopt_zhang2017_actHfKfAp')
    aD_hip = aD;
    kneePeakTorque = parameter(1, torqueParams.knee_peak_torque.idx);
    anklePeakTorque = parameter(1, torqueParams.ankle_peak_torque.idx);
    aD_knee = getTorqueControlFromParameters(kneePeakTorque, ...
        peakTime, riseTime, fallTime, numColPoints);
    aD_ankle = getTorqueControlFromParameters(anklePeakTorque, ...
        peakTime, riseTime, fallTime, numColPoints);
    
elseif strcmp(auxdata.mod_name,'fitreopt_zhang2017_actHeKe')
    aD_hip = aD;
    kneePeakTorque = parameter(1, torqueParams.knee_peak_torque.idx);
    aD_knee = getTorqueControlFromParameters(kneePeakTorque, ...
        peakTime, riseTime, fallTime, numColPoints);
    
elseif strcmp(auxdata.mod_name,'fitreopt_zhang2017_actHfKfAp_multControls')
    aD_hip = aD;
    
    peakTorque2 = parameter(1, torqueParams.peak_torque_2.idx);
    peakTime2 = parameter(1, torqueParams.peak_time_2.idx);
    riseTime2 = parameter(1, torqueParams.rise_time_2.idx);
    fallTime2 = parameter(1, torqueParams.fall_time_2.idx);
    
    aD_knee = getTorqueControlFromParameters(peakTorque2, peakTime2, riseTime2, ... 
            fallTime2, numColPoints);
    
    peakTorque3 = parameter(1, torqueParams.peak_torque_3.idx);
    peakTime3 = parameter(1, torqueParams.peak_time_3.idx);
    riseTime3 = parameter(1, torqueParams.rise_time_3.idx);
    fallTime3 = parameter(1, torqueParams.fall_time_3.idx);
    
    aD_ankle = getTorqueControlFromParameters(peakTorque3, peakTime3, riseTime3, ... 
            fallTime3, numColPoints);

elseif strcmp(auxdata.mod_name,'fitreopt_zhang2017_actHfKf_multControls')
    aD_hip = aD;
    
    peakTorque2 = parameter(1, torqueParams.peak_torque_2.idx);
    peakTime2 = parameter(1, torqueParams.peak_time_2.idx);
    riseTime2 = parameter(1, torqueParams.rise_time_2.idx);
    fallTime2 = parameter(1, torqueParams.fall_time_2.idx);
    
    aD_knee = getTorqueControlFromParameters(peakTorque2, peakTime2, riseTime2, ... 
            fallTime2, numColPoints);
    
elseif strcmp(auxdata.mod_name,'fitreopt_zhang2017_actHfAp_multControls')
    aD_hip = aD;
    
    peakTorque2 = parameter(1, torqueParams.peak_torque_2.idx);
    peakTime2 = parameter(1, torqueParams.peak_time_2.idx);
    riseTime2 = parameter(1, torqueParams.rise_time_2.idx);
    fallTime2 = parameter(1, torqueParams.fall_time_2.idx);
    
    aD_ankle = getTorqueControlFromParameters(peakTorque2, peakTime2, riseTime2, ... 
            fallTime2, numColPoints);
        
elseif strcmp(auxdata.mod_name,'fitreopt_zhang2017_actKfAp_multControls')
    aD_knee = aD;
    
    peakTorque2 = parameter(1, torqueParams.peak_torque_2.idx);
    peakTime2 = parameter(1, torqueParams.peak_time_2.idx);
    riseTime2 = parameter(1, torqueParams.rise_time_2.idx);
    fallTime2 = parameter(1, torqueParams.fall_time_2.idx);
    
    aD_ankle = getTorqueControlFromParameters(peakTorque2, peakTime2, riseTime2, ... 
            fallTime2, numColPoints);
        
elseif strcmp(auxdata.mod_name,'fitreopt_zhang2017_actHeKe_multControls')
    aD_hip = aD;
    
    peakTorque2 = parameter(1, torqueParams.peak_torque_2.idx);
    peakTime2 = parameter(1, torqueParams.peak_time_2.idx);
    riseTime2 = parameter(1, torqueParams.rise_time_2.idx);
    fallTime2 = parameter(1, torqueParams.fall_time_2.idx);
    
    aD_knee = getTorqueControlFromParameters(peakTorque2, peakTime2, riseTime2, ... 
            fallTime2, numColPoints);
        
end

% Exosuit torques
Texo_act_hip = auxdata.Tmax_act_hip.*aD_hip.*auxdata.signMoment.hip;
Texo_act_knee = auxdata.Tmax_act_knee.*aD_knee.*auxdata.signMoment.knee;
Texo_act_ankle = auxdata.Tmax_act_ankle.*aD_ankle.*auxdata.signMoment.ankle;

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