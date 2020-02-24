function phaseout = continous_Ftilde_vAExoTopology_Met_ActParam(input)

% Get input data
auxdata         = input.auxdata;
NMuscles        = input.auxdata.NMuscles;
Ndof            = input.auxdata.Ndof;
tauAct          = input.auxdata.tauAct;
tauDeact        = input.auxdata.tauDeact;
params          = input.auxdata.params;
% metabolicParams = input.auxdata.metabolicParams;
splinestruct    = input.auxdata.splinestruct;
numColPoints    = size(input.phase.state,1);

% Get controls
vA      = 100*input.phase.control(:,1:NMuscles);
aT      = input.phase.control(:,NMuscles+1:NMuscles+Ndof);
dFtilde = 10*input.phase.control(:,NMuscles+Ndof+1:NMuscles+Ndof+NMuscles);

% Get states
a      = input.phase.state(:,1:NMuscles);
Ftilde = input.phase.state(:,NMuscles+1:NMuscles+NMuscles);

% Initialize torque arrays
aD_hip = zeros(numColPoints,1);
aD_knee = zeros(numColPoints,1);
aD_ankle = zeros(numColPoints,1);

% Convert parameters to the correct range
paramsLower = auxdata.paramsLower;
paramsUpper = auxdata.paramsUpper;
parameter = 0.5*(paramsUpper-paramsLower).*(input.phase.parameter+1) + paramsLower;
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

% PATH CONSTRAINTS
% Activation dynamics - De Groote et al. (2009)
act1 = vA + a./(ones(size(a,1),1)*tauDeact);
act2 = vA + a./(ones(size(a,1),1)*tauAct);

% Hill-equilibrium constraint
muscleData = DeGroote2016Muscle_FtildeState(a,Ftilde,dFtilde,splinestruct.LMT,splinestruct.VMT,params,input.auxdata.Fvparam,input.auxdata.Fpparam,input.auxdata.Faparam);

% Exosuit torques
Texo_act_hip = auxdata.Tmax_act.*aD_hip.*auxdata.signMoment.hip;
Texo_act_knee = auxdata.Tmax_act.*aD_knee.*auxdata.signMoment.knee;
Texo_act_ankle = auxdata.Tmax_act.*aD_ankle.*auxdata.signMoment.ankle;

% Moments constraint
Topt = 150;
Tdiff = zeros(numColPoints,Ndof);
for dof = 1:Ndof
    T_exp=splinestruct.ID(:,dof);
    index_sel=(dof-1)*(NMuscles)+1:(dof-1)*(NMuscles)+NMuscles;
    T_sim=sum(muscleData.FT.*splinestruct.MA(:,index_sel),2) + Topt*aT(:,dof);

    if dof==input.auxdata.hip_DOF
        T_sim = T_sim + Texo_act_hip;
    end
    if dof==input.auxdata.knee_DOF
        T_sim = T_sim + Texo_act_knee;
    end
    if dof==input.auxdata.ankle_DOF
        T_sim = T_sim + Texo_act_ankle;
    end
   
    Tdiff(:,dof) = (T_exp-T_sim);
end

% PATH CONSTRAINTS
phaseout.path = [Tdiff muscleData.err act1 act2];

% DYNAMIC CONSTRAINTS
% Activation dynamics is implicit
% Contraction dynamics is implicit
phaseout.dynamics = [vA dFtilde];

% OBJECTIVE FUNCTION
% Calculate metabolic rate from Minetti & Alexander (1997) model
vmax = params(5,:);  
Fo = params(1,:);   
Edot = zeros(numColPoints,NMuscles);
for m = 1:NMuscles
    v = vmax(1,m)*muscleData.vMtilde(:,m);
    Edot(:,m) = calcMinettiAlexanderProbe(v,vmax(1,m),Fo(1,m),a(:,m));
end

% Calculate metabolic rate from Umberger(2010) model
% Edot = zeros(numColPoints, NMuscles);
% for m = 1:NMuscles
%     Lce = muscleData.lMtilde(:,m)*params(2,m);
%     Vce = muscleData.vMtilde(:,m)*params(5,m);
%     paramStruct = [metabolicParams(1,m), metabolicParams(2,m), ...
%                    metabolicParams(3,m), metabolicParams(4,m), ...
%                    metabolicParams(5,m)];
%     Edot(:,m) = calcUmbergerCost2010(max(0, Lce), ...
%                                      Vce, ...
%                                      max(0, muscleData.Fce(:,m)), ...
%                                      max(0, muscleData.FMltilde(:,m)), ...
%                                      min(max(0, e(:,m)), 1), ... 
%                                      min(max(0, a(:,m)), 1), ...
%                                      paramStruct);
% end

% Overwriting first time point of metabolics to avoid effects that initial spikes
% in fiber powers may have on the cost.
% Edot(1,:) = Edot(2,:);

w_Edot = 1/input.auxdata.model_mass;
w_Res = 1000;
w_Reg = 1e-6;
phaseout.integrand = w_Edot*sum(Edot, 2) + w_Res.*sum(aT.^2,2) + ...
                     w_Reg.*(sum(dFtilde.^2,2) + sum(e.^2,2));


