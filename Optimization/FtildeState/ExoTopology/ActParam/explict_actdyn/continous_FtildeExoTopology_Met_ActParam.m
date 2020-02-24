function phaseout = continous_FtildeExoTopology_Met_ActParam(input)

% Get input data
auxdata         = input.auxdata;
NMuscles        = input.auxdata.NMuscles;
Ndof            = input.auxdata.Ndof;
tauAct          = input.auxdata.tauAct;
tauDeact        = input.auxdata.tauDeact;
params          = input.auxdata.params;
metabolicParams = input.auxdata.metabolicParams;
splinestruct    = input.auxdata.splinestruct;
speed           = input.auxdata.speed;
numColPoints    = size(input.phase.state,1);

% Get controls
e       = input.phase.control(:,1:NMuscles);
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
parameter = input.phase.parameter;
paramsLower = auxdata.paramsLower;
paramsUpper = auxdata.paramsUpper;
parameter = (paramsUpper-paramsLower).*parameter + 0.5*(paramsLower+paramsUpper);
torqueParams = auxdata.active.params;
peakTorque = parameter(1, torqueParams.peak_torque.idx);
peakTime = parameter(1, torqueParams.peak_time.idx);
riseTime = parameter(1, torqueParams.rise_time.idx);
fallTime = parameter(1, torqueParams.fall_time.idx);
paramControl = getTorqueControlFromParameters(peakTorque, peakTime, riseTime, ...
        fallTime, numColPoints);
    
if contains(auxdata.mod_name,'actHf_') || contains(auxdata.mod_name,'actHe_') 
    aD_hip = paramControl;

elseif contains(auxdata.mod_name,'actKf_') || contains(auxdata.mod_name,'actKe_')
    aD_knee = paramControl;        

elseif contains(auxdata.mod_name,'actAp_') || contains(auxdata.mod_name,'actAd_')
    aD_ankle = paramControl;        

elseif strcmp(auxdata.mod_name,'fitreopt_zhang2017_actHfAp')
    anklePeakTorque = parameter(1, torqueParams.ankle_peak_torque.idx);
    ankleParamControl = getTorqueControlFromParameters(anklePeakTorque, ...
        peakTime, riseTime, fallTime, numColPoints);
    
    aD_hip = paramControl;
    aD_ankle = ankleParamControl;
    
elseif strcmp(auxdata.mod_name,'fitreopt_zhang2017_actHfKf')
    kneePeakTorque = parameter(1, torqueParams.knee_peak_torque.idx);
    kneeParamControl = getTorqueControlFromParameters(kneePeakTorque, ...
        peakTime, riseTime, fallTime, numColPoints);
    
    aD_hip = paramControl;
    aD_knee = kneeParamControl;

elseif strcmp(auxdata.mod_name,'fitreopt_zhang2017_actKfAp')
    anklePeakTorque = parameter(1, torqueParams.ankle_peak_torque.idx);
    ankleParamControl = getTorqueControlFromParameters(anklePeakTorque, ...
        peakTime, riseTime, fallTime, numColPoints);
    
    aD_knee = paramControl;
    aD_ankle = ankleParamControl;
    
elseif strcmp(auxdata.mod_name,'fitreopt_zhang2017_actHfKfAp')
    kneePeakTorque = parameter(1, torqueParams.knee_peak_torque.idx);
    anklePeakTorque = parameter(1, torqueParams.ankle_peak_torque.idx);
    kneeParamControl = getTorqueControlFromParameters(kneePeakTorque, ...
        peakTime, riseTime, fallTime, numColPoints);
    ankleParamControl = getTorqueControlFromParameters(anklePeakTorque, ...
        peakTime, riseTime, fallTime, numColPoints);
    
    aD_hip = paramControl;
    aD_knee = kneeParamControl;
    aD_ankle = ankleParamControl;
    
elseif strcmp(auxdata.mod_name,'fitreopt_zhang2017_actHeKe')
    kneePeakTorque = parameter(1, torqueParams.knee_peak_torque.idx);
    kneeParamControl = getTorqueControlFromParameters(kneePeakTorque, ...
        peakTime, riseTime, fallTime, numColPoints);
    
    aD_hip = paramControl;
    aD_knee = kneeParamControl;
            
elseif strcmp(auxdata.mod_name,'fitreopt_zhang2017_actHfKfAp_multControls')
    peakTorque2 = parameter(1, torqueParams.peak_torque_2.idx);
    peakTime2 = parameter(1, torqueParams.peak_time_2.idx);
    riseTime2 = parameter(1, torqueParams.rise_time_2.idx);
    fallTime2 = parameter(1, torqueParams.fall_time_2.idx);
    
    kneeParamControl = getTorqueControlFromParameters(peakTorque2, peakTime2, riseTime2, ... 
            fallTime2, numColPoints);
    
    peakTorque3 = parameter(1, torqueParams.peak_torque_3.idx);
    peakTime3 = parameter(1, torqueParams.peak_time_3.idx);
    riseTime3 = parameter(1, torqueParams.rise_time_3.idx);
    fallTime3 = parameter(1, torqueParams.fall_time_3.idx);
    
    ankleParamControl = getTorqueControlFromParameters(peakTorque3, peakTime3, riseTime3, ... 
            fallTime3, numColPoints);
        
    aD_hip = paramControl;
    aD_knee = kneeParamControl;
    aD_ankle = ankleParamControl;
        
elseif strcmp(auxdata.mod_name,'fitreopt_zhang2017_actHfKf_multControls')
    peakTorque2 = parameter(1, torqueParams.peak_torque_2.idx);
    peakTime2 = parameter(1, torqueParams.peak_time_2.idx);
    riseTime2 = parameter(1, torqueParams.rise_time_2.idx);
    fallTime2 = parameter(1, torqueParams.fall_time_2.idx);
    
    kneeParamControl = getTorqueControlFromParameters(peakTorque2, peakTime2, riseTime2, ... 
            fallTime2, numColPoints);
        
    aD_hip = paramControl;
    aD_knee = kneeParamControl;
    
elseif strcmp(auxdata.mod_name,'fitreopt_zhang2017_actHfAp_multControls')
    peakTorque2 = parameter(1, torqueParams.peak_torque_2.idx);
    peakTime2 = parameter(1, torqueParams.peak_time_2.idx);
    riseTime2 = parameter(1, torqueParams.rise_time_2.idx);
    fallTime2 = parameter(1, torqueParams.fall_time_2.idx);
    
    ankleParamControl = getTorqueControlFromParameters(peakTorque2, peakTime2, riseTime2, ... 
            fallTime2, numColPoints);
        
    aD_hip = paramControl;
    aD_ankle = ankleParamControl;
        
elseif strcmp(auxdata.mod_name,'fitreopt_zhang2017_actKfAp_multControls')
    peakTorque2 = parameter(1, torqueParams.peak_torque_2.idx);
    peakTime2 = parameter(1, torqueParams.peak_time_2.idx);
    riseTime2 = parameter(1, torqueParams.rise_time_2.idx);
    fallTime2 = parameter(1, torqueParams.fall_time_2.idx);
    
    ankleParamControl = getTorqueControlFromParameters(peakTorque2, peakTime2, riseTime2, ... 
            fallTime2, numColPoints);
        
    aD_knee = paramControl;
    aD_ankle = ankleParamControl;
        
elseif strcmp(auxdata.mod_name,'fitreopt_zhang2017_actHeKe_multControls')
    peakTorque2 = parameter(1, torqueParams.peak_torque_2.idx);
    peakTime2 = parameter(1, torqueParams.peak_time_2.idx);
    riseTime2 = parameter(1, torqueParams.rise_time_2.idx);
    fallTime2 = parameter(1, torqueParams.fall_time_2.idx);
    
    kneeParamControl = getTorqueControlFromParameters(peakTorque2, peakTime2, riseTime2, ... 
            fallTime2, numColPoints);
        
    aD_hip = paramControl;
    aD_knee = kneeParamControl;
end

% PATH CONSTRAINTS
% Hill-equilibrium constraint
muscleData = DeGroote2016Muscle_FtildeState(a,Ftilde,dFtilde,splinestruct.LMT, ... 
        splinestruct.VMT,params,input.auxdata.Fvparam,input.auxdata.Fpparam, ...
        input.auxdata.Faparam);

% Exosuit torques
Texo_act_hip = auxdata.Tmax_act_hip.*aD_hip.*auxdata.signMoment.hip;
Texo_act_knee = auxdata.Tmax_act_knee.*aD_knee.*auxdata.signMoment.knee;
Texo_act_ankle = auxdata.Tmax_act_ankle.*aD_ankle.*auxdata.signMoment.ankle;

% Moments constraint
Topt = 150; % Reserve actuator strength [N/m] 
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
phaseout.path = [Tdiff muscleData.err];

% DYNAMIC CONSTRAINTS
% Activation dynamics
dadt = ones(numColPoints,NMuscles);
for m = 1:NMuscles
    dadt(:,m) = ActivationDynamics(e(:,m),a(:,m),tauAct(m),tauDeact(m),input.auxdata.b);
end

% Contraction dynamics is implicit
phaseout.dynamics = [dadt dFtilde];

% OBJECTIVE FUNCTION
% Calculate metabolic rate from Minetti & Alexander (1997) model
% vmax = params(5,:);  
% Fo = params(1,:);   
% Edot = zeros(numColPoints,NMuscles);
% for m = 1:NMuscles
%     v = vmax(1,m)*muscleData.vMtilde(:,m);
%     mass = metabolicParams(4,m);
%     Edot(:,m) = calcMinettiAlexanderProbe(v,vmax(1,m),Fo(1,m),a(:,m)); % / mass;
% end

% Calculate metabolic rate from Koelewijn et al. 2018 smooth approximation of
% Umberger model
Edot = zeros(numColPoints, NMuscles);
for m = 1:NMuscles
    Lce = muscleData.lMtilde(:,m)*params(2,m);
    Vce = muscleData.vMtilde(:,m)*params(5,m);
    mass = metabolicParams(4,m);
    paramStruct = [metabolicParams(1,m), metabolicParams(2,m), ...
                   metabolicParams(3,m), mass, ...
                   metabolicParams(5,m)];
               
    Edot(:,m) = mass * calcUmbergerKoelewijn2018Cost(Lce, Vce, ...
            muscleData.Fce(:,m), muscleData.FMltilde(:,m), e(:,m), a(:,m), ...
            paramStruct);
end

w_Edot = 1/input.auxdata.model_mass;
w_Res = 1e3;
w_dFtilde = 1;
w_e = 1e-3;
phaseout.integrand = w_Edot*sum(Edot, 2) + w_Res*sum(aT.^2,2) + ...
                     w_dFtilde*sum((dFtilde/10).^2,2) + w_e*sum(e.^2,2); 


