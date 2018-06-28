function phaseout = continous_Ftilde_vAExoTopology_Met_ActParam(input)

% Get input data
NMuscles        = input.auxdata.NMuscles;
Ndof            = input.auxdata.Ndof;
tauAct          = input.auxdata.tauAct;
tauDeact        = input.auxdata.tauDeact;
params          = input.auxdata.params;
metabolicParams = input.auxdata.metabolicParams;
splinestruct    = input.auxdata.splinestruct;
numColPoints    = size(input.phase.state,1);

% Get controls
vA      = 100*input.phase.control(:,1:NMuscles);
aT      = input.phase.control(:,NMuscles+1:NMuscles+Ndof);
dFtilde = 10*input.phase.control(:,NMuscles+Ndof+1:NMuscles+Ndof+NMuscles);

% Get states
a      = input.phase.state(:,1:NMuscles);
Ftilde = input.phase.state(:,NMuscles+1:NMuscles+NMuscles);

% Get moment arms and DOF controls
signMoment_hip = 1;
signMoment_knee = 1;
signMoment_ankle = 1;
aD_hip = zeros(numColPoints,1);
aD_knee = zeros(numColPoints,1);
aD_ankle = zeros(numColPoints,1);

torqueParams = input.auxdata.active.params;
peakTorque = input.phase.parameter(1, torqueParams.peak_torque.idx);
peakTime = input.phase.parameter(1, torqueParams.peak_time.idx);
riseTime = input.phase.parameter(1, torqueParams.rise_time.idx);
fallTime = input.phase.parameter(1, torqueParams.fall_time.idx);

peakTorque = 0.5*(torqueParams.peak_torque.upper-torqueParams.peak_torque.lower)*(peakTorque+1) + torqueParams.peak_torque.lower;
peakTime = 0.5*(torqueParams.peak_time.upper-torqueParams.peak_time.lower)*(peakTime+1) + torqueParams.peak_time.lower;
riseTime = 0.5*(torqueParams.rise_time.upper-torqueParams.rise_time.lower)*(riseTime+1) + torqueParams.rise_time.lower;
fallTime = 0.5*(torqueParams.fall_time.upper-torqueParams.fall_time.lower)*(fallTime+1) + torqueParams.fall_time.lower;

aD = getTorqueControlFromParameters(peakTorque, peakTime, riseTime, fallTime, numColPoints);

if input.auxdata.active.hip
    signMoment_hip = input.auxdata.signMoment.hip;
    if input.auxdata.numActiveDOFs > 1
        % TODO
    else
        aD_hip = aD;
    end
end
if input.auxdata.active.knee
    signMoment_knee = input.auxdata.signMoment.knee;
    if input.auxdata.numActiveDOFs > 1
        % TODO
    else
        aD_knee = aD;
    end
end
if input.auxdata.active.ankle
    signMoment_ankle = input.auxdata.signMoment.ankle;
    if input.auxdata.numActiveDOFs > 1
        % TODO
    else
        aD_ankle = aD;
    end
end

% PATH CONSTRAINTS
% Activation dynamics - De Groote et al. (2009)
act1 = vA + a./(ones(size(a,1),1)*tauDeact);
act2 = vA + a./(ones(size(a,1),1)*tauAct);

% Hill-equilibrium constraint
muscleData = DeGroote2016Muscle_FtildeState(a,Ftilde,dFtilde,splinestruct.LMT,splinestruct.VMT,params,input.auxdata.Fvparam,input.auxdata.Fpparam,input.auxdata.Faparam);

% Exosuit torques
Texo_act_hip = input.auxdata.Tmax_act.*aD_hip.*signMoment_hip;
Texo_act_knee = input.auxdata.Tmax_act.*aD_knee.*signMoment_knee;
Texo_act_ankle = input.auxdata.Tmax_act.*aD_ankle.*signMoment_ankle;

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

if isfield(input.auxdata, 'powerMatchType')
    q = splinestruct.IK;
    dq = diff(q);
    t = input.phase.time;
    dt = repmat(diff(t), [1, size(q,2)]);
    dqdt = dq./dt;
        
    Texo = [Texo_act_hip Texo_act_knee Texo_act_ankle];
    Pexo = Texo(1:(end-1),:).*dqdt;
    
    switch input.auxdata.powerMatchType
        case 'avg_pos'
            Pexo_pos = Pexo;
            Pexo_pos(Pexo_pos<0) = 0;
            Pexo_match_val = sum(mean(Pexo_pos)) / input.auxdata.model_mass;
        case 'avg_net'
            Pexo_match_val = sum(mean(Pexo)) / input.auxdata.model_mass;
    end
   
    phaseout.path = [Tdiff muscleData.err act1 act2 Pexo_match_val*ones(numColPoints,1)];
else
    phaseout.path = [Tdiff muscleData.err act1 act2];
end

% DYNAMIC CONSTRAINTS
% Activation dynamics is implicit
% Contraction dynamics is implicit
phaseout.dynamics = [vA dFtilde];

% OBJECTIVE FUNCTION
Edot = zeros(numColPoints, NMuscles);
for m = 1:NMuscles
    Lce = muscleData.lMtilde(:,m)*params(2,m);
    Vce = muscleData.vMtilde(:,m)*params(5,m);
    u = computeExcitationRaasch(a(:,m), vA(:,m), tauDeact(m), tauAct(m));
    paramStruct = [metabolicParams(1,m), metabolicParams(2,m), ...
                   metabolicParams(3,m), metabolicParams(4,m), ...
                   metabolicParams(5,m)];
    Edot(:,m) = calcUmbergerCost2010(max(0, Lce), ...
                                     Vce, ...
                                     max(0, muscleData.Fce(:,m)), ...
                                     max(0, muscleData.FMltilde(:,m)), ...
                                     min(max(0, u), 1), ...
                                     min(max(0, a(:,m)), 1), ...
                                     paramStruct);
end
w_aT = 1000;
w_a = 0.05;
w_vA = 0.05;
w_Edot = 1/(input.auxdata.model_mass*9.81*1.25);
% phaseout.integrand = w_Edot*sum(Edot, 2) + w_aT.*sum(aT.^2,2) + w_a*sum(a.^2,2) + w_vA*sum((vA/100).^2,2);
phaseout.integrand = w_Edot*sum(Edot, 2) + w_aT.*sum(aT.^2,2) + sum((vA/100).^2,2);


