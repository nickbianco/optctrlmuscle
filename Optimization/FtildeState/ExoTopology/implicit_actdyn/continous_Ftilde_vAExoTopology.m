function phaseout = continous_Ftilde_vAExoTopology(input)

% Get input data
NMuscles        = input.auxdata.NMuscles;
Ndof            = input.auxdata.Ndof;
tauAct          = input.auxdata.tauAct;
tauDeact        = input.auxdata.tauDeact;
params          = input.auxdata.params;
splinestruct    = input.auxdata.splinestruct;
numColPoints    = size(input.phase.state,1);

% Get controls
vA      = 100*input.phase.control(:,1:NMuscles);
aT      = input.phase.control(:,NMuscles+1:NMuscles+Ndof);
dFtilde = 10*input.phase.control(:,NMuscles+Ndof+1:NMuscles+Ndof+NMuscles);
aD      = input.phase.control(:,end-(input.auxdata.numDeviceDOFs-1):end);

% Get states
a      = input.phase.state(:,1:NMuscles);
Ftilde = input.phase.state(:,NMuscles+1:NMuscles+NMuscles);

% Get moment arms and DOF controls
exoMomentArms = zeros(numColPoints,3);
aD_hip = zeros(numColPoints,1);
aD_knee = zeros(numColPoints,1);
aD_ankle = zeros(numColPoints,1);
if input.auxdata.device.hip
    exoMomentArms(:,1) = input.phase.parameter(:,input.auxdata.device.hip);
    if input.auxdata.numDeviceDOFs > 1
        aD_hip = aD(:,input.auxdata.device.hip);
    else
        aD_hip = aD;
    end
end
if input.auxdata.device.knee
    exoMomentArms(:,2) = input.phase.parameter(:,input.auxdata.device.knee);
    if input.auxdata.numDeviceDOFs > 1
        aD_knee = aD(:,input.auxdata.device.knee);
    else
        aD_knee = aD;
    end
end
if input.auxdata.device.ankle
    exoMomentArms(:,3) = input.phase.parameter(:,input.auxdata.device.ankle);
    if input.auxdata.numDeviceDOFs > 1
        aD_ankle = aD(:,input.auxdata.device.ankle);
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
Texo_act_hip = input.auxdata.Tmax_act.*aD_hip.*exoMomentArms(:,1);
Texo_act_knee = input.auxdata.Tmax_act.*aD_knee.*exoMomentArms(:,2).*input.auxdata.kneeAngleSign;
Texo_act_ankle = input.auxdata.Tmax_act.*aD_ankle.*exoMomentArms(:,3);

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

phaseout.path = [Tdiff muscleData.err act1 act2];

% DYNAMIC CONSTRAINTS
% Activation dynamics is implicit
% Contraction dynamics is implicit
phaseout.dynamics = [vA dFtilde];

% OBJECTIVE FUNCTION
w1 = 1000;
w2 = 0.01;
phaseout.integrand = sum(a.^2,2)+ w1.*sum(aT.^2,2)+ w2*sum((vA/100).^2,2);


