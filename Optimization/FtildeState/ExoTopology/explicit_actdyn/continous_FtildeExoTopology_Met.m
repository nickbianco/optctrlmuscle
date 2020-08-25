function phaseout = continous_FtildeExoTopology_Met(input)

% Get input data
auxdata         = input.auxdata;
NMuscles        = auxdata.NMuscles;
Ndof            = auxdata.Ndof;
tauAct          = auxdata.tauAct;
tauDeact        = auxdata.tauDeact;
params          = auxdata.params;
metabolicParams = auxdata.metabolicParams;
splinestruct    = auxdata.splinestruct;
speed           = input.auxdata.speed;
numColPoints    = size(input.phase.state,1);
NStates         = size(input.phase.state,2);
NControls       = size(input.phase.state,2);
T_exp           = splinestruct.ID;

% Get controls
e       = input.phase.control(:,1:NMuscles);
aT      = input.phase.control(:,NMuscles+1:NMuscles+Ndof);
dFtilde = 10*input.phase.control(:,NMuscles+Ndof+1:NMuscles+Ndof+NMuscles);
aD      = input.phase.control(:,end-(auxdata.numDeviceDOFs-1):end);

% Get states
a      = input.phase.state(:,1:NMuscles);
Ftilde = input.phase.state(:,NMuscles+1:NMuscles+NMuscles);


% Convert parameters to the correct range
parameter = input.phase.parameter;
paramsLower = auxdata.paramsLower;
paramsUpper = auxdata.paramsUpper;
parameter = (paramsUpper-paramsLower).*parameter + 0.5*(paramsLower+paramsUpper);

% Get moment arms and DOF controls
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
end

% PATH CONSTRAINTS
% Hill-equilibrium constraint
muscleData = DeGroote2016Muscle_FtildeState(a,Ftilde,dFtilde, ... 
    splinestruct.LMT,splinestruct.VMT,params, ... 
    auxdata.Fvparam,auxdata.Fpparam,auxdata.Faparam);

% Exosuit torques
Texo_act_hip = auxdata.Tmax_act_hip.*aD_hip.*exoMomentArms(:,1);
Texo_act_knee = auxdata.Tmax_act_knee.*aD_knee.*exoMomentArms(:,2).*auxdata.kneeAngleSign;
Texo_act_ankle = auxdata.Tmax_act_ankle.*aD_ankle.*exoMomentArms(:,3);

% Moments constraint
Topt = 150;
Tdiff = zeros(numColPoints,Ndof);
for dof = 1:Ndof
    index_sel=(dof-1)*(NMuscles)+1:(dof-1)*(NMuscles)+NMuscles;
    T_sim=sum(muscleData.FT.*splinestruct.MA(:,index_sel),2) + Topt*aT(:,dof);

    if dof==auxdata.hip_DOF
        T_sim = T_sim + Texo_act_hip;
    end
    if dof==auxdata.knee_DOF
        T_sim = T_sim + Texo_act_knee;
    end
    if dof==auxdata.ankle_DOF
        T_sim = T_sim + Texo_act_ankle;
    end
   
    Tdiff(:,dof) = (T_exp(:,dof)-T_sim);
end

phaseout.path = [Tdiff muscleData.err];

% DYNAMIC CONSTRAINTS
% Activation dynamics
dadt = ones(numColPoints,NMuscles);
daDdt = ones(numColPoints,auxdata.numDeviceDOFs);
for m = 1:NMuscles
    dadt(:,m) = ActivationDynamics(e(:,m),a(:,m),tauAct(m),tauDeact(m),input.auxdata.b);
end

% Contraction dynamics is implicit
phaseout.dynamics = [dadt dFtilde];

% OBJECTIVE FUNCTION
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

w_Edot = 1/(NMuscles*input.auxdata.model_mass);
w_Res = 1e3/Ndof;
% w_Control = 1/(3*NMuscles);
w_Reg = auxdata.regularizationWeight;

% dFtilde_diff = [zeros(1, NMuscles); dFtilde(2:end,:) - dFtilde(1:end-1,:)];
% e_diff = [zeros(1, NMuscles); e(2:end,:) - e(1:end-1,:)];
% Ftilde_diff = [zeros(1, NMuscles); Ftilde(2:end,:) - Ftilde(1:end-1,:)];
% a_diff = [zeros(1, NMuscles); a(2:end,:) - a(1:end-1,:)];
% aT_diff = [zeros(1, Ndof); aT(2:end,:) - aT(1:end-1,:)];
% aD_diff = [zeros(1, auxdata.numDeviceDOFs); aD(2:end,:) - aD(1:end-1,:)];

phaseout.integrand = w_Edot*sum(Edot, 2) + w_Res*sum(aT.^2,2) + ...
                     w_Reg*(sum((dFtilde/10).^2,2)/NMuscles + sum(e.^2,2)/NMuscles + ...
                     sum(a.^2,2)/NMuscles); 
                        