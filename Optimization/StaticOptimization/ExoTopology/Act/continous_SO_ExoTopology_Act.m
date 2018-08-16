function phaseout = continous_SO_ExoTopology_Act(input)

% Get input data
NMuscles        = input.auxdata.NMuscles;
Ndof            = input.auxdata.Ndof;
tauAct          = input.auxdata.tauAct;
tauDeact        = input.auxdata.tauDeact;
params          = input.auxdata.params;
splinestruct    = input.auxdata.splinestruct;
numColPoints    = size(input.phase.state,1);

% Get controls
e  = input.phase.control(:,1:NMuscles);
aT = input.phase.control(:,NMuscles+1:NMuscles+Ndof);
aD = input.phase.control(:,end-(input.auxdata.numActiveDOFs-1):end);

% Get states
a = input.phase.state(:, 1:NMuscles);

% Convert parameters to the correct range
paramsLower = input.auxdata.paramsLower;
paramsUpper = input.auxdata.paramsUpper;
parameter = 0.5*(paramsUpper-paramsLower).*(input.phase.parameter+1) + paramsLower;

% Get moment arms and DOF controls
exoMomentArms = zeros(numColPoints,3);
aD_hip = zeros(numColPoints,1);
aD_knee = zeros(numColPoints,1);
aD_ankle = zeros(numColPoints,1);
if input.auxdata.active.hip
    exoMomentArms(:,1) = parameter(:,input.auxdata.active.hip);
    if input.auxdata.numActiveDOFs > 1
        aD_hip = aD(:,input.auxdata.active.hip);
    else
        aD_hip = aD;
    end
end
if input.auxdata.active.knee
    exoMomentArms(:,2) = parameter(:,input.auxdata.active.knee);
    if input.auxdata.numActiveDOFs > 1
        aD_knee = aD(:,input.auxdata.active.knee);
    else
        aD_knee = aD;
    end
end
if input.auxdata.active.ankle
    exoMomentArms(:,3) = parameter(:,input.auxdata.active.ankle);
    if input.auxdata.numActiveDOFs > 1
        aD_ankle = aD(:,input.auxdata.active.ankle);
    else
        aD_ankle = aD;
    end
end

% PATH CONSTRAINTS

% Get muscle forces
[F, ~, ~, ~, ~, ~] = HillModel_RigidTendon(a, splinestruct.LMT, ... 
                                              splinestruct.VMT, ... 
                                              params, ... 
                                              input.auxdata.Fvparam, ...
                                              input.auxdata.Fpparam, ...
                                              input.auxdata.Faparam);
% Exosuit torques
% Active device
Texo_act_hip = input.auxdata.Tmax_act*aD_hip.*exoMomentArms(:,1);
Texo_act_knee = input.auxdata.Tmax_act*aD_knee.*exoMomentArms(:,2).*input.auxdata.kneeAngleSign;
Texo_act_ankle = input.auxdata.Tmax_act*aD_ankle.*exoMomentArms(:,3);

% Moments constraint
Topt = 150;
Tdiff = zeros(numColPoints,Ndof);
for dof = 1:Ndof
    T_exp=splinestruct.ID(:,dof);
    index_sel=(dof-1)*(NMuscles)+1:(dof-1)*(NMuscles)+NMuscles;
    T_sim=sum(F.*splinestruct.MA(:,index_sel),2) + Topt*aT(:,dof);

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

phaseout.path = Tdiff;

% DYNAMIC CONSTRAINTS
% Activation dynamics
dadt = ones(numColPoints,NMuscles);
for m = 1:NMuscles
    dadt(:,m) = ActivationDynamics(e(:,m),a(:,m),tauAct(m),tauDeact(m),input.auxdata.b);
end

phaseout.dynamics = dadt;

% OBJECTIVE FUNCTION
w1 = 1000;
phaseout.integrand = sum(e.^2 + a.^2, 2) + w1.*sum(aT.^2,2);


