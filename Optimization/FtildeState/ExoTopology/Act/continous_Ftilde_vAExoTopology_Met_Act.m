function phaseout = continous_Ftilde_vAExoTopology_Met_Act(input)

% Get input data
auxdata         = input.auxdata;
NMuscles        = auxdata.NMuscles;
Ndof            = auxdata.Ndof;
tauAct          = auxdata.tauAct;
tauDeact        = auxdata.tauDeact;
params          = auxdata.params;
metabolicParams = auxdata.metabolicParams;
splinestruct    = auxdata.splinestruct;
numColPoints    = size(input.phase.state,1);
T_exp           = splinestruct.ID;

% Get controls
vA      = 100*input.phase.control(:,1:NMuscles);
aT      = input.phase.control(:,NMuscles+1:NMuscles+Ndof);
dFtilde = 10*input.phase.control(:,NMuscles+Ndof+1:NMuscles+Ndof+NMuscles);
aD      = input.phase.control(:,end-(auxdata.numActiveDOFs-1):end);

% Get states
a      = input.phase.state(:,1:NMuscles);
Ftilde = input.phase.state(:,NMuscles+1:NMuscles+NMuscles);
if auxdata.shift_exo_peaks
    aDshift = input.phase.state(:,end);
end

% Convert parameters to the correct range
paramsLower = auxdata.paramsLower;
paramsUpper = auxdata.paramsUpper;
parameter = 0.5*(paramsUpper-paramsLower).*(input.phase.parameter+1) + paramsLower;

% Get moment arms and DOF controls
exoMomentArms = zeros(numColPoints,3);
aD_hip = zeros(numColPoints,1);
aD_knee = zeros(numColPoints,1);
aD_ankle = zeros(numColPoints,1);
if auxdata.active.hip
    exoMomentArms(:,1) = parameter(:,auxdata.active.hip);
    if auxdata.numActiveDOFs > 1
        aD_hip = aD(:,auxdata.active.hip);
    else
        aD_hip = aD;
    end
end
if auxdata.active.knee
    exoMomentArms(:,2) = parameter(:,auxdata.active.knee);
    if auxdata.numActiveDOFs > 1
        aD_knee = aD(:,auxdata.active.knee);
    else
        aD_knee = aD;
    end
end
if auxdata.active.ankle
    exoMomentArms(:,3) = parameter(:,auxdata.active.ankle);
    if auxdata.numActiveDOFs > 1
        aD_ankle = aD(:,auxdata.active.ankle);
    else
        aD_ankle = aD;
    end
end

% PATH CONSTRAINTS
% Activation dynamics - De Groote et al. (2009)
act1 = vA + a./(ones(size(a,1),1)*tauDeact);
act2 = vA + a./(ones(size(a,1),1)*tauAct);

% Hill-equilibrium constraint
muscleData = DeGroote2016Muscle_FtildeState(a,Ftilde,dFtilde, ... 
    splinestruct.LMT,splinestruct.VMT,params, ... 
    auxdata.Fvparam,auxdata.Fpparam,auxdata.Faparam);

% Exosuit torques
Texo_act_hip = auxdata.Tmax_act.*aD_hip.*exoMomentArms(:,1);
Texo_act_knee = auxdata.Tmax_act.*aD_knee.*exoMomentArms(:,2).*auxdata.kneeAngleSign;
Texo_act_ankle = auxdata.Tmax_act.*aD_ankle.*exoMomentArms(:,3);

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

phaseout.path = [Tdiff muscleData.err act1 act2];

% DYNAMIC CONSTRAINTS
% Muscle activation dynamics is implicit
% Contraction dynamics is implicit
phaseout.dynamics = [vA dFtilde];

if auxdata.shift_exo_peaks
    % Device activation dynamics is explicit
    tauActExo = parameter(1,end);
    dadt = ActivationDynamics(aD, aDshift, tauActExo, 0.01, auxdata.b);
    phaseout.dynamics = [phaseout.dynamics dadt];
end

% OBJECTIVE FUNCTION
Edot = zeros(numColPoints, NMuscles);
for m = 1:NMuscles
    Lce = muscleData.lMtilde(:,m)*params(2,m);
    Vce = muscleData.vMtilde(:,m)*params(5,m);
%     u = computeExcitationRaasch(a(:,m), vA(:,m), tauDeact(m), tauAct(m));
    paramStruct = [metabolicParams(1,m), metabolicParams(2,m), ...
                   metabolicParams(3,m), metabolicParams(4,m), ...
                   metabolicParams(5,m)];
    Edot(:,m) = calcUmbergerCost2010(max(0, Lce), ...
                                     Vce, ...
                                     max(0, muscleData.Fce(:,m)), ...
                                     max(0, muscleData.FMltilde(:,m)), ...
                                     min(max(0, a(:,m)), 1), ... % TODO: keep using activation in place of excitation?
                                     min(max(0, a(:,m)), 1), ...
                                     paramStruct);
end

w_aT = 1000;
w_Edot = 1/(auxdata.model_mass*9.81*1.25);
% Overwriting first time point of metabolics to avoid effects that initial spikes
% in fiber powers may have on the cost.
Edot(1,:) = Edot(2,:);
phaseout.integrand = w_Edot*sum(Edot, 2) + w_aT.*sum(aT.^2,2) + sum((vA/100).^2,2);


