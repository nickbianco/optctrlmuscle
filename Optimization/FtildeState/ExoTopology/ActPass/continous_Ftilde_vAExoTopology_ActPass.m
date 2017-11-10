function phaseout = continous_Ftilde_vAExoTopology_ActPass(input)

% Get input data
NMuscles        = input.auxdata.NMuscles;
Ndof            = input.auxdata.Ndof;
tauAct          = input.auxdata.tauAct;
tauDeact        = input.auxdata.tauDeact;
params          = input.auxdata.params;
splinestruct    = input.auxdata.splinestruct;
numColPoints    = size(input.phase.state,1);

% Get controls
vA   = 100*input.phase.control(:,1:NMuscles);
aT  = input.phase.control(:,NMuscles+1:NMuscles+Ndof);
dFtilde  = 10*input.phase.control(:,NMuscles+Ndof+1:NMuscles+Ndof+NMuscles);
aD = input.phase.control(:,end-(input.auxdata.numActiveDOFs-1):end);

% Get states
a       = input.phase.state(:,1:NMuscles);
Ftilde = input.phase.state(:,NMuscles+1:NMuscles+NMuscles);

% Get moment arms
exoMomentArms = zeros(numColPoints,6);
aD_hip = zeros(numColPoints,1);
aD_knee = zeros(numColPoints,1);
aD_ankle = zeros(numColPoints,1);
if input.auxdata.active.hip
    exoMomentArms(:,1) = input.phase.parameter(:,input.auxdata.active.hip);
    if input.auxdata.numActiveDOFs > 1
        aD_hip = aD(:,input.auxdata.active.hip);
    else
        aD_hip = aD;
    end
end
if input.auxdata.active.knee
    exoMomentArms(:,2) = input.phase.parameter(:,input.auxdata.active.knee);
    if input.auxdata.numActiveDOFs > 1
        aD_knee = aD(:,input.auxdata.active.knee);
    else
        aD_knee = aD;
    end
end
if input.auxdata.active.ankle
    exoMomentArms(:,3) = input.phase.parameter(:,input.auxdata.active.ankle);
    if input.auxdata.numActiveDOFs > 1
        aD_ankle = aD(:,input.auxdata.active.ankle);
    else
        aD_ankle = aD;
    end
end
if input.auxdata.passive.hip
    exoMomentArms(:,4) = input.phase.parameter(:,input.auxdata.passive.hip);
end
if input.auxdata.passive.knee
    exoMomentArms(:,5) = input.phase.parameter(:,input.auxdata.passive.knee);
end
if input.auxdata.passive.ankle
    exoMomentArms(:,6) = input.phase.parameter(:,input.auxdata.passive.ankle);
end

% Slack length of passive elastic device
exoSlackLength = input.phase.parameter(:,end);

% Exosuit path length
Lexo = ones(numColPoints,1);
IK_knee = input.auxdata.kneeAngleSign*splinestruct.IK(:,input.auxdata.knee_DOF);
for dof = 1:Ndof
    if input.auxdata.passive.hip && (dof==input.auxdata.hip_DOF)
        Lexo = Lexo + -exoMomentArms(:,4).*splinestruct.IK(:,input.auxdata.hip_DOF);
    end
    if input.auxdata.passive.knee && (dof==input.auxdata.knee_DOF)
        Lexo = Lexo + -exoMomentArms(:,5).*IK_knee;
    end
    if input.auxdata.passive.ankle && (dof==input.auxdata.ankle_DOF)
        Lexo = Lexo + -exoMomentArms(:,6).*splinestruct.IK(:,input.auxdata.ankle_DOF);
    end
end

% PATH CONSTRAINTS
% Activation dynamics - De Groote et al. (2009)
act1 = vA + a./(ones(size(a,1),1)*tauDeact);
act2 = vA + a./(ones(size(a,1),1)*tauAct);

% Hill-equilibrium constraint
[Hilldiff,F,~,~,~] = DeGroote2016Muscle_FtildeState(a,Ftilde,dFtilde,splinestruct.LMT,splinestruct.VMT,params,input.auxdata.Fvparam,input.auxdata.Fpparam,input.auxdata.Faparam);

% Exosuit torques
% Active device
Texo_act_hip = input.auxdata.Fmax_act*aD_hip.*exoMomentArms(:,1);
Texo_act_knee = input.auxdata.Fmax_act*aD_knee.*exoMomentArms(:,2)*input.auxdata.kneeAngleSign;
Texo_act_ankle = input.auxdata.Fmax_act*aD_ankle.*exoMomentArms(:,3);

% Calculate passive force based on normalized exo path length
k = input.auxdata.passiveStiffness;
positiveStiffnessAboveLslack = (1 ./ (1 + exp(100 * ((exoSlackLength+0.075) - Lexo))));
Fexo_pass = k*(Lexo - exoSlackLength) .* positiveStiffnessAboveLslack;
Texo_pass_hip = Fexo_pass.*exoMomentArms(:,4);
Texo_pass_knee = Fexo_pass.*exoMomentArms(:,5)*input.auxdata.kneeAngleSign;
Texo_pass_ankle = Fexo_pass.*exoMomentArms(:,6);

% Moments constraint
Topt = 150;
Tdiff = zeros(numColPoints,Ndof);
for dof = 1:Ndof
    T_exp=splinestruct.ID(:,dof);
    index_sel=(dof-1)*(NMuscles)+1:(dof-1)*(NMuscles)+NMuscles;
    T_sim=sum(F.*splinestruct.MA(:,index_sel),2) + Topt*aT(:,dof);

    if dof==input.auxdata.hip_DOF
        T_sim = T_sim + Texo_act_hip + Texo_pass_hip;
    end
    if dof==input.auxdata.knee_DOF
        T_sim = T_sim + Texo_act_knee + Texo_pass_knee;
    end
    if dof==input.auxdata.ankle_DOF
        T_sim = T_sim + Texo_act_ankle + Texo_pass_ankle;
    end
   
    Tdiff(:,dof) = (T_exp-T_sim);
end

phaseout.path = [Tdiff Hilldiff act1 act2];

% DYNAMIC CONSTRAINTS
% Activation dynamics is implicit
% Contraction dynamics is implicit
phaseout.dynamics = [vA dFtilde];

% OBJECTIVE FUNCTION
w1 = 1000;
w2 = 0.01;
phaseout.integrand = sum(a.^2,2) + w1*sum(aT.^2,2) + w2*sum((vA/100).^2,2);


