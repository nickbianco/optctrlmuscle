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
aD = input.phase.control(:,end-1);
s = input.phase.control(:,end);

% Get states
a       = input.phase.state(:,1:NMuscles);
Ftilde = input.phase.state(:,NMuscles+1:NMuscles+NMuscles);

% Get moment arms
exoMomentArms = zeros(numColPoints,6);
if input.auxdata.active.hip
    exoMomentArms(:,1) = input.phase.parameter(:,input.auxdata.active.hip);
end
if input.auxdata.active.knee
    exoMomentArms(:,2) = input.phase.parameter(:,input.auxdata.active.knee);
end
if input.auxdata.active.ankle
    exoMomentArms(:,3) = input.phase.parameter(:,input.auxdata.active.ankle);
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
Lexo = zeros(numColPoints,1);
for dof = 1:Ndof
    if input.auxdata.passive.hip && (dof==input.auxdata.hip_DOF)
        Lexo = Lexo + -exoMomentArms(:,4).*splinestruct.IK(:,input.auxdata.hip_DOF);
    end
    if input.auxdata.passive.knee && (dof==input.auxdata.knee_DOF)
        Lexo = Lexo + -exoMomentArms(:,5).*splinestruct.IK(:,input.auxdata.knee_DOF);
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
[Hilldiff,F,~,~] = DeGroote2016Muscle_FtildeState(a,Ftilde,dFtilde,splinestruct.LMT,splinestruct.VMT,params,input.auxdata.Fvparam,input.auxdata.Fpparam,input.auxdata.Faparam);

% Exosuit torques
% Calculate max active force based on subject mass
Fmax_act = 15*input.auxdata.model_mass; %  N/kg * kg
Texo_act_hip = Fmax_act*aD.*exoMomentArms(:,1);
Texo_act_knee = Fmax_act*aD.*exoMomentArms(:,2);
Texo_act_ankle = Fmax_act*aD.*exoMomentArms(:,3);

% Calculate passive force based on normalized exo path length
k = 100; % kN/m - van den Bogert 2013
Fexo_pass = k*(Lexo - exoSlackLength) + s;
Texo_pass_hip = Fexo_pass.*exoMomentArms(:,4);
Texo_pass_knee = Fexo_pass.*exoMomentArms(:,5);
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

phaseout.path = [Tdiff Hilldiff act1 act2 Fexo_pass];

% DYNAMIC CONSTRAINTS
% Activation dynamics is implicit
% Contraction dynamics is implicit
phaseout.dynamics = [vA dFtilde];

% OBJECTIVE FUNCTION
w1 = 1000;
w2 = 0.01;
phaseout.integrand = sum(a.^2,2)+ w1.*sum(aT.^2,2)+ w2*sum((vA/100).^2,2);


