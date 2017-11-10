function phaseout = continous_SO_ExoTopology_Pass(input)

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

% Get states
a = input.phase.state;

% Get moment arms
exoMomentArms = zeros(numColPoints,3);
if input.auxdata.passive.hip
    exoMomentArms(:,1) = input.phase.parameter(:,input.auxdata.passive.hip);
end
if input.auxdata.passive.knee
    exoMomentArms(:,2) = input.phase.parameter(:,input.auxdata.passive.knee);
end
if input.auxdata.passive.ankle
    exoMomentArms(:,3) = input.phase.parameter(:,input.auxdata.passive.ankle);
end

% Slack length of passive elastic device
exoSlackLength = input.phase.parameter(:,end);

% Exosuit path length
Lexo = ones(numColPoints,1);
splinestruct.IK(:,input.auxdata.knee_DOF) = input.auxdata.kneeAngleSign*splinestruct.IK(:,input.auxdata.knee_DOF);
for dof = 1:Ndof
    if input.auxdata.passive.hip && (dof==input.auxdata.hip_DOF)
        Lexo = Lexo + -exoMomentArms(:,1).*splinestruct.IK(:,input.auxdata.hip_DOF);
    end
    if input.auxdata.passive.knee && (dof==input.auxdata.knee_DOF)
        Lexo = Lexo + -exoMomentArms(:,2).*splinestruct.IK(:,input.auxdata.knee_DOF);
    end
    if input.auxdata.passive.ankle && (dof==input.auxdata.ankle_DOF)
        Lexo = Lexo + -exoMomentArms(:,3).*splinestruct.IK(:,input.auxdata.ankle_DOF);
    end
end

% PATH CONSTRAINTS

% Get muscle forces
[F, ~, ~, ~, ~, ~] = HillModel_RigidTendon(e, splinestruct.LMT, ... 
                                              splinestruct.VMT, ... 
                                              params, ... 
                                              input.auxdata.Fvparam, ...
                                              input.auxdata.Fpparam, ...
                                              input.auxdata.Faparam);

% Exosuit torques
% Calculate passive force based on normalized exo path length
k = input.auxdata.passiveStiffness;
positiveStiffnessAboveLslack = (1 ./ (1 + exp(100 * ((exoSlackLength+0.075) - Lexo))));
Fexo_pass = k*(Lexo - exoSlackLength) .* positiveStiffnessAboveLslack;
Texo_pass_hip = Fexo_pass.*exoMomentArms(:,1);
Texo_pass_knee = Fexo_pass.*exoMomentArms(:,2);
Texo_pass_ankle = Fexo_pass.*exoMomentArms(:,3);

% Moments constraint
Topt = 150;
Tdiff = zeros(numColPoints,Ndof);
for dof = 1:Ndof
    T_exp=splinestruct.ID(:,dof);
    index_sel=(dof-1)*(NMuscles)+1:(dof-1)*(NMuscles)+NMuscles;
    T_sim=sum(F.*splinestruct.MA(:,index_sel),2) + Topt*aT(:,dof);

    if dof==input.auxdata.hip_DOF
        T_sim = T_sim + Texo_pass_hip;
    end
    if dof==input.auxdata.knee_DOF
        T_sim = T_sim + Texo_pass_knee;
    end
    if dof==input.auxdata.ankle_DOF
        T_sim = T_sim + Texo_pass_ankle;
    end
   
    Tdiff(:,dof) = (T_exp-T_sim);
end

phaseout.path = Tdiff;

% DYNAMIC CONSTRAINTS
% Solve activation dynamics for one muscle so GPOPS is happy
dadt = ActivationDynamics(e(:,1),a,tauAct(1),tauDeact(1),input.auxdata.b);

phaseout.dynamics = dadt;

% OBJECTIVE FUNCTION
w1 = 1000;
phaseout.integrand = sum(a.^2,2)+ w1*sum(aT.^2,2);


