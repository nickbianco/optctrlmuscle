function phaseout = continous_lMtildeSynHipAnkle_Exc_Act(input)

% Get input data
NMuscles        = input.auxdata.NMuscles;
Ndof            = input.auxdata.Ndof;
tauAct          = input.auxdata.tauAct;
tauDeact        = input.auxdata.tauDeact;
params          = input.auxdata.params;
splinestruct    = input.auxdata.splinestruct;
numColPoints    = size(input.phase.state,1);
Fopt_exo        = input.auxdata.Fopt_exo;
numControls     = input.auxdata.numControls;

% Get controls
u       = input.phase.control(:,1:numControls);
aT      = input.phase.control(:,numControls+1:numControls+Ndof);
% vMtilde = input.phase.control(:,numControls+Ndof+1:end-1);
aD      = input.phase.control(:,end);


% Get states
a       = input.phase.state;
% lMtilde = input.phase.state(:,NMuscles+1:end);

% Get tradeoff parameter info
alpha    = input.phase.parameter(:,1);
tradeoff = input.auxdata.tradeoff;

% Synergy control structure
synActivations = u;
synVecs = input.phase.parameter(:,2:end);
synVectors = reshape(synVecs(1,:),numControls,NMuscles);

e = synActivations*synVectors;

synUnitMag = zeros(numColPoints,numControls);
for syn = 1:numControls
    idx = NMuscles*(syn-1) + 1;
    synVec = synVecs(:,idx:(idx+NMuscles-1));
    synUnitMag(:,syn) = sum(synVec,2);
end

% synVectors = input.phase.parameter(:,2:end);
% synUnitMag = zeros(numColPoints,numControls);
% e = zeros(numColPoints,NMuscles);
% for syn = 1:numControls
%     idx = NMuscles*(syn-1) + 1;
%     synVec = synVectors(:,idx:(idx+NMuscles-1));
%     synUnitMag(:,syn) = sum(synVec,2);
%     synAct = u(:,syn);
% 
%     for m = 1:NMuscles
%         e(:,m) = e(:,m) + synAct(:,1).*synVec(:,m);
%     end
% end

% PATH CONSTRAINTS
% Hill-equilibrium constraint
[F, ~, ~, ~, ~, ~] = HillModel_RigidTendon(e, splinestruct.LMT, ... 
                                              splinestruct.VMT, ... 
                                              params, ... 
                                              input.auxdata.Fvparam, ...
                                              input.auxdata.Fpparam, ...
                                              input.auxdata.Faparam);
                                    
% [Hilldiff, F] = DeGroote2016Muscle_lMtildeState(a,lMtilde,vMtilde,splinestruct.LMT,params,input.auxdata.Fvparam,input.auxdata.Fpparam,input.auxdata.Faparam);

% Moments constraint
Topt = 150;
r = 0.1;
Tdiff = zeros(numColPoints,Ndof);
for dof = 1:Ndof
    T_exp=splinestruct.ID(:,dof);
    T_exo= r*Fopt_exo(dof)*aD.*(ones(numColPoints,1)+tradeoff(dof)*alpha);
    index_sel=(dof-1)*(NMuscles)+1:(dof-1)*(NMuscles)+NMuscles;
    T_sim=sum(F.*splinestruct.MA(:,index_sel),2) + Topt*aT(:,dof) + T_exo;
    Tdiff(:,dof) =  (T_exp-T_sim);
end

% phaseout.path = [Tdiff Hilldiff synUnitMag];
phaseout.path = [Tdiff synUnitMag];

% DYNAMIC CONSTRAINTS
% Activation dynamics
% dadt = ones(numColPoints,NMuscles);
% for m = 1:NMuscles
%     dadt(:,m) = ActivationDynamics(e(:,m),a(:,m),tauAct(m),tauDeact(m),input.auxdata.b);
% end

% Solve activation dynamics for one muscle so GPOPS is happy
dadt = ActivationDynamics(e(:,1),a,tauAct(1),tauDeact(1),input.auxdata.b);

% Contraction dynamics is implicit
% dlMtildedt = 10*vMtilde;
% 
% phaseout.dynamics = [dadt dlMtildedt];
phaseout.dynamics = dadt;

% OBJECTIVE FUNCTION
w1 = 1000;
% phaseout.integrand = sum(e.^2,2) + sum(a.^2,2) + w1.*sum(aT.^2,2);
phaseout.integrand = sum(e.^2,2) + w1.*sum(aT.^2,2);