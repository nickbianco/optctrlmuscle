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
aT      = input.phase.control(:,NMuscles+1:NMuscles+Ndof);
vMtilde = input.phase.control(:,NMuscles+Ndof+1:end-1);
aD      = input.phase.control(:,end);

% Get states
a       = input.phase.state(:,1:NMuscles);
lMtilde = input.phase.state(:,NMuscles+1:end);

% Get tradeoff parameter info
alpha    = input.phase.parameter(:,1);
tradeoff = input.auxdata.tradeoff;

% Synergy control structure
synVectors = input.phase.parameter(:,2:end);
synUnitMag = zeros(numColPoints,numControls);
e = zeros(numColPoints,NMuscles);
for syn = 1:numControls
    idx = NMuscles*(syn-1) + 1;
    synVec = synVectors(:,idx:(idx+NMuscles-1));
    synUnitMag(:,syn) = sum(synVec,2);
    synAct = u(:,syn);
    
    e = e + times(synAct,synVec);   
end

% PATH CONSTRAINTS
% Hill-equilibrium constraint
[Hilldiff, F] = DeGroote2016Muscle_lMtildeState(a,lMtilde,vMtilde,splinestruct.LMT,params,input.auxdata.Fvparam,input.auxdata.Fpparam,input.auxdata.Faparam);

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

phaseout.path = [Tdiff Hilldiff synUnitMag];

% DYNAMIC CONSTRAINTS
% Activation dynamics
dadt = ones(numColPoints,NMuscles);
for m = 1:NMuscles
    dadt(:,m) = ActivationDynamics(e(:,m),a(:,m),tauAct(m),tauDeact(m),input.auxdata.b);
end

% Contraction dynamics is implicit
dlMtildedt = 10*vMtilde;

phaseout.dynamics = [dadt dlMtildedt];

% OBJECTIVE FUNCTION
w1 = 1000;
phaseout.integrand = sum(e.^2,2) + sum(a.^2,2) + w1.*sum(aT.^2,2);