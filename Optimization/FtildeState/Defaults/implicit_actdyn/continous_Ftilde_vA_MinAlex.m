function phaseout = continous_Ftilde_vA_MinAlex(input)

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
dFtilde  = 10*input.phase.control(:,NMuscles+Ndof+1:end);

% Get states
a       = input.phase.state(:,1:NMuscles);
Ftilde = input.phase.state(:,NMuscles+1:end);

% PATH CONSTRAINTS
% Activation dynamics - De Groote et al. (2009)
act1 = vA + a./(ones(size(a,1),1)*tauDeact);
act2 = vA + a./(ones(size(a,1),1)*tauAct);

% Hill-equilibrium constraint
[Hilldiff,F,~,~,vMtilde] = DeGroote2016Muscle_FtildeState(a,Ftilde,dFtilde,splinestruct.LMT,splinestruct.VMT,params,input.auxdata.Fvparam,input.auxdata.Fpparam,input.auxdata.Faparam);

% Moments constraint
Topt = 150;
Tdiff = zeros(numColPoints,Ndof);
for dof = 1:Ndof
    T_exp=splinestruct.ID(:,dof);
    index_sel=(dof-1)*(NMuscles)+1:(dof-1)*(NMuscles)+NMuscles;
    T_sim=sum(F.*splinestruct.MA(:,index_sel),2) + Topt*aT(:,dof);
    Tdiff(:,dof) =  (T_exp-T_sim);
end

phaseout.path = [Tdiff Hilldiff act1 act2];

% DYNAMIC CONSTRAINTS
% Activation dynamics is implicit
% Contraction dynamics is implicit
phaseout.dynamics = [vA dFtilde];

% OBJECTIVE FUNCTION
% Calculate metabolic rate from Minetti & Alexander (1997) model
vmax = params(5,:);  
Fo = params(1,:);   
Edot = zeros(numColPoints,NMuscles);
for m = 1:NMuscles
    v = vmax(1,m)*vMtilde(:,m);
    Edot(:,m) = calcMinettiAlexanderProbe(v,vmax(1,m),Fo(1,m),a(:,m));
end

w1 = 10000;
w2 = 0.01;
phaseout.integrand = sum(Edot,2)+ w1.*sum(aT.^2,2)+ w2*sum((vA/100).^2,2);





