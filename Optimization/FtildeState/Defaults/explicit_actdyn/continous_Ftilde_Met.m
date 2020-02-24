function phaseout = continous_Ftilde_Met(input)

% Get input data
auxdata         = input.auxdata;
NMuscles        = auxdata.NMuscles;
Ndof            = auxdata.Ndof;
tauAct          = auxdata.tauAct;
tauDeact        = auxdata.tauDeact;
params          = auxdata.params;
metabolicParams = auxdata.metabolicParams;
splinestruct    = auxdata.splinestruct;
speed           = auxdata.speed;
numColPoints    = size(input.phase.state,1);
NStates         = size(input.phase.state,2);
NControls       = size(input.phase.state,2);

% Get controls
e   = input.phase.control(:,1:NMuscles);
aT  = input.phase.control(:,NMuscles+1:NMuscles+Ndof);
dFtilde  = 10*input.phase.control(:,NMuscles+Ndof+1:end);

% Get states
a       = input.phase.state(:,1:NMuscles);
Ftilde = input.phase.state(:,NMuscles+1:end);

% PATH CONSTRAINTS
% Hill-equilibrium constraint
muscleData = DeGroote2016Muscle_FtildeState(a,Ftilde,dFtilde,splinestruct.LMT,... 
                splinestruct.VMT,params,auxdata.Fvparam,...
                auxdata.Fpparam,auxdata.Faparam);

% Moments constraint
Topt = 150;
Tdiff = zeros(numColPoints,Ndof);
for dof = 1:Ndof
    T_exp=splinestruct.ID(:,dof);
    index_sel=(dof-1)*(NMuscles)+1:(dof-1)*(NMuscles)+NMuscles;
    T_sim=sum(muscleData.FT.*splinestruct.MA(:,index_sel),2) + Topt*aT(:,dof);
    Tdiff(:,dof) =  (T_exp-T_sim);
end

phaseout.path = [Tdiff muscleData.err];

% DYNAMIC CONSTRAINTS
% Activation dynamics
dadt = ones(numColPoints,NMuscles);
for m = 1:NMuscles
    dadt(:,m) = ActivationDynamics(e(:,m),a(:,m),tauAct(m),tauDeact(m),auxdata.b);
end

% Contraction dynamics is implicit
phaseout.dynamics = [dadt dFtilde];

% OBJECTIVE FUNCTION
% Calculate metabolic rate from Minetti & Alexander (1997) model
% vmax = params(5,:);  
% Fo = params(1,:);   
% Edot = zeros(numColPoints,NMuscles);
% for m = 1:NMuscles
%     v = vmax(1,m)*muscleData.vMtilde(:,m);
%     Edot(:,m) = calcMinettiAlexanderProbe(v,vmax(1,m),Fo(1,m),a(:,m));
% end

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

w_Edot = 1/(NMuscles*auxdata.model_mass);
w_Res = 1e3/Ndof;
w_Control = 1/(3*NMuscles);
% w_Reg = auxdata.regularizationWeight/(NMuscles + Ndof);

% dFtilde_diff = [zeros(1, NMuscles); dFtilde(2:end,:) - dFtilde(1:end-1,:)];
% e_diff = [zeros(1, NMuscles); e(2:end,:) - e(1:end-1,:)];
% Ftilde_diff = [zeros(1, NMuscles); Ftilde(2:end,:) - Ftilde(1:end-1,:)];
% a_diff = [zeros(1, NMuscles); a(2:end,:) - a(1:end-1,:)];
% aT_diff = [zeros(1, Ndof); aT(2:end,:) - aT(1:end-1,:)];

phaseout.integrand = w_Edot*sum(Edot, 2) + w_Res*sum(aT.^2,2) + ...
                     w_Control*(sum((dFtilde/10).^2,2) + sum(e.^2,2) + sum(a.^2,2));

% phaseout.integrand = w_Edot*sum(Edot, 2) + w_Res*sum(aT.^2,2) + ...
%                      w_Reg*(sum((dFtilde/10).^2,2) + sum(e.^2,2) + ...
%                             sum(a_diff.^2,2) + sum(Ftilde_diff.^2,2) + ...
%                             sum(aT_diff.^2, 2)); 