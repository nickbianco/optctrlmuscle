function phaseout = continous_lMtildeISBHOBL_Exc_Act(input)

% Get input data
NMuscles        = input.auxdata.NMuscles;
Ndof            = input.auxdata.Ndof;
tauAct          = input.auxdata.tauAct;
tauDeact        = input.auxdata.tauDeact;
params          = input.auxdata.params;
splinestruct    = input.auxdata.splinestruct;
numColPoints    = size(input.phase.state,1);

% Get controls
e       = input.phase.control(:,1:NMuscles);
aT      = input.phase.control(:,NMuscles+1:NMuscles+Ndof);
vMtilde = input.phase.control(:,NMuscles+Ndof+1:end);

% Get states
a       = input.phase.state(:,1:NMuscles);
lMtilde = input.phase.state(:,NMuscles+1:end);

% Get parameters
normSpringStiff = input.phase.parameter(:, 1);
deviceRestLength = input.phase.parameter(:, 2);

% PATH CONSTRAINTS
% Hill-equilibrium constraint
[Hilldiff, FT, ~, ~] = DeGroote2016Muscle_lMtildeState(a,lMtilde,vMtilde,splinestruct.LMT,params,input.auxdata.Fvparam,input.auxdata.Fpparam,input.auxdata.Faparam);

% Moments constraint
Topt = 150;
maxSpringStiff = 500; % N/m
Tdiff = zeros(numColPoints,Ndof);
for dof = 1:Ndof
    T_exp = splinestruct.ID(:,dof);
    T_exo_ma = splinestruct.MAEXO(:,dof);
    delta_x = splinestruct.LenEXO-deviceRestLength;
    delta_x = max(delta_x,0);
    T_exo = maxSpringStiff*normSpringStiff.*T_exo_ma.*delta_x;
    index_sel=(dof-1)*(NMuscles)+1:(dof-1)*(NMuscles)+NMuscles;
    T_sim=sum(FT.*splinestruct.MA(:,index_sel),2) + Topt*aT(:,dof) + T_exo;
    Tdiff(:,dof) =  (T_exp-T_sim);
end
phaseout.path = [Tdiff Hilldiff];

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








