function phaseout = continous_Ftilde_vA_Met(input)

% Get input data
NMuscles        = input.auxdata.NMuscles;
Ndof            = input.auxdata.Ndof;
tauAct          = input.auxdata.tauAct;
tauDeact        = input.auxdata.tauDeact;
params          = input.auxdata.params;
metabolicParams = input.auxdata.metabolicParams;
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
muscleData = DeGroote2016Muscle_FtildeState(a,Ftilde,dFtilde,splinestruct.LMT,splinestruct.VMT,params,input.auxdata.Fvparam,input.auxdata.Fpparam,input.auxdata.Faparam);

% Moments constraint
Topt = 150;
Tdiff = zeros(numColPoints,Ndof);
for dof = 1:Ndof
    T_exp=splinestruct.ID(:,dof);
    T_exo=splinestruct.EXO(:,dof);
    index_sel=(dof-1)*(NMuscles)+1:(dof-1)*(NMuscles)+NMuscles;
    T_sim=sum(muscleData.FT.*splinestruct.MA(:,index_sel),2) + Topt*aT(:,dof) + T_exo;
    Tdiff(:,dof) =  (T_exp-T_sim);
end

phaseout.path = [Tdiff muscleData.err act1 act2];

% DYNAMIC CONSTRAINTS
% Activation dynamics is implicit
% Nick has some different stuff going on here!!!
% Contraction dynamics is implicit
phaseout.dynamics = [vA dFtilde];

% OBJECTIVE FUNCTION
% Calculate metabolic rate from Minetti & Alexander (1997) model
vmax = params(5,:);  
Fo = params(1,:);   
Edot = zeros(numColPoints,NMuscles);
for m = 1:NMuscles
    v = vmax(1,m)*muscleData.vMtilde(:,m);
    Edot(:,m) = calcMinettiAlexanderProbe(v,vmax(1,m),Fo(1,m),a(:,m));
end

w_Edot = 1/input.auxdata.model_mass;
w_Res = 5e3;
w_Reg = 1e-4;



phaseout.integrand = w_Edot*sum(Edot, 2) + w_Res*sum(aT.^2,2) + ...
                     w_Reg*(sum(dFtilde.^2,2) + sum(vA.^2,2)); % + sum(e.^2,2));


% calculate the metabolic rate from umberger model
% Edot = zeros(numColPoints, NMuscles);
% for m = 1:NMuscles
%     Lce = muscleData.lMtilde(:,m)*params(2,m);
%     Vce = muscleData.vMtilde(:,m)*params(5,m);
%     % u = computeExcitationRaasch(a(:,m), vA(:,m), tauDeact(m), tauAct(m));
%     paramStruct = [metabolicParams(1,m), metabolicParams(2,m), ...
%                    metabolicParams(3,m), metabolicParams(4,m), ...
%                    metabolicParams(5,m)];    
%     Edot(:,m) = calcUmbergerCost2010(max(0, Lce), ...
%                                      Vce, ...
%                                      max(0, muscleData.Fce(:,m)), ...
%                                      max(0, muscleData.FMltilde(:,m)), ...
%                                      min(max(0, a(:,m)), 1), ... % replaced excitation input w/ activation
%                                      min(max(0, a(:,m)), 1), ...
%                                      paramStruct);
% end

