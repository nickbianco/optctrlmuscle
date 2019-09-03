function phaseout = continous_Ftilde(input)

% Get input data
NMuscles        = input.auxdata.NMuscles;
Ndof            = input.auxdata.Ndof;
tauAct          = input.auxdata.tauAct;
tauDeact        = input.auxdata.tauDeact;
params          = input.auxdata.params;
splinestruct    = input.auxdata.splinestruct;
numColPoints    = size(input.phase.state,1);
metabolicParams = input.auxdata.metabolicParams;

% Get controls
e   = input.phase.control(:,1:NMuscles);
aT  = input.phase.control(:,NMuscles+1:NMuscles+Ndof);
dFtilde  = 10*input.phase.control(:,NMuscles+Ndof+1:end);

% Get states
a       = input.phase.state(:,1:NMuscles);
Ftilde = input.phase.state(:,NMuscles+1:end);

% PATH CONSTRAINTS
% Hill-equilibrium constraint
% [Hilldiff,F,~,~,~] = DeGroote2016Muscle_FtildeState(a,Ftilde,dFtilde,...
%     splinestruct.LMT,splinestruct.VMT,params,input.auxdata.Fvparam,...
%     input.auxdata.Fpparam,input.auxdata.Faparam);
muscleData = DeGroote2016Muscle_FtildeState(a,Ftilde,dFtilde,...
    splinestruct.LMT,splinestruct.VMT,params,input.auxdata.Fvparam,...
    input.auxdata.Fpparam,input.auxdata.Faparam);


% Moments constraint
Topt = 150;
Tdiff = zeros(numColPoints,Ndof);
% issue is here!!

for dof = 1:Ndof
    T_exp=splinestruct.ID(:,dof);
    index_sel=(dof-1)*(NMuscles)+1:(dof-1)*(NMuscles)+NMuscles;
    T_sim=sum(muscleData.FT.*splinestruct.MA(:,index_sel),2) + Topt*aT(:,dof);
    Tdiff(:,dof) =  (T_exp-T_sim);
end

% phaseout.path = [Tdiff Hilldiff];
phaseout.path = [Tdiff muscleData.err];

% DYNAMIC CONSTRAINTS
% Activation dynamics

dadt = ones(numColPoints,NMuscles);
for m = 1:NMuscles
    dadt(:,m) = ActivationDynamics(e(:,m),a(:,m),tauAct(m),tauDeact(m),input.auxdata.b);
end

% Contraction dynamics is implicit
phaseout.dynamics = [dadt dFtilde];

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
w_Reg = 1e-3;
phaseout.integrand = w_Edot*sum(Edot, 2) + w_Res*sum(aT.^2,2) + ...
                     w_Reg*(sum(dFtilde.^2,2) + sum(e.^2,2));
     

% Calculate metabolic rate from Koelewijn et al. 2018 smooth approximation of
% Umberger model
% gotimag = false;
% Edot = zeros(numColPoints, NMuscles);
% for m = 1:NMuscles
%     
%     Lce = muscleData.lMtilde(:,m)*params(2,m);
%     Vce = muscleData.vMtilde(:,m)*params(5,m);
%     mass = metabolicParams(4,m);
%     paramStruct = [metabolicParams(1,m), metabolicParams(2,m), ...
%                    metabolicParams(3,m), mass, ...
%                    metabolicParams(5,m)];
%                
%     Edot(:,m) = mass * calcUmbergerKoelewijn2018Cost(Lce, Vce, ...
%             muscleData.Fce(:,m), muscleData.FMltilde(:,m), e(:,m), a(:,m), ...
%             paramStruct);
%     for i=1:length(Edot(:,m))
%         if imag(Edot(i,m))
%             gotimag = true;
%         end
%     end
%     if gotimag
%         fprintf('\nWe found it!!!\n')
%         keyboard
%     end
%     
% end

% w_aT = 1000; % reserve actuator term weight
% w_a = 0.05;
% w_vA = 0.05;
% % since it is normalized by muscle mass coming out of the function try not
% % weighting it!
% w_Edot = .1; % 1/(input.auxdata.model_mass*9.81*1.25); % metabolic cost term weight #TODO why different from nick's???
% 
% w_reg = 0.001; % regularization term weight
% 
% phaseout.integrand = w_Edot*sum(Edot, 2) + w_aT.*sum(aT.^2,2)+ w_a*sum(a.^2,2) + ...
%     w_vA*sum((vA/100).^2,2) + w_reg*(sum(dFtilde.^2, 2)); % + sum(e.^2, 2));
 

