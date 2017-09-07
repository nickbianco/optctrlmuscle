function phaseout = continous_Ftilde_vAExoTopology(input)

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

% Get parameters
if isfield(input.auxdata,'passive')
   exoMomentArms = input.phase.parameter(:,1:end-1);
   exoSlackLength = input.pahse.parameter(:,end); 
else
   exoMomentArms = input.phase.parameter; 
end

% PATH CONSTRAINTS
% Activation dynamics - De Groote et al. (2009)
act1 = vA + a./(ones(size(a,1),1)*tauDeact);
act2 = vA + a./(ones(size(a,1),1)*tauAct);

% Hill-equilibrium constraint
[Hilldiff,F,~,~] = DeGroote2016Muscle_FtildeState(a,Ftilde,dFtilde,splinestruct.LMT,splinestruct.VMT,params,input.auxdata.Fvparam,input.auxdata.Fpparam,input.auxdata.Faparam);

% Moments constraint
Topt = 150;
Tdiff = zeros(numColPoints,Ndof);
Texo = zeros(numColPoints,Ndof);
for dof = 1:Ndof
    T_exp=splinestruct.ID(:,dof);
    index_sel=(dof-1)*(NMuscles)+1:(dof-1)*(NMuscles)+NMuscles;
    T_sim=sum(F.*splinestruct.MA(:,index_sel),2) + Topt*aT(:,dof);
    
    Fact = 15*auxdata.model_mass; %  N/kg * kg
    if isfield(input.auxdata,'active')
        active = input.auxdata.active;
        if active.hip && (dof==auxdata.hip_DOF)
            
        end
        if active.knee && (dof==auxdata.knee_DOF)
            
        end
        if active.ankle && (dof==auxdata.ankle_DOF)
            
        end
    end
   
    Fpass = 15*auxdata.model_mass; % N/kg * kg   
    if isfield(input.auxdata,'passive')
        passive = input.auxdata.passive;
        if passive.hip && (dof==auxdata.hip_DOF)
            
        end
        if passive.knee && (dof==auxdata.knee_DOF)
            
        end
        if passive.ankle && (dof==auxdata.ankle_DOF)
            
        end
    end
    
    Tdiff(:,dof) =  (T_exp-T_sim);
end

phaseout.path = [Tdiff Hilldiff act1 act2];

% DYNAMIC CONSTRAINTS
% Activation dynamics is implicit
% Contraction dynamics is implicit
phaseout.dynamics = [vA dFtilde];

% OBJECTIVE FUNCTION
w1 = 1000;
w2 = 0.01;
phaseout.integrand = sum(a.^2,2)+ w1.*sum(aT.^2,2)+ w2*sum((vA/100).^2,2);


