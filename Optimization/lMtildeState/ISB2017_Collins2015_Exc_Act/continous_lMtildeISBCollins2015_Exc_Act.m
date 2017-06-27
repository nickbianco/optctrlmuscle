function phaseout = continous_lMtildeISBCollins2015_Exc_Act(input)

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
springRestAngle = input.phase.parameter(:, 2);

% PATH CONSTRAINTS
% Hill-equilibrium constraint
[Hilldiff, FT, ~, ~] = DeGroote2016Muscle_lMtildeState(a,lMtilde,vMtilde,splinestruct.LMT,params,input.auxdata.Fvparam,input.auxdata.Fpparam,input.auxdata.Faparam);

% Use composition of logistic functions to turn off the spring when not 
% active
%
% The first logistic function transitions from "off-to-on" when the 
% first peak in the ankle ankle is reached, which is when the spring begins 
% to stretch. The second logistic function transitions from "on-to-off"
% when the rest length in the spring is reached after the spring recoils to
% deliver assistance before pushoff, which is when the string attached to 
% the spring begins to take on slack, ensuring a uni-directional 
% assistance.
first_peak = input.auxdata.rest_length_first_peak;
after_recoil = input.auxdata.rest_length_after_recoil;
beginSpringStretching = 1 ./ (1 + exp(100 * (first_peak  - input.phase.time)));
restLengthReached = 1 ./ (1 + exp(100 * (input.phase.time - after_recoil)));

isSpringActive = beginSpringStretching .* restLengthReached;

% Moments constraint
Topt = 150;
maxSpringStiff = 400; % N-m/rad.
Tdiff = zeros(numColPoints,Ndof);
for dof = 1:Ndof
    T_exp=splinestruct.ID(:,dof);
    index_sel=(dof-1)*(NMuscles)+1:(dof-1)*(NMuscles)+NMuscles;
    T_sim=sum(FT.*splinestruct.MA(:,index_sel),2) + Topt*aT(:,dof);
    if any(dof == input.auxdata.clutched_spring_dofs)
        springStretch = -(splinestruct.IK(:,dof) - springRestAngle);
        T_sim = T_sim + maxSpringStiff * normSpringStiff .* springStretch .* isSpringActive;
    end
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








