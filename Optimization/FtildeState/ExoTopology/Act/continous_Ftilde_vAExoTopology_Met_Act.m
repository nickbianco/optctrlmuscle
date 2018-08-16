function phaseout = continous_Ftilde_vAExoTopology_Met_Act(input)

% Get input data
NMuscles        = input.auxdata.NMuscles;
Ndof            = input.auxdata.Ndof;
tauAct          = input.auxdata.tauAct;
tauDeact        = input.auxdata.tauDeact;
params          = input.auxdata.params;
metabolicParams = input.auxdata.metabolicParams;
splinestruct    = input.auxdata.splinestruct;
numColPoints    = size(input.phase.state,1);
T_exp           = splinestruct.ID;

% Get controls
vA      = 100*input.phase.control(:,1:NMuscles);
aT      = input.phase.control(:,NMuscles+1:NMuscles+Ndof);
dFtilde = 10*input.phase.control(:,NMuscles+Ndof+1:NMuscles+Ndof+NMuscles);
aD      = input.phase.control(:,end-(input.auxdata.numActiveDOFs-1):end);

% if false
%     time = input.phase.time;
%     starttime = time(1) + (time(end)-time(1))*0.15;
%     endtime = time(1) + (time(end)-time(1))*0.85;
%     up = (1 ./ (1 + exp(50 * ((starttime - time)))));
%     down = (1 ./ (1 + exp(50 * ((time - endtime)))));
%     aD = aD.*up.*down;
% end

% Get states
a      = input.phase.state(:,1:NMuscles);
Ftilde = input.phase.state(:,NMuscles+1:NMuscles+NMuscles);
if input.auxdata.shift_exo_peaks
    aDshift = input.phase.state(:,end);
end

% Convert parameters to the correct range
paramsLower = input.auxdata.paramsLower;
paramsUpper = input.auxdata.paramsUpper;
parameter = 0.5*(paramsUpper-paramsLower).*(input.phase.parameter+1) + paramsLower;

% Get moment arms and DOF controls
exoMomentArms = zeros(numColPoints,3);
aD_hip = zeros(numColPoints,1);
aD_knee = zeros(numColPoints,1);
aD_ankle = zeros(numColPoints,1);
% shift_ref = [];
if input.auxdata.active.hip
    exoMomentArms(:,1) = parameter(:,input.auxdata.active.hip);
    if input.auxdata.numActiveDOFs > 1
        aD_hip = aD(:,input.auxdata.active.hip);
    else
        aD_hip = aD;
    end
%     if input.auxdata.shift_exo_peaks
%         shift_ref = 'hip';
%     end 
end
if input.auxdata.active.knee
    exoMomentArms(:,2) = parameter(:,input.auxdata.active.knee);
    if input.auxdata.numActiveDOFs > 1
        aD_knee = aD(:,input.auxdata.active.knee);
    else
        aD_knee = aD;
    end
%     if input.auxdata.shift_exo_peaks
%         if isempty(shift_ref)
%             shift_ref = 'knee';
%         else
%             T_ref = sign(exoMomentArms(1,1))*T_exp(:,input.auxdata.hip_DOF);
%             T_knee = sign(exoMomentArms(1,2))*T_exp(:,input.auxdata.knee_DOF);
%             [~,refIdx] = max(T_ref);
%             [~,kneeIdx] = max(T_knee);
%             shift_factor = kneeIdx - refIdx;
%             aD_knee = ShiftCurve(aD_knee, shift_factor);
%         end
%     end
end
if input.auxdata.active.ankle
    exoMomentArms(:,3) = parameter(:,input.auxdata.active.ankle);
    if input.auxdata.numActiveDOFs > 1
        aD_ankle = aD(:,input.auxdata.active.ankle);
    else
        aD_ankle = aD;
    end
%     if input.auxdata.shift_exo_peaks
%         switch shift_ref
%             case 'knee'
%                 T_ref = sign(exoMomentArms(1,1))*T_exp(:,input.auxdata.hip_DOF);
%             case 'hip'
%                 T_ref = sign(exoMomentArms(1,2))*T_exp(:,input.auxdata.knee_DOF);
%         end
%         T_ankle = sign(exoMomentArms(1,3))*T_exp(:,input.auxdata.ankle_DOF);
%         [~,refIdx] = max(T_ref);
%         [~,ankleIdx] = max(T_ankle);
%         shift_factor = ankleIdx - refIdx;
%         aD_ankle = ShiftCurve(aD_ankle, shift_factor);
%     end
    if input.auxdata.shift_exo_peaks
       aD_ankle = aDshift;
    end
end

% PATH CONSTRAINTS
% Activation dynamics - De Groote et al. (2009)
act1 = vA + a./(ones(size(a,1),1)*tauDeact);
act2 = vA + a./(ones(size(a,1),1)*tauAct);

% Hill-equilibrium constraint
muscleData = DeGroote2016Muscle_FtildeState(a,Ftilde,dFtilde,splinestruct.LMT,splinestruct.VMT,params,input.auxdata.Fvparam,input.auxdata.Fpparam,input.auxdata.Faparam);

% Exosuit torques
Texo_act_hip = input.auxdata.Tmax_act.*aD_hip.*exoMomentArms(:,1);
Texo_act_knee = input.auxdata.Tmax_act.*aD_knee.*exoMomentArms(:,2).*input.auxdata.kneeAngleSign;
Texo_act_ankle = input.auxdata.Tmax_act.*aD_ankle.*exoMomentArms(:,3);

% Moments constraint
Topt = 150;
Tdiff = zeros(numColPoints,Ndof);
for dof = 1:Ndof
    index_sel=(dof-1)*(NMuscles)+1:(dof-1)*(NMuscles)+NMuscles;
    T_sim=sum(muscleData.FT.*splinestruct.MA(:,index_sel),2) + Topt*aT(:,dof);

    if dof==input.auxdata.hip_DOF
        T_sim = T_sim + Texo_act_hip;
    end
    if dof==input.auxdata.knee_DOF
        T_sim = T_sim + Texo_act_knee;
    end
    if dof==input.auxdata.ankle_DOF
        T_sim = T_sim + Texo_act_ankle;
    end
   
    Tdiff(:,dof) = (T_exp(:,dof)-T_sim);
end

phaseout.path = [Tdiff muscleData.err act1 act2];

% DYNAMIC CONSTRAINTS
% Muscle activation dynamics is implicit
% Contraction dynamics is implicit
phaseout.dynamics = [vA dFtilde];

if input.auxdata.shift_exo_peaks
    % Device activation dynamics is explicit
    tauActExo = parameter(1,end);
    dadt = ActivationDynamics(aD,aDshift,tauActExo,0.05,input.auxdata.b);
    phaseout.dynamics = [phaseout.dynamics dadt];
end

% OBJECTIVE FUNCTION
Edot = zeros(numColPoints, NMuscles);
for m = 1:NMuscles
    Lce = muscleData.lMtilde(:,m)*params(2,m);
    Vce = muscleData.vMtilde(:,m)*params(5,m);
%     u = computeExcitationRaasch(a(:,m), vA(:,m), tauDeact(m), tauAct(m));
    paramStruct = [metabolicParams(1,m), metabolicParams(2,m), ...
                   metabolicParams(3,m), metabolicParams(4,m), ...
                   metabolicParams(5,m)];
    Edot(:,m) = calcUmbergerCost2010(max(0, Lce), ...
                                     Vce, ...
                                     max(0, muscleData.Fce(:,m)), ...
                                     max(0, muscleData.FMltilde(:,m)), ...
                                     min(max(0, a(:,m)), 1), ... % TODO: keep using activation in place of excitation?
                                     min(max(0, a(:,m)), 1), ...
                                     paramStruct);
end
w_aT = 1000;
w_a = 0.05;
w_vA = 0.05;
w_Edot = 1/(input.auxdata.model_mass*9.81*1.25);
% Overwriting first time point of metabolics to avoid effects that initial spikes
% in fiber powers may have on the cost.
Edot(1,:) = Edot(2,:);
phaseout.integrand = w_Edot*sum(Edot, 2) + w_aT.*sum(aT.^2,2) + sum((vA/100).^2,2);


