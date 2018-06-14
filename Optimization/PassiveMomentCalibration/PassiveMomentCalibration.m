function [calibratedModifiers] = PassiveMomentCalibration(model_fpath, auxdata, DatStore, paramsToCal)

import org.opensim.modeling.*
model = Model(model_fpath);

% Get coordinate names -- we will need to set kinematics of all coordinates
% for muscle length and moment arm calculations.
coordSet = model.getCoordinateSet();
coordNames = ArrayStr();
coordSet.getNames(coordNames);
allCoords = cell(0);
for i = 1:coordNames.getSize()
   allCoords{i} = char(coordNames.get(i-1));
end

% Load digitized joint angle and moment curves
load DigitizedPassiveMoments_Silder2007.mat
q_Silder2007 = JAngles;
M_Silder2007 = PassiveM;
M = createMomentArray(M_Silder2007);
[q, u] = createJointAngleVelocityArrays(q_Silder2007, allCoords);
[B, lMT] = getLengthsMomentArms(model, q, u, DatStore.DOFNames, DatStore.MuscleNames);

% Get info about muscles to calibrate
muscles = fieldnames(paramsToCal);
lMo = [];
lTs = [];
e0 = [];
for m = 1:length(muscles)
   idx = find(contains(DatStore.MuscleNames, muscles{m}));
   params = paramsToCal.(muscles{m}).params;   
   for p = 1:length(params)
      switch params{p}
          case 'optimal_fiber_length'
              lMo = [lMo idx];
          case 'tendon_slack_length'
              lTs = [lTs idx];
          case 'muscle_strain'
              e0 = [e0 idx];
      end
   end 
end
indices = struct();
indices.lMo = lMo;
indices.lTs = lTs;
indices.e0 = e0;
numlMo = length(lMo);
numlTs = length(lTs);
nume0 = length(e0);
indices.numlMo = numlMo;
indices.numlTs = numlTs;
indices.nume0 = nume0;
auxdata.indices = indices;

% Set up calibration optimization problem
numParams = numlMo + numlTs + nume0;
x0 = ones(numParams, 1);
lb = 0.75*x0;
ub = 1.25*x0;

% Rectus femoris passive muscle properties tend to be too stiff, even
% after calibration -- these modified bounds attempt to prevent this.
RFidx = find(contains(muscles,'rect_fem_r'));
lb(RFidx) = 1.0;
lb(RFidx+length(muscles)) = 1.0;
lb(RFidx+(2*length(muscles))) = 1.0;

options = optimoptions('fmincon', ...
                       'Display','iter', ...
                       'Algorithm','sqp', ...
                       'MaxFunctionEvaluations', 100000, ...
                       'Hessian', 'lbfgs');

% Muscle passive fiber force is governed by the expression seen in
% getPassiveForce, which depends on the muscle strain parameter e0. 
% For optimization variable scaling purposes, this variable is
% determined by the following linear equation for a given muscle:
%
% e0 = m*x + b, where x is the associated optimization variable for e0
%
% Below are the ranges of e0 for the same range of x, [0.75 1.25], and the 
% necessary linear equation coefficient values.
% e0: [0.6 0.75]    --> m = 0.3, b = 0.375
% e0: [0.6 0.8]     --> m = 0.4, b = 0.3
% e0: [0.6 0.7]     --> m = 0.2, b = 0.45
% e0: [0.55 0.65]   --> m = 0.2, b = 0.4
% e0: [0.525 0.675] --> m = 0.3, b = 0.3
% e0: [0.5 0.7]     --> m = 0.4, b = 0.2
m = 0.3;
b = 0.375;
auxdata.e0LinCoefs.m = m;
auxdata.e0LinCoefs.b = b;
 
% Choose whether or not to use rigid tendon assumption
% (compliant tendon optimization doesn't work yet, so not really a choice...)
isRigidTendon = true;

% Optimize
[x, res] = fmincon(@(x) obj(x, M, B, lMT, auxdata, isRigidTendon), x0, ...
                   [],[],[],[], lb, ub, ...
                   @(x) nonlincon(x, lMT, B, M, auxdata, isRigidTendon), options);

% Get results
TendonForce = getForce(x, lMT, auxdata, isRigidTendon);
MuscleMoments = zeros(size(M));
for dof = 1:size(M,2)
    MuscleMoments(:,dof) = sum(TendonForce.*squeeze(B(:,dof,:)), 2);
end
M_nz = M(M~=0);
MuscleMoments_nz = MuscleMoments(M~=0);

% create muscle strain modifiers
cal_e0 = m*x((numlMo+numlTs + 1):(numParams)) + b*ones(nume0,1);
mod_e0 = ones(nume0,1) + ((cal_e0-0.6)/0.6);

% Return calibrated parameter modifiers
calibratedModifiers.indices = indices;
calibratedModifiers.lMo = x(1:numlMo);
calibratedModifiers.lTs = x((numlMo+1):(numlMo+numlTs)); 
calibratedModifiers.e0 = mod_e0;

end

function plotMuscle(m, x, lMT, auxdata, DatStore)

[~, auxdata] = getForce(x, lMT, auxdata, true);
[FT, Fpe, lMtilde, cos_alpha] = getForceRigidTendon(lMT, auxdata.params);

mname = DatStore.MuscleNames{m}
lMo = auxdata.params(9,m)
lTs = auxdata.params(10,m)
e0 = auxdata.params(7,m)

% tendon force
figure;
plot(FT(:,m))
title(sprintf('%s: tendon force', mname))

% passive force
figure;
plot(Fpe(:,m))
title(sprintf('%s: passive force', mname))

% norm fiber length
figure;
plot(lMtilde(:,m))
title(sprintf('%s: norm fiber length', mname))

end

function [f] = obj(x, M, B, lMT, auxdata, isRigidTendon)

% Get force along tendon
[FT,~] = getForce(x, lMT, auxdata, isRigidTendon);

% Compute muscle moments
Tmuscs = zeros(size(M));
for dof = 1:size(M,2)
    Tmuscs(:,dof) = sum(FT.*squeeze(B(:,dof,:)), 2);
end

% Remove elements for which no data exists
% Note: this assumes that when matching a specific joint, the
% passive moments created at other joints can be ignored. This is the 
% approach used by Meyer et al. 2016.
M_nz = M(M~=0);
Tmuscs_nz = Tmuscs(M~=0);

% Moments matching objective
Tdiff = M_nz - Tmuscs_nz;

% Return scalar objective
f = sum(Tdiff(:).^2) + 1000*sum((x-1).^2);

end

function [c,ceq] = nonlincon(x, lMT, B, M, auxdata, isRigidTendon)

ceq = [];

auxdata = updateParams(x, auxdata);
params = auxdata.params;

% Muscle-tendon properties
optimalFiberLengthModifier = params(9,:);
tendonSlackLengthModifier = params(10,:);
pennationAngleModifier = params(11,:);
lMo = ones(size(lMT,1),1)*(params(2,:).*optimalFiberLengthModifier);
lTs = ones(size(lMT,1),1)*(params(3,:).*tendonSlackLengthModifier);
alphao = ones(size(lMT,1),1)*(params(4,:).*pennationAngleModifier);

% Hill-type muscle model: geometric relationships
[~, lMtilde] = getMuscleLength(lMT, lTs, lMo, alphao);

% Soft limits on normalized fiber length
softMax_lMtildeCon = log(sum(exp(1000*(lMtilde-1.8))))/1000;
softMin_lMtildeCon = log(sum(exp(1000*(0.1-lMtilde))))/1000;
c = [softMax_lMtildeCon, softMin_lMtildeCon];
c(c<-1000) = -1;
c(c>1000) = 1;

% Limits on peak moments
[FT,~] = getForce(x, lMT, auxdata, isRigidTendon);
Tmuscs = zeros(size(M));
for dof = 1:size(M,2)
    Tmuscs(:,dof) = sum(FT.*squeeze(B(:,dof,:)), 2);
end
M_nz = M(M~=0);
Tmuscs_nz = Tmuscs(M~=0);

[XMAX,IMAX,XMIN,IMIN] = extrema(M_nz);
max_peaks = Tmuscs_nz(IMAX) - XMAX;
min_peaks = XMIN - Tmuscs_nz(IMIN);

c = [c max_peaks' min_peaks'];

end

function [FT, auxdata] = getForce(x, lMT, auxdata, isRigidTendon)

auxdata = updateParams(x, auxdata);

% Passive moment matching objective term
if isRigidTendon
    [FT,~,~,~] = getForceRigidTendon(lMT, auxdata.params);
else
    FT = getForceCompliantTendon(lMT, auxdata.params, auxdata.Fpparam);
end

end

function [auxdata] = updateParams(x, auxdata)

indices = auxdata.indices;
numlMo = indices.numlMo;
numlTs = indices.numlTs;
nume0 = indices.nume0;

% Optimal fiber length modifier
auxdata.params(9,indices.lMo) = x(1:numlMo);

% Tendon slack length 
auxdata.params(10,indices.lTs) = x((numlMo+1):(numlMo+numlTs));

% Muscle strain
m = auxdata.e0LinCoefs.m;
b = auxdata.e0LinCoefs.b;
auxdata.params(7,indices.e0) = m*x((numlMo+numlTs + 1):(numlMo+numlTs+nume0)) + b*ones(nume0,1);
% auxdata.params(7,indices.e0) = m*x + b*ones(nume0,1);


end

function [FT, Fpe, lMtilde, cos_alpha] = getForceRigidTendon(lMT, params)

% Muscle-tendon properties
muscleStrainModifier = ones(size(lMT,1),1)*params(7,:);
muscleShapeFactModifier = ones(size(lMT,1),1)*params(8,:);
optimalFiberLengthModifier = params(9,:);
tendonSlackLengthModifier = params(10,:);
pennationAngleModifier = params(11,:);
FMo = ones(size(lMT,1),1)*params(1,:);
lMo = ones(size(lMT,1),1)*(params(2,:).*optimalFiberLengthModifier);
lTs = ones(size(lMT,1),1)*(params(3,:).*tendonSlackLengthModifier);
alphao = ones(size(lMT,1),1)*(params(4,:).*pennationAngleModifier);

% Hill-type muscle model: geometric relationships
[lM, lMtilde] = getMuscleLength(lMT, lTs, lMo, alphao);

% Passive muscle force-length characteristic
e0 = 0.6*muscleStrainModifier;
kpe = 4*muscleShapeFactModifier;
[~, Fpe] = getPassiveForce(lMtilde, e0, kpe, FMo);

% Pennation angle
cos_alpha = (lMT-lTs)./lM;

% Force along tendon is just passive muscle force times cos(alpha)
FT = Fpe.*cos_alpha;

end

function [FT] = getForceCompliantTendon(lMT, params, Fpparam)

x0 = [ones(size(lMT)) ones(size(lMT))];
lMtilde_min = 0.2; lMtilde_max = 1.8;  
lTtilde_min = 0.8; lTtilde_max = 1.5;
lb = [lMtilde_min*ones(size(lMT)) lTtilde_min*ones(size(lMT))];
ub = [lMtilde_max*ones(size(lMT)) lTtilde_max*ones(size(lMT))];
options = optimoptions('lsqnonlin','Display','iter','MaxFunctionEvaluations',10000);

x = lsqnonlin(@(x) muscleTendonEquilibrium(x, lMT, params, Fpparam), x0, lb, ub, options);

% Get tendon force from solution
numMuscles = size(lMT,2);
lTtilde = x(:,(numMuscles+1):(2*numMuscles));
tendonStiffnessModifier = ones(size(lMT,1),1)*params(6,:);
kT = Fpparam(3).*tendonStiffnessModifier;
FMo = ones(size(lMT,1),1)*params(1,:);
[~, FT] = getTendonForce(lTtilde, kT, FMo);

end

function [f] = muscleTendonEquilibrium(x, lMT, params, Fpparam)

% Muscle-tendon properties
muscleStrainModifier = ones(size(lMT,1),1)*params(7,:);
muscleShapeFactModifier = ones(size(lMT,1),1)*params(8,:);
optimalFiberLengthModifier = params(9,:);
tendonSlackLengthModifier = params(10,:);
FMo = ones(size(lMT,1),1)*params(1,:);
lMo = ones(size(lMT,1),1)*(params(2,:).*optimalFiberLengthModifier);
lTs = ones(size(lMT,1),1)*(params(3,:).*tendonSlackLengthModifier);
tendonStiffnessModifier = ones(size(lMT,1),1)*params(6,:);

% Unpack variables
numMuscles = size(lMT,2);
lMtilde = x(:, 1:numMuscles);
lTtilde = x(:, (numMuscles+1):(2*numMuscles));

% Tendon force-length characteristic
kT = Fpparam(3).*tendonStiffnessModifier;
[~, FT] = getTendonForce(lTtilde, kT, FMo);

% Muscle passive force-length characteristic
e0 = 0.6*muscleStrainModifier;
kpe = 4*muscleShapeFactModifier;
[~, Fpe] = getPassiveForce(lMtilde, e0, kpe, FMo);

% Pennation angle
lM = lMo.*lMtilde;
cos_alpha = (lMT-lTs)./lM;

% Equilibrium conditions
lT = lTs.*lTtilde;
kinEq = lMT - (lM + lT);
forceEq = FT - Fpe.*cos_alpha;
eq = [kinEq forceEq];
f = eq(:);

end

function [lM, lMtilde] = getMuscleLength(lMT, lT, lMo, alphao)

lM = real(sqrt((lMo.*sin(alphao)).^2+(lMT-lT).^2));
lMtilde = lM./lMo;

end

function [fpe, Fpe] = getPassiveForce(lMtilde, e0, kpe, FMo)

t50 = exp(kpe .* (0.2 - 0.10e1) ./ e0);
pp1 = (t50 - 0.10e1); 
t7 = exp(kpe); 
pp2 = (t7 - 0.10e1);

t5 = exp(kpe .* (lMtilde - 0.10e1) ./ e0);
fpe = ((t5 - 0.10e1) - pp1) ./ pp2;
Fpe = FMo.*fpe;

end

function [fT, FT] = getTendonForce(lTtilde, kT, FMo)

fT = (exp(kT.*(lTtilde - 0.995)))/5-0.25;
FT = FMo.*fT;

end

function [B, lMT] = getLengthsMomentArms(model, q, u, coords, muscles)

import org.opensim.modeling.*

% Construct 3D moment arm array
N = size(q, 1);
numCoords = length(coords);
numMuscles = length(muscles);
B = zeros(N, numCoords, numMuscles); 
lMT = zeros(N, numMuscles);
s = model.initSystem();
for i = 1:N
    
    fprintf('Computing muscle lengths and moment arms... (%i/%i)\n',i, N); 
    
    Q = Vector(size(q, 2), 0);
    U = Vector(size(q, 2), 0);
    for j = 1:size(q, 2)
        Q.set(j-1, q(i, j));
        U.set(j-1, u(i, j));
    end
    s.setQ(Q);
    s.setU(U);
    
    % Realize dynamics
    model.computeStateVariableDerivatives(s);

    coordSet = model.getCoordinateSet();
    muscleSet = model.getMuscles();
    for j = 1:numCoords
        coord = coordSet.get(coords{j});
        for k = 1:numMuscles
            muscle = muscleSet.get(muscles{k});
            B(i,j,k) = muscle.computeMomentArm(s, coord);
            if j == 1
                lMT(i,k) = muscle.getLength(s);
            end 
        end
    end
end

end

function [M] = createMomentArray(M_Silder2007)

% Remove high knee flexion and high hip flexion angle from hip moment matching
M_Silder2007{1} = reshape(M_Silder2007{1},numel(M_Silder2007{1})/4,4);
M_Silder2007{1}(:,4) = [];
M_Silder2007{1}(85:96,:) = [];
M_Silder2007{1} = M_Silder2007{1}(:);

% Remove high knee flexion from knee moment matching
M_Silder2007{2} = reshape(M_Silder2007{2},numel(M_Silder2007{2})/4,4);
M_Silder2007{2}(92:121,:) = [];
% M_Silder2007{2}(:,4) = [];
M_Silder2007{2} = M_Silder2007{2}(:);

% Remove high knee flexion from ankle moment matching
M_Silder2007{3} = reshape(M_Silder2007{3},numel(M_Silder2007{3})/4,4);
M_Silder2007{3}(:,1) = [];
M_Silder2007{3} = M_Silder2007{3}(:);

% Combine the passive moments
M_combined = [M_Silder2007{1} zeros(length(M_Silder2007{1}),4);...
              zeros(length(M_Silder2007{2}),2) M_Silder2007{2} zeros(length(M_Silder2007{2}),2);...
              zeros(length(M_Silder2007{3}),3) M_Silder2007{3} zeros(length(M_Silder2007{3}),1)];
      
M = M_combined(:,[1 3 4]);

end

function [q, u] = createJointAngleVelocityArrays(q_Silder2007, allCoords)

% Remove high knee flexion and high hip flexion angle from hip moment matching
q_Silder2007{1} = reshape(q_Silder2007{1},numel(q_Silder2007{1})/(4*6),4,6);
q_Silder2007{1}(:,4,:) = [];
q_Silder2007{1}(85:96,:,:) = [];
q_Silder2007{1} = reshape(q_Silder2007{1},numel(q_Silder2007{1})/(6),6);

% Remove high knee flexion from knee moment matching
q_Silder2007{2} = reshape(q_Silder2007{2},numel(q_Silder2007{2})/(4*6),4,6);
q_Silder2007{2}(92:121,:,:) = [];
% q_Silder2007{2}(:,4,:) = [];
q_Silder2007{2} = reshape(q_Silder2007{2},numel(q_Silder2007{2})/(6),6);

% Remove high knee flexion from ankle moment matching
q_Silder2007{3} = reshape(q_Silder2007{3},numel(q_Silder2007{3})/(4*6),4,6);
q_Silder2007{3}(:,1,:) = [];
q_Silder2007{3} = reshape(q_Silder2007{3},numel(q_Silder2007{3})/(6),6);

% Combine the joint angles
q_combined = [q_Silder2007{1}; q_Silder2007{2}; q_Silder2007{3}]*pi/180;

% Create model states for all coordinates 
zero_q = zeros(size(q_combined,1), 1);
q = zeros(size(q_combined,1), length(allCoords));
for i = 1:length(allCoords)
    switch allCoords{i}
        case 'hip_flexion_r'
            q(:,i) = q_combined(:,1);
        case 'knee_angle_r'
            q(:,i) = q_combined(:,3);
        case 'ankle_angle_r'
            q(:,i) = q_combined(:,4);
        otherwise
            q(:,i) = zero_q;
    end
end

u = zeros(size(q));

end
