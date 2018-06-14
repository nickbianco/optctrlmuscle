% Hill-type muscle model: equilibrium between muscle and tendon forces
% All muscle-tendon characteristics are fully described in the publication
% and its online supplement

function [muscleData] = DeGroote2016Muscle_lMtildeState(a,lMtilde,vMtilde,lMT,params,Fvparam,Fpparam,Faparam)

FMo = ones(size(a,1),1)*params(1,:);
lMo = ones(size(a,1),1)*params(2,:);
lTs = ones(size(a,1),1)*params(3,:);
alphao = ones(size(a,1),1)*params(4,:);
vMmax = ones(size(a,1),1)*params(5,:);
tendonStiffnessModifier = ones(size(a,1),1)*params(6,:);
muscleStrainModifier = ones(size(a,1),1)*params(7,:);
muscleShapeFactModifier = ones(size(a,1),1)*params(8,:);

% Hill-type muscle model: geometric relationships
lM = lMtilde.*lMo;
w = lMo.*sin(alphao);
lT = lMT - real(sqrt((lM.^2 - w.^2)));
lTtilde = lT./lTs;

% Tendon force-length characteristic
tendonStiffness = Fpparam(3).*tendonStiffnessModifier;
fse = (exp(tendonStiffness.*(lTtilde - 0.995)))/5-0.25;

% Active muscle force-length characteristic
b11 = Faparam(1);
b21 = Faparam(2);
b31 = Faparam(3);
b41 = Faparam(4);
b12 = Faparam(5);
b22 = Faparam(6);
b32 = Faparam(7);
b42 = Faparam(8);

b13 = 0.1;
b23 = 1;
b33 = 0.5*sqrt(0.5);
b43 = 0;
num3 = lMtilde-b23;
den3 = b33+b43*lMtilde;
FMtilde3 = b13*exp(-0.5*num3.^2./den3.^2);

num1 = lMtilde-b21;
den1 = b31+b41*lMtilde;
FMtilde1 = b11*exp(-0.5*num1.^2./den1.^2);

num2 = lMtilde-b22;
den2 = b32+b42*lMtilde;
FMtilde2 = b12*exp(-0.5*num2.^2./den2.^2);

FMltilde = FMtilde1+FMtilde2+FMtilde3;

% Active muscle force-velocity characteristic
e1 = Fvparam(1);
e2 = Fvparam(2);
e3 = Fvparam(3);
e4 = Fvparam(4);
FMvtilde = e1*log((e2*vMtilde+e3)+sqrt((e2*vMtilde+e3).^2+1))+e4;

% Active muscle force
fce = a.*FMltilde.*FMvtilde;

% Passive muscle force-length characteristic
e0 = 0.6*muscleStrainModifier;
kpe = 4*muscleShapeFactModifier;
t5 = exp(kpe .* (lMtilde - 0.10e1) ./ e0);
fpe = ((t5 - 0.10e1) - Fpparam(1,:)) ./ Fpparam(2,:);

% Muscle force
Fce = FMo.*fce;
Fpe = FMo.*fpe;
FM = Fce + Fpe;

% Tendon force
FT = fse .* FMo;

% Equilibrium between muscle and tendon forces
% Fm*cos(alpha) = Ft
cos_alpha = (lMT-lT)./lM;
err =  FM.*cos_alpha-FT;

% Set outputs
muscleData.err = err;
muscleData.fce = fce;
muscleData.Fce = Fce;
muscleData.fpe = fpe;
muscleData.Fpe = Fpe;
muscleData.FM = FM;
muscleData.lTtilde = lTtilde;
muscleData.lT = lT;
muscleData.fse = fse;
muscleData.FT = FT;
muscleData.FMltilde = FMltilde;
muscleData.lMtilde = lMtilde;
muscleData.lM = lM;
muscleData.vMtilde = vMtilde;
muscleData.vM = vMmax.*vMtilde;
muscleData.cos_alpha = cos_alpha;

end
