
% This function computes the muscle fiber length from the normalized tendon
% force

function [lM,lMtilde,vM,vMtilde] = FiberLengthVelocity_Ftilde(Ftilde,dFtilde,params,lMT,vMT,Fpparam)

lMo = ones(size(Ftilde,1),1)*params(2,:);
lTs = ones(size(Ftilde,1),1)*params(3,:);
alphao = ones(size(Ftilde,1),1)*params(4,:);
vMmax = ones(size(Ftilde,1),1)*params(5,:);
tendonStiffnessModifier = ones(size(Ftilde,1),1)*params(6,:);
tendonStiffness = Fpparam(3)*tendonStiffnessModifier;

% Non-linear tendon
lTtilde = real(log(5*(Ftilde + 0.25))./tendonStiffness + 0.995);

% Hill-model relationship
lM = sqrt((lMo.*sin(alphao)).^2+(lMT-lTs.*lTtilde).^2);
lMtilde = lM./lMo;

% Tendon velocity
vT = lTs.*dFtilde./(7*exp(tendonStiffness.*(lTtilde-0.995)));

% Muscle velocity
cos_alpha = (lMT-lTs.*lTtilde)./lM;
vM = (vMT-vT).*cos_alpha;
vMtilde = vM./vMmax;

end
