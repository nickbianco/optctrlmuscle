function Edot = calcMinettiAlexanderProbe(vMtilde,vmax,Fo,a)
% INPUTS:
%  vMtilde: normalized shortening speed 
%     vmax: max shortening speed [m/s]
%       Fo: max isometric force [N]
%        a: activation 

% OUTPUTS:
%     Edot: metabolic rate

% Minetti & Alexander (1997) model parameters
c1 = 0.054;
c2 = 0.506;
c3 = 2.46;
c4 = 1.13;
c5 = 12.8;
c6 = 1.64;

phi = (c1*ones(size(vMtilde)) + c2*vMtilde + c3*(vMtilde.^2))./ ... 
      (ones(size(vMtilde)) - c4*vMtilde + c5*(vMtilde.^2) - c6*(vMtilde.^3)); 

% Total metabolic rate
Edot = Fo*vmax*a.*phi;

