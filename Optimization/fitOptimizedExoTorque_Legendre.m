function [output] = fitOptimizedExoTorque_Legendre(time, q, T_exo, Tmax, startTime, N, orders)

% Calculate angular velocity
q = (pi/180)*q; % convert to radians
dq = diff(q);
t = linspace(time(1), time(end), size(q,1));
dt = diff(t)';
dqdt = dq./dt; % rad/s

% Calculate normalized device torque
T_exo_norm = T_exo / Tmax; % unitless
T_exo_norm = interp1(time, T_exo_norm, t);
T_exo_norm(end,:) = [];

% Calculate "normalized" positive power
[Popt, Popt_pos, Popt_neg, Popt_pos_avg, Popt_neg_avg] = calcDevicePower(T_exo_norm, dqdt);

% Shift time to start at zero and plot optimized device control profile
startTime = startTime - time(1);
time = linspace(0, t(end)-t(1), length(t)-1);
control = T_exo_norm(:,1);
plot(time, control)

% Map controls to degrees-of-freedom (which direction is the torque)
mapControl2DOFs = sign(mean(T_exo));

% Problem bounds
lb = -1*ones(1,N);
ub = 1*ones(1,N);

% Initial guess
x0 = (lb+ub)/2;
    
% Pre-construct the Legendre polynomial basis
[~, startIdx] = min(abs(time-startTime));
tpoly = linspace(-1, 1, length(time(startIdx:end)));
L = zeros(length(tpoly), N);
for i = 1:N
    L(:,i) = legendreP(i-1,tpoly);
end

ordersIdx = ones(1,N);
if ~isempty(orders)    
    ordersIdx = zeros(1,N);
    ordersIdx(orders+1) = 1;
end

options = optimoptions('fmincon', 'Display', 'iter', ...
                        'OptimalityTolerance', 1e-6, ...
                        'MaxFunctionEvaluations', 100000, ...
                        'MaxIterations', 100000, ...
                        'HessianApproximation', 'lbfgs');
auxdata = struct();
auxdata.control = control;
auxdata.time = time;
auxdata.startTime = startTime;
auxdata.startIdx = startIdx;
auxdata.dqdt = dqdt;
auxdata.mapControl2DOFs = mapControl2DOFs;
auxdata.Popt = Popt;
auxdata.Popt_pos_avg = Popt_pos_avg;
auxdata.Popt_neg_avg = Popt_neg_avg;
auxdata.L = L;
auxdata.ordersIdx = ordersIdx;
x = fmincon(@(x) fitfunc(x, auxdata), x0, ...
              [], [], [], [], lb, ub, @(x) nonlincon(x, auxdata), options);
          
controlFit = buildControlFit(x, auxdata);

T = controlFit*mapControl2DOFs;
[P, P_pos, P_neg, P_pos_avg, P_neg_avg] = calcDevicePower(T, dqdt);

% if isempty(orders)
%     coefCutoffs = 0.01:0.002:0.05;
%     for i = 1:length(coefCutoffs)
%         fprintf(['Legendre polynomial orders with optimized coefficients \n' ...
%             'magnitudes greater than %0.3f: \n'], coefCutoffs(i))
%         orders = find((abs(x)>coefCutoffs(i)))-1
%     end
% end

output.time = linspace(t(1), t(end), length(t)-1);
output.opt.control = control;
output.opt.P = Popt;
output.opt.P_pos = Popt_pos;
output.opt.P_neg = Popt_neg;
output.opt.P_pos_avg = Popt_pos_avg;
output.opt.P_neg_avg = Popt_neg_avg;

output.fit.control = controlFit;
output.fit.P = P;
output.fit.P_pos = P_pos;
output.fit.P_neg = P_neg;
output.fit.P_pos_avg = P_pos_avg;
output.fit.P_neg_avg = P_neg_avg;
output.x = x;

end

function [f] = fitfunc(x, auxdata)

control = auxdata.control;
dqdt = auxdata.dqdt;
mapControl2DOFs = auxdata.mapControl2DOFs;

controlFit = buildControlFit(x, auxdata);

torqueFit = controlFit*mapControl2DOFs;
[P, P_pos, P_neg, P_pos_avg, P_neg_avg] = calcDevicePower(torqueFit, dqdt);

% corrVal = corr2(controlFit, control);
% if isnan(corrVal)
%     corrVal = -1; 
% end
% corrCost = 1 - corrVal;

controlRMS = rms(controlFit - control) / max(control);
powerRMS = sum(rms(P - auxdata.Popt)) / sum(max(auxdata.Popt));

f = controlRMS;

end

function [fit] = buildControlFit(x, auxdata)

startIdx = auxdata.startIdx;
ordersIdx = auxdata.ordersIdx;
L = auxdata.L;
time = auxdata.time;
control = auxdata.control;

t = linspace(-1, 1, length(time(startIdx:end)));

% linCombLegendre = zeros(size(t));
% for i = 1:length(x)
%    n = orders(i);
%    linCombLegendre = linCombLegendre + x(i)*legendre(n, t);
% end

x(~ordersIdx) = 0;
linCombLegendre = L*x';

fit = zeros(size(control));
fit(startIdx:end) = linCombLegendre;
fit(fit<0) = 0;

end

function [c,ceq] = nonlincon(x, auxdata)

% mapControl2DOFs = auxdata.mapControl2DOFs;
% dqdt = auxdata.dqdt;

controlFit = buildControlFit(x, auxdata);
% T = controlFit*mapControl2DOFs;
% [~, ~, ~, P_pos_avg, P_neg_avg] = calcDevicePower(T, dqdt);

c = max(controlFit) - max(auxdata.control);
% ceq = [P_neg_avg - auxdata.Popt_neg_avg;
%        P_pos_avg - auxdata.Popt_pos_avg];
ceq = [];

end

function [P, P_pos, P_neg, P_pos_avg, P_neg_avg] = calcDevicePower(T, dqdt)

P = T.*dqdt;
P_pos = P;
P_neg = P;
P_pos(P_pos < 0) = 0;
P_neg(P_neg > 0) = 0;
P_pos_avg = sum(mean(P_pos));
P_neg_avg = sum(mean(P_neg));

end
