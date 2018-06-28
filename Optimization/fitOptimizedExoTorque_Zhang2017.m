function [output] = fitOptimizedExoTorque_Zhang2017(time, q, T_exo, Tmax, startTime, X0)

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
timeShift = linspace(0, t(end)-t(1), length(t)-1);
control = T_exo_norm(:,1);
% plot(timeShift, control)

% Map controls to degrees-of-freedom (which direction is the torque)
mapControl2DOFs = sign(mean(T_exo));

% Problem bounds
N = 4;
lb = zeros(1,N);
ub = ones(1,N);

% peak torque
lb(1) = 0;
ub(1) = 1;

% peak time
lb(2) = 0.05;
ub(2) = 0.95;

% rise time
lb(3) = 0.05;
ub(3) = 0.5;

% fall time
lb(4) = 0.05;
ub(4) = 0.5;

% Initial guess
if isempty(X0)
    x0 = (lb+ub)/2;
else
    x0 = X0;
end

[~, startIdx] = min(abs(timeShift-startTime));

options = optimoptions('fmincon', 'Display', 'iter', ...
                        'OptimalityTolerance', 1e-6, ...
                        'MaxFunctionEvaluations', 100000, ...
                        'MaxIterations', 100000, ...
                        'HessianApproximation', 'lbfgs');
auxdata = struct();
auxdata.control = control;
auxdata.timeShift = timeShift;
auxdata.startTime = startTime;
auxdata.startIdx = startIdx;
auxdata.dqdt = dqdt;
auxdata.mapControl2DOFs = mapControl2DOFs;
auxdata.Popt = Popt;
auxdata.Popt_pos_avg = Popt_pos_avg;
auxdata.Popt_neg_avg = Popt_neg_avg;

[x,f,flag] = fmincon(@(x) fitfunc(x, auxdata), x0, ...
    [], [], [], [], lb, ub, @(x) nonlincon(x, auxdata), options);

[controlFit, T, Y] = buildControlFit(x, auxdata);

torqueFit = controlFit*mapControl2DOFs;
[P, P_pos, P_neg, P_pos_avg, P_neg_avg] = calcDevicePower(torqueFit, dqdt);

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
output.nodes.T = T + time(1);
output.nodes.Y = Y;
output.x = x;
output.solution.peakTorque = x(1)*Tmax;
output.solution.peakTime = output.nodes.T(2);
output.solution.riseTime = output.nodes.T(2)-output.nodes.T(1);
output.solution.fallTime = output.nodes.T(3)-output.nodes.T(2);

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

f = controlRMS + powerRMS;
end

function [fit, T, Y] = buildControlFit(x, auxdata)

startTime = auxdata.startTime;
startIdx = auxdata.startIdx;
timeShift = auxdata.timeShift;
control = auxdata.control;

numPts = length(timeShift(startIdx:end));
numZeros = length(control) - numPts;

% extract variables
peakTorque = x(1);
peakTime = x(2);
riseTime = x(3);
fallTime = x(4);

% set spline nodes based on parameters
node1_t = peakTime - riseTime;
node2_t = peakTime;
node3_t = peakTime + fallTime;

node1_y = 0.01*peakTorque;
node2_y = peakTorque;
node3_y = 0;

% first segment: linear
pre_slope = node1_y/node1_t;
if isinf(pre_slope) || (pre_slope < 0) || isnan(pre_slope) || (pre_slope > 1)
   pre_slope = 0; 
end
x = linspace(0, node1_t, round(node1_t*numPts));
pre_rise = x(1:end-1)*pre_slope;

% rising torque: cubic spline
T = [node1_t node2_t];
Y = [node1_y node2_y];
TT = linspace(node1_t, node2_t, round((node2_t-node1_t)*numPts));
rising_torque = spline(T,[pre_slope Y 0],TT);

% falling torque: cubic spline
T = [node2_t node3_t];
Y = [node2_y node3_y];
TT = linspace(node2_t, node3_t, round((node3_t-node2_t)*numPts));
post_slope = -peakTorque/fallTime;
if isinf(post_slope) || (post_slope < 0) || isnan(post_slope) || (post_slope > 1)
   post_slope = 0; 
end
falling_torque = spline(T,[0 Y post_slope],TT);

% last segment: zeros
post_fall = zeros(1, round((1-node3_t)*numPts));

% concatenate segments
v = [pre_rise rising_torque falling_torque post_fall];

% re-interpolate to ensure the correct number of collocation points
t = linspace(0, 1, length(v));
tq = linspace(0, 1, numPts);

fit = [zeros(1,numZeros) interp1(t, v, tq)]';
T = [node1_t node2_t node3_t]*(timeShift(end)-startTime) + startTime;
Y = [node1_y node2_y node3_y];

end

function [c,ceq] = nonlincon(x, auxdata)

% mapControl2DOFs = auxdata.mapControl2DOFs;
% dqdt = auxdata.dqdt;

controlFit = buildControlFit(x, auxdata);
% T = controlFit*mapControl2DOFs;
% [~, ~, ~, P_pos_avg, P_neg_avg] = calcDevicePower(T, dqdt);

c = max(controlFit) - max(auxdata.control);

% P_net_avg_err = (P_neg_avg + P_pos_avg) - (auxdata.Popt_neg_avg + auxdata.Popt_pos_avg);
% P_neg_avg_err = abs(P_neg_avg - auxdata.Popt_neg_avg);
% P_pos_avg_err = abs(P_pos_avg - auxdata.Popt_pos_avg);
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





