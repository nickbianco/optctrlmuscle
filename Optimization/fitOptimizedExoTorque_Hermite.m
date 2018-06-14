function [output] = fitOptimizedExoTorque_Hermite(time, q, T_exo, Tmax, startTime, N, X0)

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
lb_t = zeros(1,N);
ub_t = timeShift(end)*ones(1,N);
lb_y = zeros(1,N);
ub_y = ones(1,N);
lb = [lb_t lb_y];
ub = [ub_t ub_y];

% Initial guess
if isempty(X0)
    h = 1/N;
    t0 = timeShift(end)*([0:h:(1-h)] + [h:h:1])/2;
    y0 = (lb_y + ub_y)/2;
    x0 = [t0 y0];
else
    x0 = X0;
end

[~, startIdx] = min(abs(timeShift-startTime));

% constraints to ensure time points don't overlap
% t2 > t1, t2 - t1 > 0, t1 - t2 < -0.01
Asq = diag(ones(1,N-1)) + diag(-ones(1,N-2),1);
A = zeros(N-1,N*2);
A(1:(N-1),1:(N-1)) = Asq;
A(N-1,N) = -1;
b = -0.01*ones(N-1,1);

options = optimoptions('fmincon', 'Display', 'iter', ...
                        'OptimalityTolerance', 1e-3, ...
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
while true
    [x,f,flag] = fmincon(@(x) fitfunc(x, auxdata), x0, ...
        A, b, [], [], lb, ub, @(x) nonlincon(x, auxdata), options);
    
    if flag>0
        break;
    else
        t0 = timeShift(end)*([0:h:(1-h)] + [h:h:1])/2;
        y0 = rand(1,N);
        x0 = [t0 y0];
    end
    
end

[controlFit, controlFit_prime, T, Y] = buildControlFit(x, auxdata);

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

end

function [f] = fitfunc(x, auxdata)

control = auxdata.control;
dqdt = auxdata.dqdt;
mapControl2DOFs = auxdata.mapControl2DOFs;

[controlFit, controlFit_prime] = buildControlFit(x, auxdata);

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

function [fit, fit_prime, T, Y] = buildControlFit(x, auxdata)

startTime = auxdata.startTime;
startIdx = auxdata.startIdx;
timeShift = auxdata.timeShift;
control = auxdata.control;

N = length(x)/2;

Tscale = startTime + (1-startTime)*x(1:N);
T = [startTime Tscale];
Y = [0 x((N+1):end)];
fit = zeros(size(control));
pp = pchip(T,Y);
fit(startIdx:end) = ppval(pp,timeShift(startIdx:end));
fit(fit<0) = 0;
fit_prime = ppval(fnder(pp,1),timeShift(startIdx:end));

end

function [c,ceq] = nonlincon(x, auxdata)

mapControl2DOFs = auxdata.mapControl2DOFs;
dqdt = auxdata.dqdt;

controlFit = buildControlFit(x, auxdata);
T = controlFit*mapControl2DOFs;
[~, ~, ~, P_pos_avg, P_neg_avg] = calcDevicePower(T, dqdt);

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





