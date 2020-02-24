function [output] = fitOptimizedExoTorque_Zhang2017(time, q, T_exo, Tmax_hip, Tmax_knee, Tmax_ankle, startTime, X0, mod_name)

%% Compute helpful terms for data fitting

% Calculate angular velocity
q = (pi/180)*q; % convert to radians
dq = diff(q);
t = linspace(time(1), time(end), size(q,1));
dt = diff(t)';
dqdt = dq./dt; % rad/s

% Calculate normalized device torque
T_exo_norm = T_exo;
T_exo_norm(:,1) = T_exo(:,1) / Tmax_hip;
T_exo_norm(:,2) = T_exo(:,2) / Tmax_knee; 
T_exo_norm(:,3) = T_exo(:,3) / Tmax_ankle; 
T_exo_norm = interp1(time, T_exo_norm, t);
T_exo_fit = T_exo_norm;

% Calculate "normalized" positive power
[Popt, Popt_pos, Popt_neg, Popt_pos_avg, Popt_neg_avg] = calcDevicePower(T_exo_norm, dqdt);
Popt_fit = Popt;

% Shift time to start at zero
startTime = startTime - time(1);
timeShift = linspace(0, t(end)-t(1), length(t)-1);

% Map controls to degrees-of-freedom (which direction is the torque)
mapControl2DOFs = sign(mean(T_exo));

%% Bounds & Initial Guesses
config = ReadYaml('C:\Users\Nick\Projects\ExoTopology\exotopology\config.yaml');
peak_torque = config.param_bounds.peak_torque;
peak_time = config.param_bounds.peak_time;
ext_peak_time = config.param_bounds.ext_peak_time;
rise_time = config.param_bounds.rise_time;
fall_time = config.param_bounds.fall_time;

N = 4;
lb = zeros(1,N);
ub = ones(1,N);

% peak torque
lb(1) = peak_torque{1};
ub(1) = peak_torque{2};

% peak time
lb(2) = peak_time{1};
ub(2) = peak_time{2};

% rise time
lb(3) = rise_time{1};
ub(3) = rise_time{2};

% fall time
lb(4) = fall_time{1};
ub(4) = fall_time{2};

% Take absolute value of exo torque (since it only acts in one direction) and
% find max values and their times.
T_exo_fit_abs = abs(T_exo_fit);
[maxVals, idx] = max(T_exo_fit_abs);
% Create percent gait cycle vector based on length of exo torque vector
perGC = linspace(0,1,size(T_exo_fit,1));

% Active hip flexion assistance
if strcmp(mod_name, 'mrsmod_actHf')
      
    X0 = (lb+ub)/2;
    X0(1) = maxVals(1); % peak torque guess based on max hip exo torque value
    X0(2) = perGC(idx(1)); % peak time guess based on max hip exo torque time
    
% Active knee flexion assistance
elseif strcmp(mod_name,'mrsmod_actKf')
    
    % TODO
    [~, locs] = findpeaks(T_exo_fit_abs(:,2));
    midIdx = round(size(T_exo_fit_abs,1)/2);
    [~,idx] = min(abs(locs - midIdx));
    midPeakLoc = locs(idx);
    X0 = (lb+ub)/2;
    X0(1) = T_exo_fit_abs(midPeakLoc,2);
    X0(2) = perGC(midPeakLoc);

% Active ankle plantarflexion
elseif strcmp(mod_name, 'mrsmod_actAp')
  
    X0 = (lb+ub)/2;
    X0(1) = maxVals(3); % peak torque guess based on max ankle exo torque value
    X0(2) = perGC(idx(3)); % peak time guess based on max ankle exo torque time

% Active coupled hip flexion + ankle plantarflexion assistance
elseif strcmp(mod_name,'mrsmod_actHfAp')
    % ankle peak torque
    lb(5) = peak_torque{1};
    ub(5) = peak_torque{2};
    
    X0 = (lb+ub)/2;
    X0(1) = maxVals(1); % peak torque guess based on max hip exo torque value
    X0(2) = perGC(idx(1)); % peak time guess based on max hip exo torque time
    X0(5) = maxVals(3); 
    
% Active coupled hip flexion + knee flexion assistance
elseif strcmp(mod_name,'mrsmod_actHfKf')
    % knee peak torque
    lb(5) = peak_torque{1};
    ub(5) = peak_torque{2};
    
    X0 = (lb+ub)/2;
    X0(1) = maxVals(1); % peak torque guess based on max hip exo torque value
    X0(2) = perGC(idx(1)); % peak time guess based on max hip exo torque time
    X0(5) = maxVals(2); 

% Active coupled knee flexion + ankle plantarflexion assistance
elseif strcmp(mod_name,'mrsmod_actKfAp')
    % ankle peak torque
    lb(5) = peak_torque{1};
    ub(5) = peak_torque{2};
    
    X0 = (lb+ub)/2;
    X0(1) = maxVals(2); % peak torque guess based on max knee exo torque value
    X0(2) = perGC(idx(2)); % peak time guess based on max knee exo torque time
    X0(5) = maxVals(3); 

% Active coupled hip flexion + knee flexion + ankle plantarflexion assistance
elseif strcmp(mod_name,'mrsmod_actHfKfAp')
    % knee peak torque
    lb(5) = peak_torque{1};
    ub(5) = peak_torque{2};
    
    % ankle peak torque
    lb(6) = peak_torque{1};
    ub(6) = peak_torque{2};
    
    X0 = (lb+ub)/2;
    X0(1) = maxVals(1); % peak torque guess based on max hip exo torque value
    X0(2) = perGC(idx(1)); % peak time guess based on max hip exo torque time
    X0(5) = maxVals(2); 
    X0(6) = maxVals(3);
    
% Active hip extension assistance
elseif contains(mod_name,'actHe_')
    
    % peak time
    lb(2) = ext_peak_time{1};
    ub(2) = ext_peak_time{2};

    % Find the first peak in hip extension torque, then find the next point in
    % in time where torque hits zero (within tolerance) and set torque to zero
    % at every time point after that. This is so we fit to the main hip 
    % extension assistive torque peak and ignore the rest, which we couldn't 
    % fit with the current parameterization anyway.
    hipFlag = false;
    for i = 1:size(T_exo_fit_abs,1)
        if ((i > idx(1)) && (T_exo_fit_abs(i,1) < 0.001)) || hipFlag
            T_exo_fit(i,1) = 0;
            hipFlag = true;
        end      
    end
%     Popt_fit = calcDevicePower(T_exo_fit, dqdt);
    
    perGC = linspace(0,1,size(T_exo_fit,1));
    X0 = (lb+ub)/2;
    X0(1) = maxVals(1); % peak torque guess based on max hip exo torque value
    X0(2) = perGC(idx(1)); % peak time guess based on max hip exo torque time
    % Set rise time guess based on peak time, since torque usually begins
    % immediately after heel strike.
    if X0(2) < ub(3)
        X0(3) = X0(2);
        ub(3) = X0(2)+0.05;
    end
    % Set fall time guess based on when the torque hit zero after peak time.
    X0(4) = perGC(nnz(abs(T_exo_fit(:,1))))-perGC(idx(1));
    
% Active knee extension assistance
elseif contains(mod_name,'actKe_')
    
    % peak time
    lb(2) = ext_peak_time{1};
    ub(2) = ext_peak_time{2};
    
    % Find the first peak in knee extension torque, then find the next point in
    % in time where torque hits zero (within tolerance) and set torque to zero
    % at every time point after that. This is so we fit to the main knee 
    % extension assistive torque peak and ignore the rest, which we couldn't 
    % fit with the current parameterization anyway.
    kneeFlag = false;
    for i = 1:size(T_exo_fit_abs,1)
        if (i > idx(2)) && (T_exo_fit_abs(i,2) < 0.001) || kneeFlag
            T_exo_fit(i,2) = 0;
            kneeFlag = true;
        end       
    end
    % Create new power trajectory to fit based on adjusted torque profile.
%     Popt_fit = calcDevicePower(T_exo_fit, dqdt);
    
    perGC = linspace(0,1,size(T_exo_fit,1));
    X0 = (lb+ub)/2;
    X0(1) = maxVals(2); % peak torque guess based on max knee exo torque value
    X0(2) = perGC(idx(2)); % peak time guess based on max knee exo torque time
    % Set rise time guess based on peak time, since torque usually begins
    % immediately after heel strike.
    if X0(2) < ub(3)
        X0(3) = X0(2);
        ub(3) = X0(2)+0.05;
    end
    % Set fall time guess based on when the torque hit zero after peak time.
    X0(4) = perGC(nnz(abs(T_exo_fit(:,2))))-perGC(idx(2));

% Active coupled hip extension + knee extension assistance
elseif strcmp(mod_name,'mrsmod_actHeKe')
    
    % peak time
    lb(2) = ext_peak_time{1};
    ub(2) = ext_peak_time{2};
    
    % knee peak torque
    lb(5) = peak_torque{1};
    ub(5) = peak_torque{2};
    
    % Find the first peak in both the knee extension and hip extension torques, 
    % then find the next point in time where each torque hits zero (within 
    % tolerance) and set torques to zero at every time point after that. This is 
    % so we fit to the main knee extension and hip extension assistive torque 
    % peaks and ignore the rest, which we couldn't fit with the current 
    % parameterization anyway.
    hipFlag = false;
    kneeFlag = false;
    for i = 1:size(T_exo_fit_abs,1)
        if ((i > idx(1)) && (T_exo_fit_abs(i,1) < 0.001)) || hipFlag
            T_exo_fit(i,1) = 0;
            hipFlag = true;
        end
        if (i > idx(2)) && (T_exo_fit_abs(i,2) < 0.001) || kneeFlag
            T_exo_fit(i,2) = 0;
            kneeFlag = true;
        end       
    end
    % Create new power trajectory to fit based on adjusted torque profile.
%     Popt_fit = calcDevicePower(T_exo_fit, dqdt);
    
    X0 = (lb+ub)/2;
    X0(1) = maxVals(1); % peak torque guess based on max hip exo torque value
    X0(2) = perGC(idx(1)); % peak time guess based on max hip exo torque time
    % Set rise time guess based on peak time, since torque usually begins
    % immediately after heel strike.
    if X0(2) < ub(3)
        X0(3) = X0(2);
        ub(3) = X0(2)+0.05;
    end
    % Set fall time guess based on average of when knee and hip torques hit zero 
    % after their respective peak times.
    X0(4) = (perGC(nnz(abs(T_exo_fit(:,1))))-perGC(idx(1)) + ... 
             perGC(nnz(abs(T_exo_fit(:,2))))-perGC(idx(2)))/2;
    X0(5) = maxVals(2);
   
elseif strcmp(mod_name, 'mrsmod_actHfKfAp_multControls')    
    
    [peakVals, peakLocs] = findpeaks(T_exo_fit_abs(:,2));
    [peakValsSort, idxSort] = sort(peakVals, 'descend');
    peakLocsSort = peakLocs(idxSort);
    
    % There are often two peaks in the knee assistive torque curve. We would 
    % prefer to fit to the one nearest to the hip & ankle peak moments near the
    % middle of the gait cycle.
    if abs(peakLocsSort(1) - 50) > abs(peakLocsSort(2) - 50)
        kneePeakTorqueGuess = peakValsSort(2);
        kneePeakTimeGuess = perGC(peakLocsSort(2));
    else
        kneePeakTorqueGuess = maxVals(2);
        kneePeakTimeGuess = perGC(idx(2));
    end
    
    % Create new power trajectory to fit based on adjusted torque profile.
%     Popt_fit = calcDevicePower(T_exo_fit, dqdt);
        
    % peak torque
    lb(5) = peak_torque{1};
    ub(5) = peak_torque{2};
    lb(9) = peak_torque{1};
    ub(9) = peak_torque{2};

    % peak time
    lb(6) = peak_time{1};
    ub(6) = peak_time{2};
    lb(10) = peak_time{1};
    ub(10) = peak_time{2};

    % rise time
    lb(7) = rise_time{1};
    ub(7) = rise_time{2};
    lb(11) = rise_time{1};
    ub(11) = rise_time{2};

    % fall time
    lb(8) = fall_time{1};
    ub(8) = fall_time{2};
    lb(12) = fall_time{1};
    ub(12) = fall_time{2};
    
    % initial guess
    X0 = (lb+ub)/2;
    X0(1) = maxVals(1); % peak torque guess based on max hip exo torque value
    X0(2) = perGC(idx(1)); % peak time guess based on max hip exo torque time
    
    X0(5) = kneePeakTorqueGuess; % peak torque guess based on max knee exo torque value
    X0(6) = kneePeakTimeGuess; % peak time guess based on max *hip* exo torque time
    
    X0(9) = maxVals(3); % peak torque guess based on max ankle exo torque value
    X0(10) = perGC(idx(3)); % peak time guess based on max ankle exo torque time
    
elseif strcmp(mod_name, 'mrsmod_actHfKf_multControls')
    
    [peakVals, peakLocs] = findpeaks(T_exo_fit_abs(:,2));
    [peakValsSort, idxSort] = sort(peakVals, 'descend');
    peakLocsSort = peakLocs(idxSort);
    
    % There are often two peaks in the knee assistive torque curve. We would 
    % prefer to fit to the one nearest to the hip & ankle peak moments near the
    % middle of the gait cycle.
    if abs(peakLocsSort(1) - 50) > abs(peakLocsSort(2) - 50)
        kneePeakTorqueGuess = peakValsSort(2);
        kneePeakTimeGuess = perGC(peakLocsSort(2));
    else
        kneePeakTorqueGuess = maxVals(2);
        kneePeakTimeGuess = perGC(idx(2));
    end
    % Create new power trajectory to fit based on adjusted torque profile.
%     Popt_fit = calcDevicePower(T_exo_fit, dqdt);
        
    % peak torque
    lb(5) = peak_torque{1};
    ub(5) = peak_torque{2};

    % peak time
    lb(6) = peak_time{1};
    ub(6) = peak_time{2};

    % rise time
    lb(7) = rise_time{1};
    ub(7) = rise_time{2};

    % fall time
    lb(8) = fall_time{1};
    ub(8) = fall_time{2};
    
    % initial guess
    X0 = (lb+ub)/2;
    X0(1) = maxVals(1); % peak torque guess based on max hip exo torque value
    X0(2) = perGC(idx(1)); % peak time guess based on max hip exo torque time
    
    X0(5) = kneePeakTorqueGuess; % peak torque guess based on max knee exo torque value
    X0(6) = kneePeakTimeGuess; % peak time guess based on max *hip* exo torque time
    
elseif strcmp(mod_name, 'mrsmod_actHfAp_multControls')
    % peak torque
    lb(5) = peak_torque{1};
    ub(5) = peak_torque{2};

    % peak time
    lb(6) = peak_time{1};
    ub(6) = peak_time{2};

    % rise time
    lb(7) = rise_time{1};
    ub(7) = rise_time{2};

    % fall time
    lb(8) = fall_time{1};
    ub(8) = fall_time{2};
    
    % initial guess
    X0 = (lb+ub)/2;
    X0(1) = maxVals(1); % peak torque guess based on max hip exo torque value
    X0(2) = perGC(idx(1)); % peak time guess based on max hip exo torque time
    X0(5) = maxVals(3); % peak torque guess based on max ankle exo torque value
    X0(6) = perGC(idx(3)); % peak time guess based on max ankle exo torque time
    
elseif strcmp(mod_name, 'mrsmod_actKfAp_multControls')
    % peak torque
    lb(5) = peak_torque{1};
    ub(5) = peak_torque{2};

    % peak time
    lb(6) = peak_time{1};
    ub(6) = peak_time{2};

    % rise time
    lb(7) = rise_time{1};
    ub(7) = rise_time{2};

    % fall time
    lb(8) = fall_time{1};
    ub(8) = fall_time{2};
    
    % initial guess
    X0 = (lb+ub)/2;
    X0(1) = maxVals(2); % peak torque guess based on max knee exo torque value
    X0(2) = perGC(idx(2)); % peak time guess based on max knee exo torque time
    X0(5) = maxVals(3); % peak torque guess based on max ankle exo torque value
    X0(6) = perGC(idx(3)); % peak time guess based on max ankle exo torque time
    
elseif strcmp(mod_name, 'mrsmod_actHeKe_multControls')
     % peak torque
    lb(5) = peak_torque{1};
    ub(5) = peak_torque{2};

    % peak time
    lb(2) = ext_peak_time{1};
    ub(2) = ext_peak_time{2};
    lb(6) = ext_peak_time{1};
    ub(6) = ext_peak_time{2};

    % rise time
    lb(7) = rise_time{1};
    ub(7) = rise_time{2};

    % fall time
    lb(8) = fall_time{1};
    ub(8) = fall_time{2};
        
    % Find the first peak in both the knee extension and hip extension torques, 
    % then find the next point in time where each torque hits zero (within 
    % tolerance) and set torques to zero at every time point after that. This is 
    % so we fit to the main knee extension and hip extension assistive torque 
    % peaks and ignore the rest, which we couldn't fit with the current 
    % parameterization anyway.
    hipFlag = false;
    kneeFlag = false;
    for i = 1:size(T_exo_fit_abs,1)
        if ((i > idx(1)) && (T_exo_fit_abs(i,1) < 0.001)) || hipFlag
            T_exo_fit(i,1) = 0;
            hipFlag = true;
        end
        if (i > idx(2)) && (T_exo_fit_abs(i,2) < 0.001) || kneeFlag
            T_exo_fit(i,2) = 0;
            kneeFlag = true;
        end       
    end
    % Create new power trajectory to fit based on adjusted torque profile.
    Popt_fit = calcDevicePower(T_exo_fit, dqdt);
    
    X0 = (lb+ub)/2;
    X0(1) = maxVals(1); % peak torque guess based on max hip exo torque value
    X0(2) = perGC(idx(1)); % peak time guess based on max hip exo torque time
    X0(5) = maxVals(2); % peak torque guess based on max knee exo torque value
    X0(6) = perGC(idx(2)); % peak time guess based on max knee exo torque time
    
    % Set rise time guess based on peak time, since torque usually begins
    % immediately after heel strike.
    if X0(2) < ub(3)
        X0(3) = X0(2);
        ub(3) = X0(2)+0.05;
    end
    if X0(6) < ub(7)
        X0(7) = X0(6);
        ub(7) = X0(6)+0.05;
    end
    % Set fall time guess based on when knee and hip torques hit zero 
    % after their respective peak times.
    X0(4) = perGC(nnz(abs(T_exo_fit(:,1))))-perGC(idx(1)); 
    X0(8) = perGC(nnz(abs(T_exo_fit(:,2))))-perGC(idx(2));
        
end

%% Problem setup

% Auxiliary data available to each data fitting problem
[~, startIdx] = min(abs(timeShift-startTime));
auxdata = struct();
auxdata.mod_name = mod_name;
auxdata.timeShift = timeShift;
auxdata.startTime = startTime;
auxdata.startIdx = startIdx;
auxdata.dqdt = dqdt;
auxdata.mapControl2DOFs = mapControl2DOFs;
auxdata.Popt = Popt;
auxdata.Popt_pos_avg = Popt_pos_avg;
auxdata.Popt_neg_avg = Popt_neg_avg;
auxdata.T_exo_norm = T_exo_norm;
auxdata.T_exo_fit = T_exo_fit;
auxdata.Popt_fit = Popt_fit;
auxdata.lb = lb;
auxdata.ub = ub;

% Set initial guess structure to pass to problem
if isempty(X0)
    x0 = (lb+ub)/2;
else
    x0 = scaleVariables(X0, auxdata);
end

% LSQNONLIN options
options = optimoptions('lsqnonlin', 'Display', 'iter', ...
                        'OptimalityTolerance', 1e-9, ...
                        'MaxFunctionEvaluations', 100000, ...
                        'MaxIterations', 100000);

% Run problem ten times and slightly perturb each solution to create an initial
% guess for the next problem. Take best solution of the ten.
x0_redo = x0;
a = 1.0;
b = -1.0;
resnorm_best = inf;
x_best = x0;
for i = 1:25    
    [x,resnorm] = lsqnonlin(@(x) fitfunc_lsq(x, auxdata), x0_redo, -ones(size(lb)), ...
        ones(size(ub)), options);
    if resnorm < resnorm_best
       x_best = x; 
       resnorm_best = resnorm;
    end
    perturb = a + (b-a).*rand(size(x_best));
    x0_redo = x_best+perturb;
    a = a*0.9;
    b = b*0.9;
end

%% Store solution in output struct
torque = calcDeviceTorque(x_best, auxdata);
[P, P_pos, P_neg, P_pos_avg, P_neg_avg] = calcDevicePower(torque, dqdt);

output.time = linspace(t(1), t(end), length(t)-1);
output.opt.P = Popt;
output.opt.P_pos = Popt_pos;
output.opt.P_neg = Popt_neg;
output.opt.P_pos_avg = Popt_pos_avg;
output.opt.P_neg_avg = Popt_neg_avg;
output.opt.torque = T_exo_norm(1:(end-1),:);

output.fit.P = P;
output.fit.P_pos = P_pos;
output.fit.P_neg = P_neg;
output.fit.P_pos_avg = P_pos_avg;
output.fit.P_neg_avg = P_neg_avg;
output.fit.torque = torque(1:(end-1),:);
output.x = unscaleVariables(x_best, auxdata);
output.lb = lb;
output.ub = ub;

end

%% SCALEVARIABLES
function [x] = scaleVariables(x, auxdata)
    x = ((2*x - auxdata.lb)./(auxdata.ub-auxdata.lb))-1;
end

%% UNSCALEVARIABLES
function [x] = unscaleVariables(x, auxdata)
    x = 0.5*(auxdata.ub-auxdata.lb).*(x+1) + auxdata.lb;
end

%% FITFUNC_LSQ
function [f] = fitfunc_lsq(x, auxdata)

% torque errors
torque = calcDeviceTorque(x, auxdata);
% perGC_torque = linspace(0, 1, size(torque, 1));
% [~, maxIdx_torque] = max(abs(torque));
% [~, maxIdx_torqueFit] = max(abs(auxdata.T_exo_fit));
% peakTimeErr_torque = abs(perGC_torque(maxIdx_torque) - perGC_torque(maxIdx_torqueFit));
torqueErr = (abs(torque) - abs(auxdata.T_exo_fit)); % ./ sum(max(abs(auxdata.T_exo_fit)));
% peakValErr_torque = torqueErr(maxIdx_torqueFit);
% torqueErrPos = torqueErr;
% torqueErrPos(torqueErrPos<0) = 0;

% power errors
[P, P_pos, ~, ~, ~] = calcDevicePower(torque, auxdata.dqdt);
% perGC_P = linspace(0, 1, size(P, 1));
% [~, maxIdx_P] = max(abs(P));
% [~, maxIdx_Pfit] = max(abs(auxdata.Popt_fit));
% peakTimeErr_P = abs(perGC_P(maxIdx_P) - perGC_P(maxIdx_Pfit));
powerErr = (P - auxdata.Popt_fit); % / sum(max(auxdata.Popt_fit));

% f = [powerErr(:)];
f = [torqueErr(:)];

end

%% BUILDPARAMCONTROL
function [paramControl, T, Y] = buildParamControl(peakTorque, peakTime, riseTime, fallTime, auxdata)

startTime = auxdata.startTime;
startIdx = auxdata.startIdx;
timeShift = auxdata.timeShift;

numPts = length(timeShift(startIdx:end));
numZeros = size(auxdata.T_exo_norm, 1) - numPts;

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

paramControl = [zeros(1,numZeros) interp1(t, v, tq)]';
T = [node1_t node2_t node3_t]*(timeShift(end)-startTime) + startTime;
Y = [node1_y node2_y node3_y];

end

%% CALCDEVICETORQUE
function [torque] = calcDeviceTorque(x, auxdata)

% create empty torque array
torque = zeros(size(auxdata.T_exo_norm));

x = unscaleVariables(x, auxdata);

% extract curve variables
peakTorque = x(1);
peakTime = x(2);
riseTime = x(3);
fallTime = x(4);

mod_name = auxdata.mod_name;

if contains(mod_name,'actHf_') || contains(auxdata.mod_name,'actHe_')
    torque(:,1) = buildParamControl(peakTorque, peakTime, riseTime, fallTime, auxdata);
    
elseif contains(mod_name,'actKf_') || contains(auxdata.mod_name,'actKe_')
    torque(:,2) = buildParamControl(peakTorque, peakTime, riseTime, fallTime, auxdata);

elseif contains(mod_name, 'actAp_') || contains(auxdata.mod_name,'actAd_')
    torque(:,3) = buildParamControl(peakTorque, peakTime, riseTime, fallTime, auxdata);

elseif strcmp(mod_name,'mrsmod_actHfAp') 
    anklePeakTorque = x(5);
    torque(:,1) = buildParamControl(peakTorque, peakTime, riseTime, fallTime, auxdata);     
    torque(:,3) = buildParamControl(anklePeakTorque, peakTime, riseTime, fallTime, auxdata);
    
elseif strcmp(mod_name,'mrsmod_actHfKf') 
    kneePeakTorque = x(5);
    torque(:,1) = buildParamControl(peakTorque, peakTime, riseTime, fallTime, auxdata);     
    torque(:,2) = buildParamControl(kneePeakTorque, peakTime, riseTime, fallTime, auxdata);
    
elseif strcmp(mod_name,'mrsmod_actKfAp') 
    anklePeakTorque = x(5);
    torque(:,2) = buildParamControl(peakTorque, peakTime, riseTime, fallTime, auxdata);     
    torque(:,3) = buildParamControl(anklePeakTorque, peakTime, riseTime, fallTime, auxdata);

elseif strcmp(mod_name,'mrsmod_actHfKfAp')
    kneePeakTorque = x(5);
    anklePeakTorque = x(6);
    torque(:,1) = buildParamControl(peakTorque, peakTime, riseTime, fallTime, auxdata);
    torque(:,2) = buildParamControl(kneePeakTorque, peakTime, riseTime, fallTime, auxdata);
    torque(:,3) = buildParamControl(anklePeakTorque, peakTime, riseTime, fallTime, auxdata);    
    
elseif strcmp(mod_name,'mrsmod_actHeKe')
    kneePeakTorque = x(5);
    torque(:,1) = buildParamControl(peakTorque, peakTime, riseTime, fallTime, auxdata);
    torque(:,2) = buildParamControl(kneePeakTorque, peakTime, riseTime, fallTime, auxdata);
    
elseif strcmp(mod_name, 'mrsmod_actHfKfAp_multControls')
    peakTorque2 = x(5);
    peakTime2 = x(6);
    riseTime2 = x(7);
    fallTime2 = x(8);
    peakTorque3 = x(9);
    peakTime3 = x(10);
    riseTime3 = x(11);
    fallTime3 = x(12);
    torque(:,1) = buildParamControl(peakTorque, peakTime, riseTime, fallTime, auxdata);
    torque(:,2) = buildParamControl(peakTorque2, peakTime2, riseTime2, fallTime2, auxdata);
    torque(:,3) = buildParamControl(peakTorque3, peakTime3, riseTime3, fallTime3, auxdata);        
    
elseif strcmp(mod_name, 'mrsmod_actHfKf_multControls')
    peakTorque2 = x(5);
    peakTime2 = x(6);
    riseTime2 = x(7);
    fallTime2 = x(8);
    torque(:,1) = buildParamControl(peakTorque, peakTime, riseTime, fallTime, auxdata);
    torque(:,2) = buildParamControl(peakTorque2, peakTime2, riseTime2, fallTime2, auxdata);
    
elseif strcmp(mod_name, 'mrsmod_actHfAp_multControls')
    peakTorque2 = x(5);
    peakTime2 = x(6);
    riseTime2 = x(7);
    fallTime2 = x(8);
    torque(:,1) = buildParamControl(peakTorque, peakTime, riseTime, fallTime, auxdata);
    torque(:,3) = buildParamControl(peakTorque2, peakTime2, riseTime2, fallTime2, auxdata);
    
elseif strcmp(mod_name, 'mrsmod_actKfAp_multControls')
    peakTorque2 = x(5);
    peakTime2 = x(6);
    riseTime2 = x(7);
    fallTime2 = x(8);
    torque(:,2) = buildParamControl(peakTorque, peakTime, riseTime, fallTime, auxdata);
    torque(:,3) = buildParamControl(peakTorque2, peakTime2, riseTime2, fallTime2, auxdata);
    
elseif strcmp(mod_name, 'mrsmod_actHeKe_multControls')
    peakTorque2 = x(5);
    peakTime2 = x(6);
    riseTime2 = x(7);
    fallTime2 = x(8);
    torque(:,1) = buildParamControl(peakTorque, peakTime, riseTime, fallTime, auxdata);
    torque(:,2) = buildParamControl(peakTorque2, peakTime2, riseTime2, fallTime2, auxdata);
end

torque = torque.*auxdata.mapControl2DOFs;

end

%% CALCDEVICEPOWER
function [P, P_pos, P_neg, P_pos_avg, P_neg_avg] = calcDevicePower(torque, dqdt)

torque(end,:) = [];
P = torque.*dqdt;
P_pos = P;
P_neg = P;
P_pos(P_pos < 0) = 0;
P_neg(P_neg > 0) = 0;
P_pos_avg = sum(mean(P_pos));
P_neg_avg = sum(mean(P_neg));

end


