data = load('DatStore.mat');

% Find approx joint angle derivatives
q_exp_temp = data.DatStore.q_exp;
dq_exp = diff(q_exp_temp);

% Lose time step with diff, adjust other variables
time = data.DatStore.time(1:end-1,:);
q_exp = data.DatStore.q_exp(1:end-1,:);
T_exp = data.DatStore.T_exp(1:end-1,:);

% Joint powers
P_exp = dq_exp .* T_exp;

figure(1)
plot(time,q_exp(:,1))
hold on
plot(time,T_exp(:,1))
plot(time,dq_exp(:,1))
plot(time,P_exp(:,1))

legend({'angle','moment','velocity','power'})

figure(2)
plot(time,dq_exp(:,1))

[max_dq,max_time] = max(-dq_exp(:,1));
time(max_time)
[max_P,max_time] = max(P_exp(:,1));
time(max_time)

DingExoCurves = load('DingExoCurves.mat');
figure(3)
plot(time,dq_exp(:,1))
hold on
plot(DingExoCurves.time,-DingExoCurves.eslp.F/100)
