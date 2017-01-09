numDOFs = DatStore.nDOF;
numMuscles = DatStore.nMuscles;

time = OptInfo.result.solution.phase.time;

% Extract experimental data.
expTime = DatStore.time;
qExp = DatStore.q_exp;
momArmsExp = DatStore.dM;
momArms = interp1(expTime, momArmsExp, time);
jointAngles = pi / 180. * interp1(expTime, qExp, time);

% Extract parts of the solution related to the device.
control = OptInfo.result.solution.phase.control;
state = OptInfo.result.solution.phase.state;

% Joint moment breakdown.
deviceIndices = strmatch('ankle_angle', DatStore.DOFNames);
assert(length(deviceIndices) == 1);

figure;
for idof = 1:numDOFs
    subplot(numDOFs, 1, idof);
    hold on;
    plot(expTime, DatStore.T_exp(:, idof), 'k', 'LineWidth', 2);
    legendEntries = {'net'};
    sumMoment = zeros(length(TForce(:, 1)), 1);
    for imusc = 1:numMuscles
        if any(momArms(:, idof, imusc)) > 0.00001
            thisMoment = TForce(:, imusc) .* momArms(:, idof, imusc);
            plot(time(1:end-1), thisMoment(1:end-1));
            legendEntries = [legendEntries {MuscleNames{imusc}}];
            sumMoment = sumMoment + thisMoment;
        end
    end
     deviceIndex = find(deviceIndices == idof);
     if ~isempty(deviceIndex)
         normSpringStiff = OptInfo.result.solution.parameter(1);
         maxSpringStiff = 400; % N-m/rad.
         rest_angle = OptInfo.result.solution.parameter(2);
         ankleAngle = -(jointAngles(:, idof) - rest_angle);
         deviceMoment = maxSpringStiff * normSpringStiff .* ankleAngle;
         plot(time, deviceMoment);
         legendEntries = [legendEntries {'device'}];
         sumMoment = sumMoment + deviceMoment;
     end
    plot(time(1:end-1), sumMoment(1:end-1), 'r', 'LineWidth', 2);
    legendEntries = [legendEntries {'sum'}];
    legend(legendEntries, 'Interpreter', 'none');
    title(DatStore.DOFNames{idof}, 'Interpreter', 'none');
    if idof == numDOFs
        xlabel('time (s)');
    end
    ylabel('moment (N-m)');
end

% TODO percent reduction in metabolic cost / sum squared activation and
% squared excitation.
% TODO left and right limbs together.





















