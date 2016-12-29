numDOFs = DatStore.nDOF;
numMuscles = DatStore.nMuscles;

time = OptInfo.result.solution.phase.time;

% Extract experimental data.
expTime = DatStore.time;
qExp = DatStore.q_exp;
momArmsExp = DatStore.dM;
momArms = interp1(expTime, momArmsExp, time);

%for idof = 1:numDOFs
%    qSpline(idof) = spline(expTime, pi / 180. * qExp(:, idof));
%end

% Extract parts of the solution related to the device.
control = OptInfo.result.solution.phase.control;
state = OptInfo.result.solution.phase.state;

% Joint moment breakdown.
deviceIndices = strmatch('ankle_angle', DatStore.DOFNames);

Topt = 2*150.0;
for idof = 1:numDOFs
    %subplot(numDOFs, 1, idof);
    figure;
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
%     deviceIndex = find(deviceIndices == idof);
%     if ~isempty(deviceIndex)
%         deviceMoment = TODO;
%         plot(time, deviceMoment);
%         legendEntries = [legendEntries {'device'}];
%         sumMoment = sumMoment + deviceMoment;
%     end
    plot(time(1:end-1), sumMoment(1:end-1), 'r', 'LineWidth', 2);
    legendEntries = [legendEntries {'sum'}];
    legend(legendEntries, 'Interpreter', 'none');
    title(DatStore.DOFNames{idof}, 'Interpreter', 'none');
end

% TODO percent reduction in metabolic cost / sum squared activation and
% squared excitation.
% TODO left and right limbs together.





















