function ExoTorques = calcExoTorques_lMtildeISBCollins2015_Exc_Act(OptInfo,DatStore)

time = OptInfo.result.solution.phase.time;
auxdata = OptInfo.result.setup.auxdata;

normSpringStiff = OptInfo.result.solution.parameter(1);
springRestAngle = OptInfo.result.solution.parameter(2);

first_peak = auxdata.rest_length_first_peak;
after_recoil = auxdata.rest_length_after_recoil;
beginSpringStretching = 1 ./ (1 + exp(100 * (first_peak  - time)));
restLengthReached = 1 ./ (1 + exp(100 * (time - after_recoil)));

isSpringActive = beginSpringStretching .* restLengthReached;

ExoTorques = zeros(length(time), length(DatStore.DOFNames));

maxSpringStiff = 400; % N-m/rad.
for dof = 1:length(DatStore.DOFNames)
    if any(dof == auxdata.clutched_spring_dofs)
        jointAngle = pi / 180. * ppval(OptInfo.result.setup.auxdata.JointIKSpline(dof),time);
        springStretch = -(jointAngle - springRestAngle);
        ExoTorques(:,dof) = maxSpringStiff * normSpringStiff .* springStretch .* isSpringActive;
    end
end
