function ExoTorques = calcExoTorques_FtildeAHE(OptInfo,DatStore)

time = OptInfo.result.solution.phase.time;
auxdata = OptInfo.result.setup.auxdata;
ExoTorques = zeros(length(time), length(DatStore.DOFNames));
for dof = 1:length(DatStore.DOFNames)
    exoTorque = ...
        ppval(OptInfo.result.setup.auxdata.JointEXOSpline(dof),time);
    ExoTorques(:,dof) = exoTorque;
end

end
