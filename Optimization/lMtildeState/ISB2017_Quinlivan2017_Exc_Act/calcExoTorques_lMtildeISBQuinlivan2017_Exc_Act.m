function ExoTorques = calcExoTorques_lMtildeISBQuinlivan2017_Exc_Act(OptInfo,DatStore)

time = OptInfo.result.solution.phase.time;
auxdata = OptInfo.result.setup.auxdata;
p = auxdata.p_linreg;
exo_force_level = OptInfo.result.solution.parameter;
ExoTorques = zeros(length(time), length(DatStore.DOFNames));
for dof = 1:length(DatStore.DOFNames)
    exoPeak = p(1,dof)*exo_force_level + p(2,dof);
    exoNormTorque = ...
        ppval(OptInfo.result.setup.auxdata.JointEXOSpline(dof),time);
    ExoTorques(:,dof) = exoPeak*exoNormTorque;
end

end
