function ExoTorques = calcExoTorques_lMtildeExoHipAnkle_Exc_Act(OptInfo,DatStore)

time = OptInfo.result.solution.phase.time;
auxdata = OptInfo.result.setup.auxdata;
Fopt_exo = auxdata.Fopt_exo;
tradeoff = auxdata.tradeoff;
alpha = OptInfo.result.solution.parameter;
aD = OptInfo.result.solution.phase.control(:,end);
r = 0.1;

ExoTorques = zeros(length(time), length(DatStore.DOFNames));
for dof = 1:length(DatStore.DOFNames)
    ExoTorques(:,dof) = r*Fopt_exo(dof)*aD.*(ones(length(time),1)+tradeoff(dof)*alpha);
end

end
