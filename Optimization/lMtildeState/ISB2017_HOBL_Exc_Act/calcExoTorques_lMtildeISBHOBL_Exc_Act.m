function ExoTorques = calcExoTorques_lMtildeISBHOBL_Exc_Act(OptInfo,DatStore)

time = OptInfo.result.solution.phase.time;
normSpringStiff = OptInfo.result.solution.parameter(1);
deviceRestLength = OptInfo.result.solution.parameter(2);

ExoTorques = zeros(length(time), length(DatStore.DOFNames));
exoLength = ppval(OptInfo.result.setup.auxdata.JointLenEXOSpline,time);
maxSpringStiff = 500; % N/m

for dof = 1:length(DatStore.DOFNames)
    exoMomentArms = ppval(OptInfo.result.setup.auxdata.JointMAEXOSpline(dof),time);
    deltaX = exoLength-deviceRestLength;
    deltaX = max(deltaX,0);
    ExoTorques(:,dof) = maxSpringStiff.*normSpringStiff.*exoMomentArms.*deltaX;    
end