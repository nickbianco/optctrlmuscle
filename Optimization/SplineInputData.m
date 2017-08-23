function sstruct = SplineInputData(t,input)

numColPoints = length(t);
NMuscles = input.auxdata.NMuscles;
Ndof = input.auxdata.Ndof;

sstruct.LMT = zeros(numColPoints,NMuscles);
sstruct.VMT = zeros(numColPoints,NMuscles);
sstruct.LenEXO = ppval(input.auxdata.JointLenEXOSpline,t);

for dof = 1:Ndof
    for m = 1:NMuscles
        index_sel=(dof-1)*(NMuscles)+m;
        sstruct.MA(:,index_sel) = ppval(input.auxdata.JointMASpline(dof).Muscle(m),t);   
    end
    sstruct.ID(:,dof) = ppval(input.auxdata.JointIDSpline(dof),t);
    sstruct.EXO(:,dof) = ppval(input.auxdata.JointEXOSpline(dof),t);
    sstruct.MAEXO(:,dof) = ppval(input.auxdata.JointMAEXOSpline(dof),t);
    sstruct.IK(:,dof) = pi / 180. * ppval(input.auxdata.JointIKSpline(dof),t);
end

for m = 1:NMuscles
    [sstruct.LMT(:,m),sstruct.VMT(:,m),~] = SplineEval_ppuval(input.auxdata.LMTSpline(m),t,1);
end
