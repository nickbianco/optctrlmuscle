function sstruct = SplineInputData_MultiPhase(t,input,p)

numColPoints = length(t);
NMuscles = input.auxdata.NMuscles;
Ndof = input.auxdata.Ndof;

sstruct.LMT = zeros(numColPoints,NMuscles);
sstruct.VMT = zeros(numColPoints,NMuscles);
sstruct.EMG = zeros(numColPoints,NMuscles);
sstruct.FL = zeros(numColPoints,NMuscles);
sstruct.FV = zeros(numColPoints,NMuscles);

for dof = 1:Ndof
    for m = 1:NMuscles
        index_sel=(dof-1)*(NMuscles)+m;
        sstruct.MA(:,index_sel) = ppval(input.auxdata.JointMASpline(dof).Muscle(p).p(m),t);   
    end
    sstruct.ID(:,dof) = ppval(input.auxdata.JointIDSpline(p).p(dof),t);
    sstruct.EXO(:,dof) = ppval(input.auxdata.JointEXOSpline(p).p(dof),t);
    sstruct.IK(:,dof) = pi / 180. * ppval(input.auxdata.JointIKSpline(p).p(dof),t);
end

for m = 1:NMuscles
    [sstruct.LMT(:,m),sstruct.VMT(:,m),~] = SplineEval_ppuval(input.auxdata.LMTSpline(p).p(m),t,1);
    if isfield(input.auxdata, 'EMGSpline')
        sstruct.EMG(:,m) = ppval(input.auxdata.EMGSpline(p).p(m),t);
    end
    if isfield(input.auxdata, 'FLSpline')
        sstruct.FL(:,m) = ppval(input.auxdata.FLSpline(p).p(m),t);
    end
    if isfield(input.auxdata, 'FVSpline')
        sstruct.FV(:,m) = ppval(input.auxdata.FVSpline(p).p(m),t);
    end
end
