function DatStore = SolveStaticOptimization_ExoTopology(DatStore, auxdata, Misc)

% Problem bounds
e_min = 0; e_max = 1;             % bounds on muscle excitation
a_min = 0; a_max = 1;

% Time bounds
t0 = DatStore.time(1); 
tf = DatStore.time(end);
bounds.phase.initialtime.lower = t0; 
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tf; 
bounds.phase.finaltime.upper = tf;
auxdata.initialtime = t0;
auxdata.finaltime = tf;

% Controls bounds
excMin = e_min*ones(1,auxdata.NMuscles); 
excMax = e_max*ones(1,auxdata.NMuscles);
aTmin = -1*ones(1,auxdata.Ndof); 
aTmax = 1*ones(1,auxdata.Ndof);
aDmin = zeros(1, auxdata.numActiveDOFs); 
aDmax = ones(1, auxdata.numActiveDOFs);
if auxdata.hasActiveDevice
    control_bounds_lower = [excMin aTmin aDmin];
    control_bounds_upper = [excMax aTmax aDmax];
else
    control_bounds_lower = [excMin aTmin];
    control_bounds_upper = [excMax aTmax];
end
bounds.phase.control.lower = control_bounds_lower; 
bounds.phase.control.upper = control_bounds_upper;

% States bounds
actMin = a_min*ones(1,auxdata.NMuscles); 
actMax = a_max*ones(1,auxdata.NMuscles);
bounds.phase.initialstate.lower = actMin; 
bounds.phase.initialstate.upper =  actMax;
bounds.phase.state.lower = actMin; 
bounds.phase.state.upper =  actMax;
bounds.phase.finalstate.lower = actMin; 
bounds.phase.finalstate.upper =  actMax;

% Integral bounds
bounds.phase.integral.lower = 0;
bounds.phase.integral.upper = 10000*(tf-t0);

% Parameter bounds
% Parameter variable in the optimization problem always lie in the range 
% [-1 1], and should be converted as necessary in the continuous function
% for to make the correct computations.
bounds.parameter.lower = -1*ones(size(auxdata.paramsLower));
bounds.parameter.upper = 1*ones(size(auxdata.paramsUpper));

% Path constraints
ID_bounds = zeros(1, auxdata.Ndof);
bounds.phase.path.lower = ID_bounds;
bounds.phase.path.upper = ID_bounds;

% Eventgroup
% Impose mild periodicity
pera_lower = -1*ones(1, auxdata.NMuscles);
pera_upper = 1*ones(1, auxdata.NMuscles);
bounds.eventgroup.lower = pera_lower; 
bounds.eventgroup.upper = pera_upper;

% Initial guess
N = length(DatStore.time);
guess.phase.time = DatStore.time;
if auxdata.hasActiveDevice 
    control_guess = [0.2*ones(N,auxdata.NMuscles) zeros(N,auxdata.Ndof) 0.2*ones(N,auxdata.numActiveDOFs)];
else
    control_guess = [0.2*ones(N,auxdata.NMuscles) zeros(N,auxdata.Ndof)];
end

guess.phase.control = control_guess;

guess.phase.state = 0.2*ones(N,auxdata.NMuscles);
guess.phase.integral = 0;
if auxdata.hasPassiveDevice
    guess.parameter = [zeros(1,auxdata.numExoParams-1) 1];
else
    guess.parameter = zeros(1,auxdata.numExoParams);
end

% Spline structures
for dof = 1:auxdata.Ndof
    for m = 1:auxdata.NMuscles       
        auxdata.JointMASpline(dof).Muscle(m) = spline(DatStore.time,auxdata.MA(dof).Joint(:,m));       
    end
    auxdata.JointIDSpline(dof) = spline(DatStore.time,DatStore.T_exp(:,dof));
    auxdata.JointIKSpline(dof) = spline(DatStore.time,DatStore.q_exp(:,dof));
end

for m = 1:auxdata.NMuscles
    auxdata.LMTSpline(m) = spline(DatStore.time,DatStore.LMT(:,m));
end

% GPOPS setup
setup.name = 'StaticOptimization_ExoTopology';
setup.auxdata = auxdata;
setup.bounds = bounds;
setup.guess = guess;
setup.nlp.solver = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.nlp.ipoptoptions.tolerance = 1e-3;
setup.nlp.ipoptoptions.maxiterations = 10000;
setup.derivatives.supplier = 'sparseCD';
setup.derivatives.derivativelevel = 'first';
setup.derivatives.dependencies = 'sparse';
setup.scales.method = 'none';
setup.mesh.method = 'hp-PattersonRao';
setup.mesh.tolerance = 1e-3;
setup.mesh.maxiterations = 20;
setup.mesh.colpointsmin = 5;
setup.mesh.colpointsmax = 10;
setup.method = 'RPM-integration';
setup.displaylevel = 2;
NMeshIntervals = round((tf-t0)*Misc.Mesh_Frequency);
setup.mesh.phase.colpoints = 5*ones(1,NMeshIntervals);
setup.mesh.phase.fraction = (1/(NMeshIntervals))*ones(1,NMeshIntervals);
setup.functions.continuous = str2func(['continous_SO_ExoTopology_' auxdata.subcase]);
setup.functions.endpoint = str2func(['endpoint_SO_ExoTopology_' auxdata.subcase]);

input.auxdata = auxdata;
tdummy = guess.phase.time;
splinestruct = SplineInputData_SO(tdummy,input);
splinenames = fieldnames(splinestruct);
for Scount = 1:length(splinenames)
    secdim = size(splinestruct.(splinenames{Scount}),2);
    splinestructad.(splinenames{Scount}) = adigatorCreateAuxInput([Inf,secdim]);
    splinestruct.(splinenames{Scount}) = zeros(0,secdim);
end
setup.auxdata.splinestruct = splinestructad;
% adigatorGenFiles4gpops2(setup)
setup.functions.continuous = str2func(['Wrap4continous_SO_ExoTopology_' auxdata.subcase]);

output = gpops2(setup);
result = output.result.solution.phase(1);
DatStore.SO_MAct = interp1(result.time, ...
    result.control(:,1:auxdata.NMuscles), DatStore.time);
DatStore.SO_RAct = interp1(result.time,  ... 
    result.control(:,auxdata.NMuscles+1:auxdata.NMuscles+auxdata.Ndof), DatStore.time);
DatStore.SO_parameter = output.result.solution.parameter;

if auxdata.hasActiveDevice && ~strcmp(auxdata.subcase, 'ActParam')
    DatStore.SO_ExoAct = interp1(result.time, result.control(:,end-(auxdata.numActiveDOFs-1):end), DatStore.time);
end

end