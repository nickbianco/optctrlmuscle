function output = test_SO

numControls = 2;

% Time bounds
t0 = 0;
tf = 1;
bounds.phase.initialtime.lower = t0; 
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tf; 
bounds.phase.finaltime.upper = tf;
auxdata.initialtime = t0;
auxdata.finaltime = tf;

% Controls bounds
bounds.phase.control.lower = zeros(1,numControls); 
bounds.phase.control.upper = ones(1,numControls);

% States bounds
bounds.phase.state.lower = 0;
bounds.phase.state.upper = 0;
bounds.phase.initialstate.lower = 0;
bounds.phase.initialstate.upper = 0;
bounds.phase.finalstate.lower = 0;
bounds.phase.finalstate.upper = 0;

% Integral bounds
bounds.phase.integral.lower = 0;
bounds.phase.integral.upper = 10000*(tf-t0);

% Path constraints
bounds.phase.path.lower = 0;
bounds.phase.path.upper = 0;

% Eventgroup
% Impose mild periodicity
% pera_lower = -1;
% pera_upper = 1;
% bounds.eventgroup.lower = [pera_lower]; 
% bounds.eventgroup.upper = [pera_upper];

% Initial guess
N = 101;
guess.phase.time = linspace(0,1,N)';
guess.phase.control = 0.5*ones(N,numControls);
guess.phase.state = zeros(N,1);
guess.phase.integral = 0;

% GPOPS setup
setup.name = 'test_SO';
setup.auxdata = auxdata;
setup.bounds = bounds;
setup.guess = guess;
% setup.nlp.solver = 'ipopt';
% setup.nlp.ipoptoptions.linear_solver = 'ma57';
% setup.derivatives.derivativelevel = 'second';
% setup.nlp.ipoptoptions.tolerance = 1e-6;
% setup.nlp.ipoptoptions.maxiterations = 10000;
% setup.derivatives.supplier = 'sparseCD';
% setup.scales.method = 'none';
% setup.mesh.method = 'hp-PattersonRao';
% setup.mesh.tolerance = 1e-3;
% setup.mesh.maxiterations = 20;
% setup.mesh.colpointsmin = 5;
% setup.mesh.colpointsmax = 10;
% setup.method = 'RPM-integration';
% setup.displaylevel = 2;
% NMeshIntervals = round((tf-t0)*Misc.Mesh_Frequency);
% setup.mesh.phase.colpoints = 5*ones(1,NMeshIntervals);
% setup.mesh.phase.fraction = (1/(NMeshIntervals))*ones(1,NMeshIntervals);
setup.functions.continuous = @continuous;
setup.functions.endpoint = @endpoint;

output = gpops2(setup);

end

function [phaseout] = continuous(input)

time = input.phase.time;
e  = input.phase.control;

Tprescribed = sin(2*pi*time);
Tmatch = 10*e(:,1) - e(:,2).^2; 
diff = Tprescribed - Tmatch;

phaseout.path = diff;

phaseout.dynamics = zeros(size(time));

phaseout.integrand = sum(e.^2,2);

end

function [output] = endpoint(input)

q = input.phase.integral;
output.objective = q;

end