function output = spring()
% Solve for the spring stiffness that moves a point mass from x = -1 to x = 1
% in the given amount of time. This is not an optimal control problem; there is
% only one solution (it is a boundary value problem). As such, there is no
% objective function.

auxdata.m = 10;

% The correct answer for the spring stiffness parameter.
k_answer = 5;

%-------------------------------------------------------------------------%
%----------------- Provide All Bounds for Problem ------------------------%
%-------------------------------------------------------------------------%
t0 = 0;
% Period of oscillation.
period = 2 * pi * sqrt(auxdata.m / k_answer);
% We simulate half a period (x = -1 to 1).
tf = 0.5 * period;

initialstate = [-1 0];
finalstate = [1 0];

state_lower = [-2 -10];
state_upper = [2 10];

%-------------------------------------------------------------------------%
%----------------------- Setup for Problem Bounds ------------------------%
%-------------------------------------------------------------------------%
bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tf;
bounds.phase.finaltime.upper = tf;
bounds.phase.initialstate.lower = initialstate;
bounds.phase.initialstate.upper = initialstate;
bounds.phase.state.lower = state_lower;
bounds.phase.state.upper = state_upper;
bounds.phase.finalstate.lower = finalstate;
bounds.phase.finalstate.upper = finalstate;
bounds.parameter.lower = [0];
bounds.parameter.upper = [10];

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
t_guess             = [t0; tf];

state_guess = [state_lower; state_upper];

guess.phase.state   = state_guess;
guess.phase.time    = [t_guess];
guess.parameter = [8];

%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method          = 'hp-PattersonRao';
mesh.tolerance       = 1e-5;
mesh.maxiterations   = 10;
mesh.colpointsmin    = 3;
mesh.colpointsmax    = 10;
N                    = 10;
mesh.phase.colpoints = 3*ones(1,N);
mesh.phase.fraction  = ones(1,N)/N;

%-------------------------------------------------------------------------%
%------------- Assemble Information into Problem Structure ---------------%
%-------------------------------------------------------------------------%
setup.mesh                            = mesh;
setup.name                            = 'spring';
setup.functions.endpoint              = @endpoint;
setup.functions.continuous            = @continuous;
setup.displaylevel                    = 2;
setup.auxdata                         = auxdata;
setup.bounds                          = bounds;
setup.guess                           = guess;
setup.nlp.solver                      = 'ipopt';
setup.nlp.ipoptoptions.linear_solver  = 'ma57';
setup.derivatives.supplier            = 'sparseCD';
setup.derivatives.derivativelevel     = 'second';
setup.method                          = 'RPM-Integration';
setup.nlp.ipoptoptions.tolerance      = 1e-7;
setup.derivatives.derivativelevel     = 'second';
setup.scales.method                   = 'automatic-guess';

%-------------------------------------------------------------------------%
%----------------------- Solve Problem Using GPOPS2 ----------------------%
%-------------------------------------------------------------------------%
output = gpops2(setup);
end

function phaseout = continuous(input)

% input
% input.phase(phasenumber).state
% input.phase(phasenumber).control
% input.phase(phasenumber).time
% input.phase(phasenumber).parameter
%
% input.auxdata = auxiliary information
%
% output
% phaseout(phasenumber).dynamics
% phaseout(phasenumber).path
% phaseout(phasenumber).integrand

m = input.auxdata.m;
x = input.phase.state(:, 1);
u = input.phase.state(:, 2);
%% GPOPS provides the same parameter for each time point (size(k) = N x 1), 
%% for numerical efficiency reasons.
k = input.phase.parameter;

phaseout.dynamics = [u, -k .* x / m];

end

function output = endpoint(input)

% Inputs
% input.phase(phasenumber).initialstate -- row
% input.phase(phasenumber).finalstate -- row
% input.phase(phasenumber).initialtime -- scalar
% input.phase(phasenumber).finaltime -- scalar
% input.phase(phasenumber).integral -- row
%
% input.parameter -- row

% input.auxdata = auxiliary information

% Output
% output.objective -- scalar
% output.eventgroup(eventnumber).event -- row

%output.objective = 0;
% GPOPS seems to require that I give an actual objective function.
output.objective = input.phase.finaltime;

end
