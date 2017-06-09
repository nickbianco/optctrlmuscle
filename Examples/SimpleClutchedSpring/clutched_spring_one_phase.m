function output = clutched_spring_one_phase()
% DOES NOT WORK (YET).
% (Read the description of clutched_spring_with_phases.m first) In this
% formulation, we use one phase, but introduce two parameters that describe
% when the clutch is engaged.

auxdata.m = 10;
auxdata.g = 10;

%-------------------------------------------------------------------------%
%----------------- Provide All Bounds for Problem ------------------------%
%-------------------------------------------------------------------------%
t0 = 0;
tf = 5;

% First state, x: the depth of the mass (x increases downwards).
% Second state, u: rate of change of depth of the mass.
initialstate = [0 0];
finalstate = [0 0];

state_lower = [-5 -100];
state_upper = [2000 100];

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

% Spring stiffness, time of clutch engagement, duration of clutch engagement.
bounds.parameter.lower = [0, t0, 0];
bounds.parameter.upper = [100, tf, (tf - t0)];

% First event group is for the phase1-phase2 boundary.
% The 3 constraints are that the phases are sequential in time and that the
% final states of phase 1 match the final states of phase 2.
% The disengagement time cannot be past tf.
bounds.eventgroup.lower = [-inf];
bounds.eventgroup.upper = [0];

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
guess.phase.time    = linspace(t0, tf, 4)';
guess.phase.state   = [0, 0;
                       1, 1;
                       1, -1;
                       0, 0];
guess.parameter = [8, t0 + 0.3 * (tf - t0), 0.3 * (tf - t0)];
%guess.parameter = [8, t0, (tf - t0)];

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
gpopsVerifySetup(setup);
output = gpops2(setup);
figure;
hold on;
sol = output.result.solution;
plot(sol.phase.time, -sol.phase.state(:, 1));
fprintf('Spring stiffness: %d\n', sol.parameter(1));
fprintf('Clutch engagement duration: %d\n', sol.parameter(3));
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
g = input.auxdata.g;
x = input.phase.state(:, 1);
u = input.phase.state(:, 2);
%% GPOPS provides the same parameter for each time point (size(k) = N x 1), 
%% for numerical efficiency reasons.
k = input.phase.parameter(:, 1);
engagement_time = input.phase.parameter(1, 2); % TODO not using these as vectors.
engagement_duration = input.phase.parameter(1, 3);
disengagement_time = engagement_time + engagement_duration;
time = input.phase.time;
N = length(time);
udot = g * ones(N, 1);
% TODO + 1 / m * -engagement_start(1)

% The rest length for the spring is the length of the 'path' when the clutch
% engages.
engagement_length = x(find(time > engagement_time, 1, 'first'));
engagement_indices = find(time > engagement_time & time < disengagement_time);
spring_accel = 1 / m * -k(engagement_indices) .* (x(engagement_indices) - engagement_length);
udot(engagement_indices) = udot(engagement_indices) + spring_accel;

phaseout.dynamics = [u, udot];

% TODO not needed disengagement_length = x(find(time == disengagement_time, 1, 'first'));
%for i = 1:N
%    if input.phase.time(i) > engagement_start
%
%end
%udot = g + (is_engaged) .* 1/m* -k .* (x - rest_length);

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

k = input.parameter(1);
engagement_duration = input.parameter(3);
output.objective = engagement_duration + 0.1 * k;

% disengagement_time <= finaltime.
engagement_time = input.parameter(2);
disengagement_time = engagement_time + engagement_duration;
output.eventgroup.event = disengagement_time - input.phase.finaltime;

end
