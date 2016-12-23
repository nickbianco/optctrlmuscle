function output = clutched_spring_with_phases()
% A mass undergoes a motion, in the presence of gravity, in which the mass
% starts and ends at the same height and at rest. We prescribe the duration of
% this motion. The mass hangs from a clutched spring, and this is how it is
% able to raise back up to its initial height. We seek the spring stiffness and
% the times at which the clutch engages and disengages in order to arrive at
% the final state after the specified duration. There are infinitely many
% combinations of spring stiffness and clutch engagement duration that can
% cause this motion. To obtain a unique solution, we choose an objective
% function that penalizes both the clutch engagement duration and the spring
% stiffness; the relative weight between these two terms affects the unique
% solution.
% We solve this problem in 3 phases:
%   1. free fall
%   2. clutch engaged
%   3. free "fall"

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

%F_lower = -50;
%F_upper = 50;

%-------------------------------------------------------------------------%
%----------------------- Setup for Problem Bounds ------------------------%
%-------------------------------------------------------------------------%
bounds.phase(1).initialtime.lower = t0;
bounds.phase(1).initialtime.upper = t0;
% The final time of the 1st phase through the initial time of the 3rd
% phase can be anywhere between t0 and tf; we do not constrain when the clutch
% is engaged.
bounds.phase(1).finaltime.lower = t0;
bounds.phase(1).finaltime.upper = tf;
bounds.phase(1).initialstate.lower = initialstate;
bounds.phase(1).initialstate.upper = initialstate;
bounds.phase(1).state.lower = state_lower;
bounds.phase(1).state.upper = state_upper;
% The final state of the 1st phase can be anything.
bounds.phase(1).finalstate.lower = state_lower;
bounds.phase(1).finalstate.upper = state_upper;

bounds.phase(2).initialtime.lower = t0 + 0.1 * (tf - t0);
bounds.phase(2).initialtime.upper = tf;
bounds.phase(2).finaltime.lower = t0;
bounds.phase(2).finaltime.upper = tf;
% The initial state of the 2nd phase can be anything.
bounds.phase(2).initialstate.lower = state_lower;
bounds.phase(2).initialstate.upper = state_upper;
bounds.phase(2).state.lower = state_lower;
bounds.phase(2).state.upper = state_upper;
% The final state of the 2nd phase can be anything.
bounds.phase(2).finalstate.lower = state_lower;
bounds.phase(2).finalstate.upper = state_upper;

bounds.phase(3).initialtime.lower = t0;
bounds.phase(3).initialtime.upper = tf;
bounds.phase(3).finaltime.lower = tf;
bounds.phase(3).finaltime.upper = tf;
% The initial state of the 3nd phase can be anything.
bounds.phase(3).initialstate.lower = state_lower;
bounds.phase(3).initialstate.upper = state_upper;
bounds.phase(3).state.lower = state_lower;
bounds.phase(3).state.upper = state_upper;
% The final state of the 3rd phase is the actual final state.
bounds.phase(3).finalstate.lower = finalstate;
bounds.phase(3).finalstate.upper = finalstate;

% The only parameter is the spring stiffness.
bounds.parameter.lower = [0];
bounds.parameter.upper = [100];

% First event group is for the phase1-phase2 boundary.
% The 3 constraints are that the phases are sequential in time and that the
% final states of phase 1 match the final states of phase 2.
bounds.eventgroup(1).lower = [0, 0, 0];
bounds.eventgroup(1).upper = [0, 0, 0];
% Similarly for the phase2-phase3 boundary.
bounds.eventgroup(2).lower = [0, 0, 0];
bounds.eventgroup(2).upper = [0, 0, 0];

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
guess.phase(1).time    = [t0; t0 + 1/3 * (tf - t0)];
guess.phase(1).state   = [0, 0; 1, 1];
guess.phase(2).time    = [t0 + 1/3 * (tf - t0); t0 + 2/3 * (tf - t0)];
guess.phase(2).state   = [1, 1; 1, -1];
guess.phase(3).time    = [t0 + 2/3 * (tf - t0); tf];
guess.phase(3).state   = [1, -1; 0, 0];
guess.parameter = [8];

%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method          = 'hp-PattersonRao';
mesh.tolerance       = 1e-5;
mesh.maxiterations   = 10; % TODO 
mesh.colpointsmin    = 3;
mesh.colpointsmax    = 10;
N                    = 10;
mesh.phase(1).colpoints = 3*ones(1,N);
mesh.phase(1).fraction  = ones(1,N)/N;
mesh.phase(2).colpoints = 3*ones(1,N);
mesh.phase(2).fraction  = ones(1,N)/N;
mesh.phase(3).colpoints = 3*ones(1,N);
mesh.phase(3).fraction  = ones(1,N)/N;

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
for i = 1:3
    plot(sol.phase(i).time, -sol.phase(i).state(:, 1));
end
fprintf('Spring stiffness: %d\n', sol.parameter);
fprintf('Clutch engagement duration: %d\n', ...
    sol.phase(2).time(end) - sol.phase(2).time(1));
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
x1 = input.phase(1).state(:, 1);
u1 = input.phase(1).state(:, 2);
x2 = input.phase(2).state(:, 1);
u2 = input.phase(2).state(:, 2);
x3 = input.phase(3).state(:, 1);
u3 = input.phase(3).state(:, 2);
%% GPOPS provides the same parameter for each time point (size(k) = N x 1), 
%% for numerical efficiency reasons.
k = input.phase(2).parameter;

phaseout(1).dynamics = [u1, g * ones(length(input.phase(1).time), 1)];
% The rest length for the spring is the length of the 'path' when the clutch
% engages.
Lrest = x2(1);
phaseout(2).dynamics = [u2, 1 / m * -k .* (x2 - Lrest) + g];
phaseout(3).dynamics = [u3, g * ones(length(input.phase(3).time), 1)];

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

clutch_engagement_duration = ...
        input.phase(2).finaltime - input.phase(2).initialtime;
output.objective = clutch_engagement_duration + 0.1 * input.parameter;

output.eventgroup(1).event = [...
    input.phase(1).finaltime - input.phase(2).initialtime, ...
    input.phase(1).finalstate - input.phase(2).initialstate];

output.eventgroup(2).event = [...
    input.phase(2).finaltime - input.phase(3).initialtime, ...
    input.phase(2).finalstate - input.phase(3).initialstate];

end
