function output = clutched_spring_smooth_clutching()
% DOES NOT WORK (YET).
% (Read the description of clutched_spring_with_phases.m first) In this
% formulation, we use an additional control variable to keep track of the rest
% length (actually, the stretch) for the clutched spring, and a control
% continuous control variable to dictate when the clutch is engaged or not
% engaged.
% Note that we DO get the solution we expect if we set control_lower to 1 (so
% that the spring is just a normal spring; "always clutched").

auxdata.m = 10;
auxdata.g = 10;
auxdata.logistic_steepness = 30;
auxdata.stretch_decay_time_const = 0.010;

%-------------------------------------------------------------------------%
%----------------- Provide All Bounds for Problem ------------------------%
%-------------------------------------------------------------------------%
t0 = 0;
tf = 5;

% First state, x: the depth of the mass (x increases downwards).
% Second state, u: rate of change of depth of the mass.
initialstate = [0 0 0];
finalstate = [0 0 0];

state_lower = [-5 -100 -50];
state_upper = [2000 100 50];

% Is the clutch on or off?
control_lower = [0];
control_upper = [1];

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
bounds.phase.control.lower = control_lower;
bounds.phase.control.upper = control_upper;

% The only parameter is the spring stiffness.
bounds.parameter.lower = [0];
bounds.parameter.upper = [100];
bounds.phase.integral.lower = -1000;
bounds.phase.integral.upper = 0;

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
guess.phase.time    = linspace(t0, tf, 4)';
guess.phase.state   = [0, 0, 0;
                       1, 1, 0;
                       1, -1, 0;
                       0, 0, 0];
guess.phase.control = [0;
                       1;
                       1;
                       0];
guess.phase.integral = -1;
guess.parameter = [8];

%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method          = 'hp-PattersonRao';
mesh.tolerance       = 1e-5;
mesh.maxiterations   = 5; 
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
subplot(2, 1, 1);
plot(sol.phase.time, -sol.phase.state(:, 1));
subplot(2, 1, 2);
plot(sol.phase.time, sol.phase.control(:, 1));

fprintf('Spring stiffness: %d\n', sol.parameter(1));
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
logistic_steep = input.auxdata.logistic_steepness;
decay = input.auxdata.stretch_decay_time_const;
x = input.phase.state(:, 1);
u = input.phase.state(:, 2);
stretch = input.phase.state(:, 3);
clutch = input.phase.control(:, 1);
%% GPOPS provides the same parameter for each time point (size(k) = N x 1), 
%% for numerical efficiency reasons.
k = input.phase.parameter(:, 1);

udot = g + 1 / m * -k .* stretch;

%clutch = 
N = length(input.phase.time);
stretchdot = zeros(N, 1);
for i = 1:N
    if clutch(i) > 0.5
        stretchdot(i) = u(i);
    else
        stretchdot(i) = - 1 / decay * stretch(i);
    end
end

phaseout.dynamics = [u, udot, stretchdot];

% The smallest penalty is for clutch = 1 and clutch = 0.
phaseout.integrand = -(clutch - 0.5).^2;

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

% Arbitrary objective function (for now).
% Penalize clutch so that it is either 0 or 1.
output.objective = input.phase.integral; %input.phase.finaltime;

end
