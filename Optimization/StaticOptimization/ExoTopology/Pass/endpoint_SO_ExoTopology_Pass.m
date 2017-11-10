function output = endpoint_SO_ExoTopology_Pass(input)

q = input.phase.integral;
output.objective = q;

% Constraints - mild periodicity
a_end = input.phase.finalstate;
a_init = input.phase.initialstate;
pera = a_end - a_init;
output.eventgroup.event = pera;