function output = musdynEndpoint_lMtildeState_Exc_ActPh(input)

output.objective = 0;
for ip = 1:length(input.phase)
    output.objective = output.objective + input.phase(ip).integral;
end

NMuscles = input.auxdata.NMuscles;

% Initial and end states
a_end = input.phase(end).finalstate(1:NMuscles);
lMtilde_end = input.phase(end).finalstate(NMuscles+1:end);

a_init = input.phase(1).initialstate(1:NMuscles);
lMtilde_init = input.phase(1).initialstate(NMuscles+1:end);

% Constraints - mild periodicity
pera = a_end - a_init;
perlMtilde = lMtilde_end - lMtilde_init;

output.eventgroup(1).event = [pera perlMtilde];
output.eventgroup(2).event = [...
    input.phase(1).finalstate - input.phase(2).initialstate];
