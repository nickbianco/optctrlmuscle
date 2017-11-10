function phaseout = Wrap4continous_SO_ExoTopology_Pass(input)

persistent splinestruct

if isempty(splinestruct) || size(splinestruct.MA,1) ~= length(input.phase.time) 
    splinestruct = SplineInputData_SO(input.phase.time,input);
end

input.auxdata.splinestruct = splinestruct;

phaseout = continous_SO_ExoTopology_Pass(input);