function phaseout = Wrap4continous_Ftilde_vAAHE_Met(input)

% input is good first time here
persistent splinestruct
% then splinestruct is not fine

if isempty(splinestruct)|| size(splinestruct.MA,1) ~= length(input.phase.time) 
    splinestruct = SplineInputData(input.phase.time,input);
end

input.auxdata.splinestruct = splinestruct;

phaseout = continous_Ftilde_vAAHE_Met(input);