function phaseout = Wrap4continous_FtildeParamCalMultiPhase(input)

global splinestruct

numPhases = input.auxdata.numPhases;

for p = 1:numPhases
    
    if isempty(splinestruct(p).p) || size(splinestruct(p).p.MA,1) ~= length(input.phase(p).time)
        splinestruct(p).p = SplineInputData_MultiPhase(input.phase(p).time,input,p);
    end
    
    input.auxdata.splinestruct(p).p = splinestruct(p).p;
end


phaseout = continous_FtildeParamCalMultiPhase(input);