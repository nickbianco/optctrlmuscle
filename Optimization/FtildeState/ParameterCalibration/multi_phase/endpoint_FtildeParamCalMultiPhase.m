function output = endpoint_FtildeParamCalMultiPhase(input)

numPhases = input.auxdata.numPhases;
obj = 0;
for p = 1:numPhases
    q = input.phase(p).integral;
    obj = obj + q;
end
output.objective = obj;


