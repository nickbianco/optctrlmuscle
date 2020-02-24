function phaseout = continous_FtildeExoTopology_Met_ActHesWrap(input)

persistent splinestruct

if isempty(splinestruct) || size(splinestruct.MA,1) ~= length(input.phase.time.f)    
    splinestruct = SplineInputData(input.phase.time.f,input);
end

input.auxdata.splinestruct = splinestruct;

phaseout = continous_FtildeExoTopology_Met_ActADiGatorHes(input);