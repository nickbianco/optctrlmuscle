function phaseout = musdynContinous_lMtildeState_Exc_ActSprGrdWrap(input)

persistent splinestruct

if isempty(splinestruct)|| size(splinestruct.MA,1) ~= length(input.phase.time.f)    
    splinestruct = SplineInputData(input.phase.time.f,input);
end

input.auxdata.splinestruct = splinestruct;

phaseout = musdynContinous_lMtildeState_Exc_ActSprADiGatorGrd(input);
