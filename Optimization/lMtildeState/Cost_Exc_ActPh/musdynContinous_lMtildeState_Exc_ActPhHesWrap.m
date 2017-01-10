function phaseout = musdynContinous_lMtildeState_Exc_ActPhHesWrap(input)

persistent splinestruct

if isempty(splinestruct)
    splinestruct = struct();
end

for ip = 1:length(input.phase)
    if ~isfield(splinestruct, 'phase') || length(splinestruct.phase) < ip || ...
            size(splinestruct.phase(ip).MA,1) ~= length(input.phase(ip).time.f)
        splinestruct.phase(ip) = SplineInputData(input.phase(ip).time.f,input);
    end
end

input.auxdata.splinestruct = splinestruct;

phaseout = musdynContinous_lMtildeState_Exc_ActPhADiGatorHes(input);
