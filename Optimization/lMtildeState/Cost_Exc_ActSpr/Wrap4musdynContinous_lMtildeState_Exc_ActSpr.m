function phaseout = Wrap4musdynContinous_lMtildeState_Exc_ActSpr(input)

persistent splinestruct
persistent isStancePhase

if isempty(splinestruct)|| size(splinestruct.MA,1) ~= length(input.phase.time) 
    splinestruct = SplineInputData(input.phase.time,input);
    isStancePhase = zeros(length(input.phase.time), 1);
    isStancePhase(find(input.phase.time < input.auxdata.pushoff_time)) = 1;
end

    disp('DEBUGwrap');
    disp(length(isStancePhase));
    % TODO detect if we should be setting isStancePhase at all.
input.auxdata.splinestruct = splinestruct;
input.auxdata.isStancePhase = isStancePhase;

phaseout = musdynContinous_lMtildeState_Exc_ActSpr(input);
