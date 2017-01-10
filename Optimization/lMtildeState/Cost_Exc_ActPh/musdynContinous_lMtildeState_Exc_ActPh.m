function phaseout = musdynContinous_lMtildeState_Exc_ActPh(input)

% Get input data
NMuscles        = input.auxdata.NMuscles;
Ndof            = input.auxdata.Ndof;
tauAct          = input.auxdata.tauAct;
tauDeact        = input.auxdata.tauDeact;
params          = input.auxdata.params;
splinestruct    = input.auxdata.splinestruct;

for ip = 1:length(input.phase)
    numColPoints    = size(input.phase(ip).state,1);
    
    % Get controls
    e       = input.phase(ip).control(:,1:NMuscles);
    aT      = input.phase(ip).control(:,NMuscles+1:NMuscles+Ndof);
    vMtilde = input.phase(ip).control(:,NMuscles+Ndof+1:end);
    
    % Get states
    a       = input.phase(ip).state(:,1:NMuscles);
    lMtilde = input.phase(ip).state(:,NMuscles+1:end);
    
    % Get parameters
    if input.auxdata.ankle_clutched_spring
        % This actually has dimensions length(time) x 1.
        normSpringStiff = input.phase(ip).parameter(:, 1);
        springRestLength = input.phase(ip).parameter(:, 2);
    end

    % PATH CONSTRAINTS
    % Hill-equilibrium constraint
    [Hilldiff, F] = ForceEquilibrium_lMtildeState_Exc_ActPh(...
            a, lMtilde, vMtilde, splinestruct.phase(ip).LMT, params, ...
            input.auxdata.Fvparam, input.auxdata.Fpparam, input.auxdata.Faparam);

    % Moments constraint
    Topt = 150;
    % Only used if ankle_clutched_spring == true.
    maxSpringStiff = 400; % N-m/rad.
    Tdiff = zeros(numColPoints,Ndof);
    for dof = 1:Ndof
        T_exp=splinestruct.phase(ip).ID(:,dof);
        index_sel=(dof-1)*(NMuscles)+1:(dof-1)*(NMuscles)+NMuscles;
        T_sim=sum(F.*splinestruct.phase(ip).MA(:,index_sel),2) + Topt*aT(:,dof);
        if input.auxdata.ankle_clutched_spring
            % TODO phase-dependent.
            if any(dof == input.auxdata.clutched_spring_dofs)
                ankleAngle = -(splinestruct.phase(ip).IK(:,dof) - springRestLength);
                T_sim = T_sim + maxSpringStiff * normSpringStiff .* ankleAngle;
            end
        end
        Tdiff(:,dof) =  (T_exp-T_sim);
    end

    phaseout(ip).path = [Tdiff Hilldiff];

    % DYNAMIC CONSTRAINTS
    % Activation dynamics
    dadt = ones(numColPoints,NMuscles);
    for m = 1:NMuscles
        dadt(:,m) = ActivationDynamics(e(:,m),a(:,m),tauAct(m),tauDeact(m),input.auxdata.b);
    end

    % Contraction dynamics is implicit
    dlMtildedt = 10*vMtilde;
    
    phaseout(ip).dynamics = [dadt dlMtildedt];
    
    % OBJECTIVE FUNCTION
    w1 = 1000;
    phaseout(ip).integrand = sum(e.^2,2) + sum(a.^2,2) + w1.*sum(aT.^2,2);
end









