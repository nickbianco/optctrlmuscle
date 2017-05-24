close all; clear all; clc;
import org.opensim.modeling.*

cost = 3;
switch cost
    case 1
        costdir = 'Exc_Act';
    case 2
        costdir = 'MinAlex';
    case 3
        costdir = 'Exc_Act_MinAlex';
end

DingExoCurves = load('DingExoCurves.mat');
for c = 1:7
    
    conds = {'slack','esep','eslp','lsep','lslp','esmp','lsmp'};
    condActual = ([0 -0.35 -0.5 -0.37 -0.42 -5.92 -5.92]+5.92)/5.92;
    condColor = [64 64 64; 241 102 69; 255 198 93; 152 204 103; 76 195 217; 233 0 255; 21 0 255]/255;
    condName = {'UNPD','ESEP','ESLP','LSEP','LSLP','ESMP','LSMP'};
    load(fullfile(costdir,[conds{c} '.mat'])) 
    
    numDOFs = DatStore.nDOF;
    numMuscles = DatStore.nMuscles;
    DOFNames = DatStore.DOFNames;
    MuscleNames = DatStore.MuscleNames;
    
    time = OptInfo.result.solution.phase.time;
    numColPoints = length(time);
    
    auxdata = OptInfo.result.setup.auxdata;
    
    % Extract experimental data.
    expTime = DatStore.time;
    qExp = DatStore.q_exp;
    T_exp = DatStore.T_exp;
    momArmsExp = DatStore.dM;
    momArms = interp1(expTime, momArmsExp, time);
    jointAngles = pi / 180. * interp1(expTime, qExp, time);
    for k = 1:numDOFs
        if strcmp(DOFNames(k),'hip_flexion_r')
            hipFlexIdx = k;
        end
    end
    
    % Extract parts of the solution related to the device.
    control = OptInfo.result.solution.phase.control;
    state = OptInfo.result.solution.phase.state;
    
    % Get controls
    e       = control(:,1:numMuscles); e(e<0)=0;
    aT      = control(:,numMuscles+1:numMuscles+numDOFs);
    if strcmp(conds{c},'DingOpt')
        vMtilde = control(:,numMuscles+numDOFs+1:end-1);
        aD = control(:,end);
    else
        vMtilde = control(:,numMuscles+numDOFs+1:end);
    end
    
    % Get states
    a       = state(:,1:numMuscles);
    lMtilde = state(:,numMuscles+1:end);
    
    % Metabolic cost
    modelApoorva = Model('Rajagopal2015.osim');
    musclesApoorva = modelApoorva.getMuscles();
    % pect_r, quad_fem_r and gem_r are not in Apoorva model
    % all replaced by omit
    MuscleNamesApoorva = {'glut_med1_r','glut_med2_r','glut_med3_r',...
        'bifemlh_r','bifemsh_r','sar_r','add_mag2_r','tfl_r','omit',...
        'grac_r','glut_max1_r','glut_max2_r','glut_max3_r','iliacus_r',...
        'psoas_r','omit','omit','peri_r','rect_fem_r','vas_int_r'...
        'med_gas_r','soleus_r','tib_post_r','tib_ant_r'};
    muscleMap = containers.Map(MuscleNames,MuscleNamesApoorva);
    
    probeSet = modelApoorva.getProbeSet();
    probe = probeSet.get('metabolic_power');
    probeUmberger = Umberger2010MuscleMetabolicsProbe.safeDownCast(probe);
    
    rho = 1059.7; % Muscle density [kg/m^3]
    maxFiberVel = 12;  % Fiber-lengths per second
    
    lMT = interp1(DatStore.time,DatStore.LMT,time);
    [F,Fiso] = calcMuscleForcesDeGroote(a,lMtilde,vMtilde,lMT,auxdata);
    
    musc_energy_rate = NaN(numColPoints,numMuscles);
    for m = 1:numMuscles
        muscleNameApoorva = muscleMap(MuscleNames{m});
        if string(muscleNameApoorva)=='omit'
            % leave as NaN
        else
            musc = musclesApoorva.get(muscleNameApoorva);
            Fmax = musc.getMaxIsometricForce;   % Max isometric force [N]
            Lceopt = musc.getOptimalFiberLength;         % Optimal fiber length [m]
            
            rST = probeUmberger.getRatioSlowTwitchFibers(muscleNameApoorva);
            param_rFT = 1 - rST;        % Proportion of fast-twitch muscle fibers
            param_Arel = 0.1 + 0.4*param_rFT;
            param_Brel = param_Arel*maxFiberVel;
            
            sigma = probeUmberger.getSpecificTension(muscleNameApoorva); % Specific tension [N/m^2]
            PCSA = Fmax/sigma;      % Physiological cross sectional area [m^2]
            mass = PCSA*rho*Lceopt; % Muscle mass [kg]
            
            paramsUmb = struct('Lceopt',Lceopt, 'Arel',param_Arel, ...
                'Brel',param_Brel, 'Fmax',Fmax, 'rFT',param_rFT, ...
                'VceMax_LceoptsPerSecond',param_Brel/param_Arel, ...
                'muscleMass',mass, 'scalingFactorS',1.0, ...
                'versionNumber',2010);
            VCEmax_mps = paramsUmb.VceMax_LceoptsPerSecond * Lceopt; % [m/s]
            
            heatRates = NaN(numColPoints,5);
            for i = 1:numColPoints
                Lce = lMtilde(i,m)*Lceopt;
                Vce = vMtilde(i,m)*VCEmax_mps;
                heatRates(i,:) = calcUmbergerProbe(Lce,Vce,F(i,m),Fiso(i,m),e(i,m),a(i,m),paramsUmb);
            end
            
            musc_energy_rate(:,m) = heatRates(:,5) * mass;
        end
    end
    
    bodyMass = 75; % kg
    wholebody_energy_rate = nansum(musc_energy_rate,2);
    norm_average_wholebody_energy_rate(c) = mean(wholebody_energy_rate) / bodyMass
    if c==1
        scale = norm_average_wholebody_energy_rate(c);
    end
    
    % Muscle activations
    h1 = figure(1);
    for m = 1:numMuscles
        subplot(6,4,m)
        title(MuscleNames(m),'interpreter', 'none')
        plot(time,a(:,m),'Color',condColor(c,:),'LineWidth',1.2)
        hold on
    end
    
    % Metabolic rate
    h2 = figure(2);
    bar(c,condActual(c),'FaceColor',condColor(c,:))
    hold on
    axis([0 8 0.85 1.15])
    
    % Normalized fiber lengths
    h3 = figure(3);
    for m = 1:numMuscles
        subplot(6,4,m)
        title(MuscleNames(m),'interpreter', 'none')
        plot(time,lMtilde(:,m),'Color',condColor(c,:),'LineWidth',1.2)
        hold on
    end
    
    % Normalized fiber velocities
    h4 = figure(4);
    for m = 1:numMuscles
        subplot(6,4,m)
        title(MuscleNames(m),'interpreter', 'none')
        plot(time,vMtilde(:,m),'Color',condColor(c,:),'LineWidth',1.2)
        hold on
    end
    
    % Reserve actuators
    h5 = figure(5);
    for k = 1:numDOFs
        subplot(3,2,k)
        title(DOFNames(k),'interpreter', 'none')
        plot(time,aT(:,k),'Color',condColor(c,:),'LineWidth',1.2)
        hold on
    end
    
    % Assistive moments
    h6 = figure(6);
    if strcmp(conds{c},'DingOpt')
        Topt_exo = DatStore.Topt_exo;
        peakHipExtMoment = Topt_exo(hipFlexIdx);        
        plot(time,-peakHipExtMoment*aD)
    else
        if ~strcmp(conds{c},'slack')
%             exoForce = DingExoCurves.(conds{c}).F;
%             exoMomentArm = DingExoCurves.(conds{c}).r;
%             exoTime = DingExoCurves.time;
            plot(DatStore.time,-DatStore.T_exo(:,hipFlexIdx),'Color',condColor(c,:),'LineWidth',1.5)
            hold on
        end
    end
end

set(h1,'Name','Muscle Activations')
set(h1,'units','normalized','position',[0 0 1 1])
print(h1,'Ding_Gait2354_MuscleActivations','-dmeta')

set(h2,'Name','Metabolic Rate')
set(h2,'units','normalized','position',[0.25 0.25 0.75 0.75])
figure(2)
plot(1:7,norm_average_wholebody_energy_rate/scale,'o--','LineWidth',1.5,...
    'Color',[72/255 0 1])
set(gca,'XTick',1:7,'XTickLabels',condName)
ylabel('Normalized Metabolic Rate')
print(h2,'Ding_Gait2354_MetabolicRate','-dmeta')

set(h3,'Name','Normalized Fiber Length')
set(h3,'units','normalized','position',[0 0 1 1])
print(h3,'Ding_Gait2354_NormalizedFiberLength','-dmeta')

set(h4,'Name','Normalized Fiber Velocity')
set(h4,'units','normalized','position',[0 0 1 1])
print(h4,'Ding_Gait2354_NormalizedFiberVelocity','-dmeta')

set(h5,'Name','Reserve Actuators')
set(h5,'units','normalized','position',[0 0 1 1])
print(h5,'Ding_Gait2354_ReserveActuators','-dmeta')

set(h6,'Name','Assistive Moment')
figure(6)
plot(expTime,-T_exp(:,hipFlexIdx),'k--','LineWidth',1.5)
dq = diff(qExp);
plot(expTime(1:end-1),T_exp(1:end-1,hipFlexIdx).*dq(:,hipFlexIdx),'--','LineWidth',1.5,'Color',[179 0 0]/255)
plot(expTime(1:end-1),-50*dq(:,hipFlexIdx),'--','LineWidth',1.5,'Color',[51 102 0]/255)
legend({'ESEP','ESLP','LSEP','LSLP','ESMP','LSMP','Hip ID Moment','Hip Power','Hip Velocity (scaled)'},'Location','Northeast');
axis([0.6 1.4 -50 100])
xlabel('Stance phase (s)')
ylabel('Moment (N/m)')
print(h6,'Ding_Gait2354_AssistiveMoments','-dmeta')











