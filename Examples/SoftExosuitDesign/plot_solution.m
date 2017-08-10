%% Calculate metabolic cost
import org.opensim.modeling.*
load ExoCurves.mat

cost=1;
switch cost
    case 1
        costdir = 'Exc_Act';
    case 2
        costdir = 'MinAlex';
    case 3 % Add mass
        costdir = 'Exc_Act';
end

cmap = zeros(11,3);
cmap(1,:) = [0 0 0];
cmap(2:end,:) = jet(10);

for x=1:11

    % add back 'HipAnkle', before costdir? 
    
    filename=strcat('forceLevel',int2str(x-1),'.mat');
    load(fullfile('HipAnkle',costdir,filename))

    numDOFs = DatStore.nDOF;
    numMuscles = DatStore.nMuscles;

    time = OptInfo.result.solution.phase.time;
    numColPoints = length(time);

    auxdata = OptInfo.result.setup.auxdata;

    % Extract experimental data.
    expTime = DatStore.time;
    qExp = DatStore.q_exp;
    momArmsExp = DatStore.dM;
    momArms = interp1(expTime, momArmsExp, time);
    jointAngles = pi / 180. * interp1(expTime, qExp, time);
    T_exp = DatStore.T_exp;
    Fopt_exo = DatStore.Fopt_exo;
    ankleIdx = strcmp('ankle_angle_r',DatStore.DOFNames);
    anklePeakForce = Fopt_exo(ankleIdx);
    ankleID = T_exp(:,ankleIdx);
    hipIdx = strcmp('hip_flexion_r',DatStore.DOFNames);
    hipPeakForce = Fopt_exo(hipIdx);
   	hipID = T_exp(:,hipIdx);
    
    % Interpolate inverse dynamics moments
    timeID = linspace(0.6,1.4,length(hipID));
    hipID = interp1(timeID,hipID,time);
    ankleID = interp1(timeID,ankleID,time);
    
    % Extract parts of the solution related to the device.
    control = OptInfo.result.solution.phase.control;
    state = OptInfo.result.solution.phase.state;

    % Get controls
    e       = control(:,1:numMuscles); e(e<0)=0; e(e>1)=1;
    aT      = control(:,numMuscles+1:numMuscles+numDOFs);
    vMtilde = control(:,numMuscles+numDOFs+1:end-1);
    aD      = control(:,end);

    % Get states
    a       = state(:,1:numMuscles);
    lMtilde = state(:,numMuscles+1:end);
    
    % Get parameter
    alpha = OptInfo.result.solution.parameter;
    tradeoff = auxdata.tradeoff;

    % Joint moment breakdown.
    deviceIndices = strmatch('ankle_angle', DatStore.DOFNames);
    assert(length(deviceIndices) == 1);

    % Metabolic cost
    modelApoorva = Model('../Rajagopal2015.osim');
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
    norm_average_wholebody_energy_rate(x) = mean(wholebody_energy_rate) / bodyMass;
     
    % Mass properties of suit
    if cost==3
        waist_mass = 4.9+0.433;
        shank_mass = 0.356;
        foot_mass = 0.364;
        
        % equations from Browning 2007
        % units are watts/kg
        waist_cost = (0.045*waist_mass);
        shank_cost = (0.076*shank_mass);
        foot_cost = (0.2*foot_mass);
        device_cost = waist_cost+shank_cost+foot_cost;
        
        % calculate new metabolic cost
        norm_average_wholebody_energy_rate(x) = mean(wholebody_energy_rate) / bodyMass + device_cost;
    end
        
    % Muscle activations
    h1 = figure(1);
    for m = 1:numMuscles
        subplot(6,4,m)
        title(MuscleNames(m),'interpreter', 'none')
        plot(time,a(:,m),'Color',cmap(x,:),'LineWidth',1.2)
        hold on
    end
    
    % Normalized fiber length
    h2 = figure(2);
    for m = 1:numMuscles
        subplot(6,4,m)
        title(MuscleNames(m),'interpreter', 'none')
        plot(time,lMtilde(:,m),'Color',cmap(x,:),'LineWidth',1.2)
        hold on
    end
    
    % Normalized fiber velocity
    h3 = figure(3);
    for m = 1:numMuscles
        subplot(6,4,m)
        title(MuscleNames(m),'interpreter', 'none')
        plot(time,vMtilde(:,m),'Color',cmap(x,:),'LineWidth',1.2)
        hold on
    end
    
    % Device control and tradeoff parameter solution
    r = 0.1;
    h4 = figure(4);
    subplot(2,2,1)
    bar(x,alpha)
    hold on
    title('Tradeoff parameter')
    axis([0 12 -1 1])
    
    subplot(2,2,3)
    plot(time,aD,'Color',cmap(x,:),'LineWidth',1.5)
    hold on
    title('Device control')
    
    subplot(2,2,2)
    plot(time,hipPeakForce*aD*r*(1+tradeoff(hipIdx)*alpha),'Color',cmap(x,:),'LineWidth',1.5)
    hold on
    plot(time,hipID,'k--')
    title('Hip Moment')
    
    subplot(2,2,4)
    plot(time,anklePeakForce*aD*r*(1+tradeoff(ankleIdx)*alpha),'Color',cmap(x,:),'LineWidth',1.5)
    hold on
    plot(time,ankleID,'k--')
    title('Ankle Moment')
end

folder = [Misc.costfun '_Hip_Shift'];

%% Plot bar graphs

for i=1:length(norm_average_wholebody_energy_rate)-1
    p_change(i)= (norm_average_wholebody_energy_rate(i+1)-...
    norm_average_wholebody_energy_rate(1))/norm_average_wholebody_energy_rate(1)*100;
end

figure;
paper_results=[-3.59,-6.45,-14.79,-22.83];
color=[[0.498, 0.396, 0.635];[0.258, 0.670, 0.784];[0.596, 0.8, 0.403];[0.968, 0.6, 0.215]];
ind=am_peak.*(75/.6695);

hold on
for i = 1:length(paper_results)
    h=bar(ind(i),paper_results(i));
    set(h,'FaceColor',color(i,:));
    set(h,'BarWidth',10);
end

plot(ind,p_change,'k--')
plot(ind,p_change,'k.','MarkerSize',15)
hold off
ylabel('Change in Metabolic Rate [%]')
xlabel('Peak Assistive Force [% BW]')
title('Reduction in Net Metabolic Rate')