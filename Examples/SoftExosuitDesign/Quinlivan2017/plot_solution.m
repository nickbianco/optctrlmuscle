%% Calculate metabolic cost
clear all; close all; clc;
import org.opensim.modeling.*

cost=3;
switch cost
    case 1
        costdir = 'Exc_Act';
    case 2
        costdir = 'MinAlex';
    case 3
        costdir = 'Exc_Act_Hip_Shift';
    case 4
        costdir = 'Exc_Act_MinAlex';
end

study = 'Quinlivan2017';

load(fullfile(costdir,'ExoCurves.mat'))

cmap = zeros(11,3);
cmap(1,:) = [0 0 0];
cmap(2:end,:) = jet(10);

peak_force=am_peak.*(75/.6695);

for x=1:11

    filename=strcat('forceLevel',int2str(x-1),'.mat');
    load(fullfile(study,costdir,filename))

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
    T_exo = DatStore.T_exo;

    % Extract parts of the solution related to the device.
    control = OptInfo.result.solution.phase.control;
    state = OptInfo.result.solution.phase.state;

    % Get controls
    e       = control(:,1:numMuscles); e(e<0)=0;
    aT      = control(:,numMuscles+1:numMuscles+numDOFs);
    vMtilde = control(:,numMuscles+numDOFs+1:end);

    % Get states
    a       = state(:,1:numMuscles);
    lMtilde = state(:,numMuscles+1:end);

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
    
    % Muscle activations
    h1 = figure(1);
    for m = 1:numMuscles
        subplot(6,4,m)
        title(MuscleNames(m),'interpreter', 'none')
        if x==1
            plot(time,a(:,m),'k--','LineWidth',1.5)
        else
            plot(time,a(:,m),'Color',cmap(x,:),'LineWidth',1.2)
        end
        hold on
    end
    
    % Desired activations
    h7 = figure(7);
    muscles = {'med_gas_r','tib_ant_r','glut_med2_r','add_mag2_r','iliacus_r'};
    muscleNamesFull = {'Medial Gastrocnemius','Tibalis Anterior', ...
                       'Gluteus Medius','Adductor Magnus','Iliacus'};
    for i = 1:length(muscles)
        for m = 1:numMuscles
            if strcmp(muscles{i},MuscleNames(m))
                subplot(5,1,i)
                %ylabel(muscleNamesFull(i),'FontAngle','italic')
                hold on
                if x==1
                    plot(time,a(:,m),'k--','LineWidth',1.5)
                else
                    plot(time,a(:,m),'Color',cmap(x,:),'LineWidth',1.2)
                end
                axis([0.6 1.4 0 1])
            end
        end
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
    
    h4 = figure(4);
    if x>1
        subplot(2,1,1)
        if x==2
            plot(expTime,T_exp(:,1),'k-','LineWidth',1.75)
        end
        hold on
        plot(expTime,T_exo(:,1),'Color',cmap(x,:),'LineWidth',1.5)
        hold on
        title('Hip Moment','FontAngle','italic')
        box on
        ax = gca;
        ax.LineWidth = 1.5;
        ax.FontSize = 14;
        
        subplot(2,1,2)
        if x==2
            plot(expTime,T_exp(:,5),'k-','LineWidth',1.75)
        end
        hold on
        plot(expTime,T_exo(:,5),'Color',cmap(x,:),'LineWidth',1.5)
        hold on
        title('Ankle Moment','FontAngle','italic')
        box on
        ax = gca;
        ax.LineWidth = 1.5;
        ax.FontSize = 14;
    end
    
end

%% Plot bar graphs

for i=1:length(norm_average_wholebody_energy_rate)-1
    p_change(i)= (norm_average_wholebody_energy_rate(i+1)-...
    norm_average_wholebody_energy_rate(1))/norm_average_wholebody_energy_rate(1)*100;
end

figure;
paper_results=[-3.59,-6.45,-14.79,-22.83];

hold on
for i = 1:length(paper_results)
    h=bar(peak_force(i),paper_results(i));
    set(h,'FaceColor',[1 1 1]);
    set(h,'BarWidth',10);
    set(h,'LineWidth',1.5);
end

plot(peak_force,p_change,'k--','LineWidth',1.5)
box on
for i = 1:length(p_change)
    plot(peak_force(i),p_change(i),'o','MarkerSize',10,'MarkerEdgeColor',cmap(i+1,:),'MarkerFaceColor',cmap(i+1,:))
end

hold off
ylabel('Change in Metabolic Rate [%]')
xlabel('Peak Assistive Force [% BW]')
ax = gca;
ax.LineWidth = 1.5;
ax.FontSize = 14;
%% Plot all activations
% 
% figure;
% %list of muscle numbers to include
% muscles=1:24;
% %these do not have muscle parameters
% muscles(20)=[];
% muscles(17)=[];
% muscles(16)=[];
% %muscle only crosses knee
% muscles(9)=[];
% 
% p_stance = (time-0.6)./(1.43-0.6);
% cmap=jet(11);
% 
% for i=1:length(muscles)
%     subplot(5,4,i)
%     % remove _r and _ from titles
%     title(strrep(strrep(MuscleNames(muscles(i)),'_r',''),'_',''));
%     hold on
%     for j=1:11
%         plot(p_stance,activations(:,muscles(i),j),'Color',cmap(j,:))
%     end
% end
% 
% %% Plot hip flexors/extensors activations
% 
% figure;
% title('Hip Muscle Activations')
% %list of muscle numbers to include
% muscles=[14,15,19,4];
% p_stance = (time-0.6)./(1.43-0.6);
% cmap=jet(11);
% 
% for i=1:length(muscles)
%     subplot(1,4,i)
%     xlabel('Percent Stance Phase')
%     if i==1
%         ylabel('Activation')
%     end
%     % remove _r and _ from titles
%     title(strrep(strrep(MuscleNames(muscles(i)),'_r',''),'_',''));
%     hold on
%     for j=1:11
%         plot(p_stance,activations(:,muscles(i),j),'Color',cmap(j,:))
%     end
% end
% 
% %% Plot ankle flexors/extensors activations
% 
% figure;
% %list of muscle numbers to include
% muscles=[21,22,23,24];
% p_stance = (time-0.6)./(1.43-0.6);
% cmap=jet(11);
% 
% for i=1:length(muscles)
%     subplot(1,4,i)
%     xlabel('Percent Stance Phase')
%     if i==1
%         ylabel('Activation')
%     end
%     % remove _r and _ from titles
%     title(strrep(strrep(MuscleNames(muscles(i)),'_r',''),'_',''));
%     hold on
%     for j=1:11
%         plot(p_stance,activations(:,muscles(i),j),'Color',cmap(j,:))
%     end
% end
% 
% %% Plot high power activations
% 
% figure;
% %list of muscle numbers to include
% muscles=[7,12];
% p_stance = (time-0.6)./(1.43-0.6);
% cmap=jet(11);
% 
% for i=1:length(muscles)
%     subplot(1,2,i)
%     % remove _r and _ from titles
%     title(strrep(strrep(MuscleNames(muscles(i)),'_r',''),'_',''));
%     xlabel('Percent Stance Phase')
%     ylabel('Activation')
%     hold on
%     for j=1:11
%         plot(p_stance,activations(:,muscles(i),j),'Color',cmap(j,:),'LineWidth',1)
%     end
% end