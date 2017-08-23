import org.opensim.modeling.*

addpath(genpath('../../../Analysis/'))

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

plotDOFs=[1];
n=1;

for idof = plotDOFs
    subplot(length(plotDOFs),1,n);
    hold on;
    plot(expTime, DatStore.T_exp(:, idof), 'k', 'LineWidth', 2);
    legendEntries = {'net'};
    sumMoment = zeros(length(TForce(:, 1)), 1);
    for imusc = 1:numMuscles
        if any(momArms(:, idof, imusc)) > 0.00001
            thisMoment = TForce(:, imusc) .* momArms(:, idof, imusc);
            plot(time(1:end-1), thisMoment(1:end-1));
            legendEntries = [legendEntries {MuscleNames{imusc}}];
            sumMoment = sumMoment + thisMoment;
        end
    end
    deviceMoment = DatStore.ExoTorques(:,idof);
    plot(time,deviceMoment, 'r', 'LineWidth', 2);
    legendEntries = [legendEntries {'device'}];
%     sumMoment = sumMoment + deviceMoment; 
%     plot(time(1:end-1), sumMoment(1:end-1), 'r', 'LineWidth', 2);
%     legendEntries = [legendEntries {'sum'}];
    legend(legendEntries, 'Interpreter', 'none');
    title(DatStore.DOFNames{idof}, 'Interpreter', 'none');
    if idof == numDOFs
        xlabel('time (s)');
    end
    ylabel('moment (N-m)');
    n=n+1;
end

% Metabolic cost
modelApoorva = Model('Rajagopal2015.osim');
musclesApoorva = modelApoorva.getMuscles();
% pect_r, quad_fem_r, gem_r and per_tert_r are not in Apoorva model
% all replaced by omit
% for 3 dof on right leg remove glut_med1_r
MuscleNamesApoorva = {'glut_med1_r','glut_med2_r','glut_med3_r','glut_min1_r','glut_min2_r'...
    'glut_min3_r','semimem_r','semiten_r','bifemlh_r','bifemsh_r','sar_r'...
    'add_long_r','add_brev_r','add_mag1_r','add_mag2_r','add_mag3_r','tfl_r'...
    'omit','grac_r','glut_max1_r','glut_max2_r','glut_max3_r','iliacus_r'...
    'psoas_r','omit','omit','peri_r','rect_fem_r','vas_med_r','vas_int_r'...
    'vas_lat_r','med_gas_r','lat_gas_r','soleus_r','tib_post_r','flex_dig_r'...
    'flex_hal_r','tib_ant_r','per_brev_r','per_long_r','omit','ext_dig_r'...
    'ext_hal_r'};
% MuscleNamesApoorva = {'glut_med2_r','glut_med3_r','glut_min1_r','glut_min2_r'...
%     'glut_min3_r','semimem_r','semiten_r','bifemlh_r','bifemsh_r','sar_r'...
%     'add_long_r','add_brev_r','add_mag1_r','add_mag2_r','add_mag3_r','tfl_r'...
%     'omit','grac_r','glut_max1_r','glut_max2_r','glut_max3_r','iliacus_r'...
%     'psoas_r','omit','omit','peri_r','rect_fem_r','vas_med_r','vas_int_r'...
%     'vas_lat_r','med_gas_r','lat_gas_r','soleus_r','tib_post_r','flex_dig_r'...
%     'flex_hal_r','tib_ant_r','per_brev_r','per_long_r','omit','ext_dig_r'...
%     'ext_hal_r','glut_med1_l','glut_med2_l','glut_med3_l','glut_min1_l','glut_min2_l'...
%     'glut_min3_l','semimem_l','semiten_l','bifemlh_l','bifemsh_l','sar_l'...
%     'add_long_l','add_brev_l','add_mag1_l','add_mag2_l','add_mag3_l','tfl_l'...
%     'omit','grac_l','glut_max1_l','glut_max2_l','glut_max3_l','iliacus_l'...
%     'psoas_l','omit','omit','peri_l','rect_fem_l','vas_med_l','vas_int_l'...
%     'vas_lat_l','med_gas_l','lat_gas_l','soleus_l','tib_post_l','flex_dig_l'...
%     'flex_hal_l','tib_ant_l','per_brev_l','per_long_l','omit','ext_dig_l'...
%     'ext_hal_l'};  

muscleMap = containers.Map(MuscleNames,MuscleNamesApoorva);

probeSet = modelApoorva.getProbeSet();
probe = probeSet.get('metabolic_power');
probeUmberger = Umberger2010MuscleMetabolicsProbe.safeDownCast(probe);

rho = 1059.7; % Muscle density [kg/m^3]
maxFiberVel = 12;  % Fiber-lengths per second

lMT = interp1(DatStore.time,DatStore.LMT,time);
[~, ~, F, Fiso] = DeGroote2016Muscle_lMtildeState(a,lMtilde,vMtilde, ...
                         lMT,auxdata.params,auxdata.Fvparam, ...
                         auxdata.Fpparam,auxdata.Faparam);

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
            heatRates(i,:) = calcUmbergerProbe(Lce,Vce,F(i,m),Fiso(i,m), ... 
                                               e(i,m),a(i,m),paramsUmb);
        end

        musc_energy_rate(:,m) = heatRates(:,5) * mass;
    end
end

bodyMass = 75; % kg
wholebody_energy_rate = nansum(musc_energy_rate,2);
norm_average_wholebody_energy_rate = mean(wholebody_energy_rate) / bodyMass

% Find top 15 muscle activations (based on sum of squares)
a_sum_sq = sum(a.^2,1);
[sortedX,sortingInd] = sort(a_sum_sq,'descend');
maxValues=sortedX(1:15);
maxInd=sortingInd(1:15);

% Muscle Activations
figure
muscles = DatStore.MuscleNames(maxInd);

for i = 1:length(muscles)
    for m = 1:numMuscles
        if strcmp(muscles{i},MuscleNames(m))
            subplot(5,3,i)
            title(muscles{i},'Interpreter','none')
            hold on
            plot(time,a(:,m),'LineWidth',1.2)
            axis([0.05 0.98 0 1])
        end
    end
end

