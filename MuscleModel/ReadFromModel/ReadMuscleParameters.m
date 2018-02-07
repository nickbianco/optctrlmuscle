function [params,lOpt,L_TendonSlack,Fiso,PennationAngle,metabolicParams]=ReadMuscleParameters(ModelPath,names)
% input= path to model and cell array with muscle names
% output= params (5xNMuscles) with  row: (1)  IsomForce (2)OptFiberLength
% 			(3) TendonSlackLength (4) PennationAngle (5) MaxFiberVelocity

% read the model
import org.opensim.modeling.*;
model = Model(ModelPath);

% read the muscle properties
nom = length(names);
params = zeros(5, nom);			% pre allocate
muscles = model.getMuscles();

% read the muscle parameters necessary for calculating metabolic cost
metabolicParams = zeros(5, nom);

probeSet = model.getProbeSet();
probe = probeSet.get('metabolic_power');
probeUmberger = Umberger2010MuscleMetabolicsProbe.safeDownCast(probe);

for i = 1:nom
   muscle = muscles.get(names{i});
   params(3,i) = muscle.getTendonSlackLength();		
   params(2,i) = muscle.getOptimalFiberLength(); 	
   params(1,i) = muscle.getMaxIsometricForce();  	
   params(4,i) = muscle.getPennationAngleAtOptimalFiberLength(); 
   params(5,i) = muscle.getMaxContractionVelocity()*params(2,i);
   
   % proportion of fast twitch fibers
   rST = probeUmberger.getRatioSlowTwitchFibers(names{i});
   metabolicParams(1,i) = 1 - rST;
   % optimal fiber length
   metabolicParams(2,i) = muscle.getOptimalFiberLength();
   % max contractile velocity (Lceopt/s)
   metabolicParams(3,i) = muscle.getMaxContractionVelocity();
   % muscle mass
   sigma = probeUmberger.getSpecificTension(names{i}); % Specific tension [N/m^2]
   PCSA = params(1,i)/sigma;      % Physiological cross sectional area [m^2]
   rho = 1059.7; % Muscle density [kg/m^3]
   mass = PCSA*rho*params(2,i); % Muscle mass [kg]
   metabolicParams(4,i) = mass;
   % aerobic vs. anaerobic scaling factor
   metabolicParams(5,i) = 1.0;
end

% create additional variables with the same information
Fiso=params(1,:);
lOpt=params(2,:);
L_TendonSlack=params(3,:);
PennationAngle=params(4,:);
end