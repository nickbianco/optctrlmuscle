function plotDingExpData
clear all; close all; clc;

dir = 'C:\Users\Nick\Projects\AssistedHipExtension\data\mat';
subj = strcat('subject',{'56','74','78','79','81','82'},'.mat');
subjMass = [85 90.9 85 80 81.2 73.2]; 

for s = 1:length(subj)
   data = load(fullfile(dir,subj{s}));
   
   % ESEP (prf_2)
   esepRForce(:,s) = data.prf_2.ahe_for_r.x.ave;
   esepRArm(:,s) = data.prf_2.hip_arm_r.x.ave;
   esepLForce(:,s) = data.prf_2.ahe_for_r.x.ave;
   esepLArm(:,s) = data.prf_2.hip_arm_r.x.ave;
   esepRHipTor(:,s) = data.prf_2.hip_tor_r.x.ave * subjMass(s);
   esepLHipTor(:,s) = data.prf_2.hip_tor_l.x.ave * subjMass(s);
   esepRHipPow(:,s) = data.prf_2.hip_pow_r.x.ave * subjMass(s);
   esepLHipPow(:,s) = data.prf_2.hip_pow_l.x.ave * subjMass(s);
   esepRHipAngVel(:,s) = data.prf_2.hip_ang_vel_r.x.ave * subjMass(s);
   esepLHipAngVel(:,s) = data.prf_2.hip_ang_vel_l.x.ave * subjMass(s);
   
   % ESLP (prf_3)
   eslpRForce(:,s) = data.prf_3.ahe_for_r.x.ave;
   eslpRArm(:,s) = data.prf_3.hip_arm_r.x.ave;
   eslpLForce(:,s) = data.prf_3.ahe_for_r.x.ave;
   eslpLArm(:,s) = data.prf_3.hip_arm_r.x.ave;
   eslpRHipTor(:,s) = data.prf_3.hip_tor_r.x.ave * subjMass(s);
   eslpLHipTor(:,s) = data.prf_3.hip_tor_l.x.ave * subjMass(s);
   eslpRHipPow(:,s) = data.prf_3.hip_pow_r.x.ave * subjMass(s);
   eslpLHipPow(:,s) = data.prf_3.hip_pow_l.x.ave * subjMass(s);
   eslpRHipAngVel(:,s) = data.prf_3.hip_ang_vel_r.x.ave * subjMass(s);
   eslpLHipAngVel(:,s) = data.prf_3.hip_ang_vel_l.x.ave * subjMass(s);
   
   % LSEP (prf_1)
   lsepRForce(:,s) = data.prf_1.ahe_for_r.x.ave;
   lsepRArm(:,s) = data.prf_1.hip_arm_r.x.ave;
   lsepLForce(:,s) = data.prf_1.ahe_for_r.x.ave;
   lsepLArm(:,s) = data.prf_1.hip_arm_r.x.ave;
   lsepRHipTor(:,s) = data.prf_1.hip_tor_r.x.ave * subjMass(s);
   lsepLHipTor(:,s) = data.prf_1.hip_tor_l.x.ave * subjMass(s);
   lsepRHipPow(:,s) = data.prf_1.hip_pow_r.x.ave * subjMass(s);
   lsepLHipPow(:,s) = data.prf_1.hip_pow_l.x.ave * subjMass(s);
   lsepRHipAngVel(:,s) = data.prf_1.hip_ang_vel_r.x.ave * subjMass(s);
   lsepLHipAngVel(:,s) = data.prf_1.hip_ang_vel_l.x.ave * subjMass(s);
   
   % LSLP (prf_4)
   lslpRForce(:,s) = data.prf_4.ahe_for_r.x.ave;
   lslpRArm(:,s) = data.prf_4.hip_arm_r.x.ave;
   lslpLForce(:,s) = data.prf_4.ahe_for_r.x.ave;
   lslpLArm(:,s) = data.prf_4.hip_arm_r.x.ave;
   lslpRHipTor(:,s) = data.prf_4.hip_tor_r.x.ave * subjMass(s);
   lslpLHipTor(:,s) = data.prf_4.hip_tor_l.x.ave * subjMass(s);
   lslpRHipPow(:,s) = data.prf_4.hip_pow_r.x.ave * subjMass(s);
   lslpLHipPow(:,s) = data.prf_4.hip_pow_l.x.ave * subjMass(s);
   lslpRHipAngVel(:,s) = data.prf_4.hip_ang_vel_r.x.ave * subjMass(s);
   lslpLHipAngVel(:,s) = data.prf_4.hip_ang_vel_l.x.ave * subjMass(s);
  
end

h1 = figure(1);
subplot(1,2,1)
plotData(1,'ESEP',esepLForce,esepLArm,esepLHipTor,esepLHipPow,esepLHipAngVel)
title('ESEP: Left Leg')
subplot(1,2,2)
plotData(1,'ESEP',esepRForce,esepRArm,esepRHipTor,esepRHipPow,esepRHipAngVel)
title('ESEP: Right Leg')

h2 = figure(2);
subplot(1,2,1)
plotData(2,'ESLP',eslpLForce,eslpLArm,eslpLHipTor,eslpLHipPow,eslpLHipAngVel)
title('ESLP: Left Leg')
subplot(1,2,2)
plotData(2,'ESLP',eslpRForce,eslpRArm,eslpRHipTor,eslpRHipPow,eslpRHipAngVel)
title('ESLP: Right Leg')

h3 = figure(3);
subplot(1,2,1)
plotData(3,'LSEP',lsepLForce,lsepLArm,lsepLHipTor,lsepLHipPow,lsepLHipAngVel)
title('LSEP: Left Leg')
subplot(1,2,2)
plotData(3,'LSEP',lsepRForce,lsepRArm,lsepRHipTor,lsepRHipPow,lsepRHipAngVel)
title('LSEP: Right Leg')

h4 = figure(4);
subplot(1,2,1)
plotData(4,'LSLP',lslpLForce,lslpLArm,lslpLHipTor,lslpLHipPow,lslpLHipAngVel)
title('LSLP: Left Leg')
subplot(1,2,2)
plotData(4,'LSLP',lslpRForce,lslpRArm,lslpRHipTor,lslpRHipPow,lslpRHipAngVel)
title('LSLP: Right Leg')

% Individual subjects
for s = 1:length(subj)
    fig = s+4;
    figure(fig);
    subplot(1,2,1)
    plotData(fig,'ESEP',esepLForce(:,s),esepLArm(:,s),esepLHipTor(:,s),esepLHipPow(:,s),esepLHipAngVel(:,s)); hold on;
    plotData(fig,'ESLP',eslpLForce(:,s),eslpLArm(:,s),eslpLHipTor(:,s),eslpLHipPow(:,s),eslpLHipAngVel(:,s)); hold on;
    plotData(fig,'LSEP',lsepLForce(:,s),lsepLArm(:,s),lsepLHipTor(:,s),lsepLHipPow(:,s),lsepLHipAngVel(:,s)); hold on;
    plotData(fig,'LSLP',lslpLForce(:,s),lslpLArm(:,s),lslpLHipTor(:,s),lslpLHipPow(:,s),lslpLHipAngVel(:,s))
    title([subj{s} ': Left Leg'])
    
    subplot(1,2,2)
    plotData(fig,'ESEP',esepRForce(:,s),esepRArm(:,s),esepRHipTor(:,s),esepRHipPow(:,s),esepRHipAngVel(:,s)); hold on;
    plotData(fig,'ESLP',eslpRForce(:,s),eslpRArm(:,s),eslpRHipTor(:,s),eslpRHipPow(:,s),eslpRHipAngVel(:,s)); hold on;
    plotData(fig,'LSEP',lsepRForce(:,s),lsepRArm(:,s),lsepRHipTor(:,s),lsepRHipPow(:,s),lsepRHipAngVel(:,s)); hold on;
    plotData(fig,'LSLP',lslpRForce(:,s),lslpRArm(:,s),lslpRHipTor(:,s),lslpRHipPow(:,s),lslpRHipAngVel(:,s))
    title([subj{s} ': Right Leg'])
    
end

end

function plotData(fig,prf,Force,Arm,HipTor,HipPow,HipAngVel)

condColors = [241 102 69; 255 198 93; 152 204 103; 76 195 217]/255;
switch prf
    case 'ESEP'
        condColor = condColors(1,:);
    case 'ESLP'
        condColor = condColors(2,:);
    case 'LSEP'
        condColor = condColors(3,:);
    case 'LSLP'
        condColor = condColors(4,:);
end


figure(fig)
phase = 1:500;
HipDeviceMom = Force.*Arm;
plotCurve(HipDeviceMom(phase,:),condColor,condColor,phase,'-')
hold on
plotCurve(HipTor(phase,:),condColor,condColor,phase,'--')
plotCurve(HipPow(phase,:),condColor,condColor,phase,'-.')
plotCurve(-HipAngVel(phase,:)/100,condColor,condColor,phase,':')
%legend('Device Moment','Hip Moment','Hip Power','Hip Angular Velocity','Location','best')

hold off

end

function plotCurve(qty,color,condColor,phase,linestyle)

meanQty = mean(qty,2); 
stdQty = std(qty,0,2);
plot(meanQty,linestyle,'LineWidth',1.5,'Color',color)
hold on
%meanInterp = interp1(1:length(phase),meanQty,1:50:length(phase));
%plot(1:50:length(phase),meanInterp,'o','MarkerSize',3,'MarkerFaceColor',condColor,'MarkerEdgeColor',condColor)

% x=(1:1:1001)';
% up = meanQty+stdQty;
% down = meanQty-stdQty;
% patch([x fliplr(x)],[up fliplr(up)],color)
% patch([x fliplr(x)],[down fliplr(down)],color)
%plot(meanQty+stdQty,'--','LineWidth',1.5,'Color',color)
%plot(meanQty-stdQty,'--','LineWidth',1.5,'Color',color)


end