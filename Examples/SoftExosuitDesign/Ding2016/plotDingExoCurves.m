DingExoCurves = load('DingExoCurves.mat');

time = DingExoCurves.time;

cond = {'esep','eslp','lsep','lslp'};
momArms = zeros(length(time),length(cond));
figure;
for c = 1:length(cond)
    momArms(:,c) = DingExoCurves.(cond{c}).r;
    plot(momArms(:,c))
    hold on
end

momArmsMean = mean(momArms,2);
momArmsStd = std(momArms,0,2);

plot(momArmsMean,'k-','LineWidth',2)
plot(momArmsMean+momArmsStd,'k--','LineWidth',1.5)
plot(momArmsMean-momArmsStd,'k--','LineWidth',1.5)