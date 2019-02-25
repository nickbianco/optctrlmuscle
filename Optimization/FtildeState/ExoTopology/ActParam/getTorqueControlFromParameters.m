function [control] = getTorqueControlFromParameters(peakTorque, peakTime, riseTime, fallTime, numColPoints)

% set spline nodes based on parameters
node1_x = peakTime - riseTime;
node2_x = peakTime;
node3_x = peakTime + fallTime;

node1_y = 0.01*peakTorque;
node2_y = peakTorque;
node3_y = 0;

% first segment: linear
pre_slope = node1_y/node1_x;
if isinf(pre_slope) || (pre_slope < 0) || isnan(pre_slope) || (pre_slope > 1)
   pre_slope = 0; 
end
x = linspace(0, node1_x, round(node1_x*numColPoints));
pre_rise = x(1:end-1)*pre_slope;

% rising torque: cubic spline
X = [node1_x node2_x];
Y = [node1_y node2_y];
XX = linspace(node1_x, node2_x, round((node2_x-node1_x)*numColPoints));
rising_torque = spline(X,[pre_slope Y 0],XX);

% falling torque: cubic spline
X = [node2_x node3_x];
Y = [node2_y node3_y];
XX = linspace(node2_x, node3_x, round((node3_x-node2_x)*numColPoints));
post_slope = -peakTorque/fallTime;
if isinf(post_slope) || (post_slope < 0) || isnan(post_slope) || (post_slope > 1)
   post_slope = 0; 
end
falling_torque = spline(X,[0 Y post_slope],XX);

% last segment: zeros
post_fall = zeros(1, round((1-node3_x)*numColPoints));

% concatenate segments
v = [pre_rise rising_torque falling_torque post_fall];

% re-interpolate to ensure the correct number of collocation points
x = linspace(0, 1, length(v));
xq = linspace(0, 1, numColPoints);
control = interp1(x, v, xq)';

end