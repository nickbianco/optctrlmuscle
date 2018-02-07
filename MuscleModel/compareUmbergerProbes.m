
muscle = 'med_gas';
switch muscle
    case 'med_gas'
        max_isometric_force = 3115.51;
        optimal_fiber_length = 0.051;
        tendon_slack_length = 0.39872;
        pennation_angle_at_optimal = 0.16568;
        max_contraction_velocity = 10.0;
        rST = 0.566;
        rFT = 1 - rST;
end

import org.opensim.modeling.*
model = Model('isolated_active_muscle');
origin = Body('origin', 1.0, Vec3(