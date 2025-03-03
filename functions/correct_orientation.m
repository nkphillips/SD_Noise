%%% correct_orientation

% take an orientation for a line
% correct it so that the angle falls along the desired quadrant in "compass axis" (ie vertical = 0°)
% assume orientation range [0, 180)

function corrected_orientation = correct_orientation(original_orientation)

    % Transform probe orientation to match the compass axis
    corrected_orientation = original_orientation + 90;
    corrected_orientation(corrected_orientation >= 180) = corrected_orientation(corrected_orientation >= 180) - 90;

end