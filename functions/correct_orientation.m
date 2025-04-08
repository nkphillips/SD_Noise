%%% correct_orientation

% Take an orientation for a line to use in Screen('DrawLine')
% Correct it so that the angle falls along the desired quadrant in "compass axis" (ie vertical = 0°)
% Assume orientation range [0, 180)

function corrected_orientation = correct_orientation(original_orientation)

    corrected_orientation = original_orientation + 90;

end