function orientation = sample_orientation(min_orientation, max_orientation, num_samples)

    % orientation = round(min_orientation + (max_orientation - min_orientation) .* rand(num_samples, 1));
    orientation = randi([min_orientation, max_orientation], num_samples, 1);   

end