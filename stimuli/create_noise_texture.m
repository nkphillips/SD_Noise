%% Create white noise texture and aperture

function noise_texture = create_noise_texture(stimuli)
    noise_texture = 2 * rand(stimuli.diameter_px, stimuli.diameter_px) - 1; % generates random noise
end
