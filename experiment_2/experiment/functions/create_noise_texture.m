%% Create white noise texture and aperture

function noise_texture = create_noise_texture(texture_height, texture_width)

    if nargin > 1
    
        noise_texture = 2 * rand(texture_height, texture_width) - 1;
    
    else

        noise_texture = 2 * rand(texture_height) - 1;
        
    end

end
