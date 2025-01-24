%%% make_grating

%{

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

tic
%% Fixed parameters

sinusoid_orientation = 90;
sinusoid_phase = 0;

%% Make a grating for every SF

grating_made = nan(p.num_sf,1);

for n = 1:p.num_sf

    %% Get SF

    sinusoid_sf = p.sf(n); % spatial frequency (cycles/deg)
    sinusoid_sf = stimuli.diameter_px/w.ppd * sinusoid_sf; % scale sf to stimulus size

    %% Generate grating

    [x, y] = meshgrid(0:stimuli.diameter_px-1, 0:stimuli.diameter_px-1); % Generate pixel coordinate space
    
    init_grating =  sin(sinusoid_sf*2*pi / stimuli.diameter_px*(x.*sin(sinusoid_orientation*(pi/180)) + y.*cos(sinusoid_orientation*(pi/180)))-sinusoid_phase);
    
    filtered_grating = init_grating * stimuli.contrast * gray + gray;

    %{
    figure('Name','Filtered Grating','Color', [1 1 1]), imshow(filtered_grating)
    %}

    %% Make grating
    
    grating_made(n) = Screen('MakeTexture',window,filtered_grating);

end

disp(['Grating generation took ~' num2str(round(toc)) ' (s)'])
