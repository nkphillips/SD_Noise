%% Prepare Psychtoolbox for Imaging 

% PsychImaging('PrepareConfiguration');
% PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');

%% Open Window

p.demo_run = 0;

if p.demo_run

    [window, p.screen_rect_sz_px] = PsychImaging('OpenWindow', w.use_screen, p.bg_color, [0 0 800 600]); % smaller screen

else
    
    [window, p.screen_rect_sz_px] = PsychImaging('OpenWindow', w.use_screen, p.bg_color, [0 0 w.screen_width_px w.screen_height_px]);

end

commandwindow; 
% HideCursor;

%% % Grab screen refresh rate 

t.monitor_refresh_dur = Screen('GetFlipInterval', window); % in seconds

%% Enable alpha blending

Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%% Load CLUT

w.DefaultCLUT = Screen('ReadNormalizedGammaTable', window);

%% Define center coordinates

w.centerX = p.screen_rect_sz_px(3)/2; 
w.centerY = p.screen_rect_sz_px(4)/2; 

%% Define patches

% Define fixation space and dot patches (how and where the fixation space and dot will be drawn)
fixation_dot_patch = CenterRectOnPoint([0 0 p.fixation_dot_px p.fixation_dot_px], w.centerX, w.centerY);
fixation_space_patch = CenterRectOnPoint([0 0 p.fixation_space_px p.fixation_space_px], w.centerX, w.centerY);

% Define stimulus patch
stimulus_patch = CenterRectOnPoint([0 0 stimuli.diameter_px stimuli.diameter_px], w.centerX, w.centerY);

%% Text settings

Screen('TextStyle', window, 1); % 0=normal, 1=bold, 2=italic
Screen('TextSize', window, 18);

%% Loading screen

loading_text = 'Loading stimuli...'; % standby text for loading experiment and stimuli
loading_text_boundary = Screen('TextBounds', window, loading_text);
loading_text_patch = CenterRectOnPoint(loading_text_boundary, w.centerX, w.centerY);

Screen('DrawText', window, loading_text, loading_text_patch(1),  loading_text_patch(2), white);
Screen('Flip', window);

%% Make fixation space

% The fixation space is a small circular area around fixation that will be drawn over stimuli to help with fixating
fixation_space_bgcolor = gray;
fixation_space = (ones(p.fixation_space_px) * fixation_space_bgcolor); % Create a gray square

% Create an aperture for the fixation space texture
[x_cart, y_cart] = meshgrid(-(p.fixation_space_px/2):(p.fixation_space_px/2) -1, -(p.fixation_space_px/2):(p.fixation_space_px/2)-1); % create coordinates for texture: -x to x, -y to y;
fixation_space_eccen_radii = sqrt(x_cart.^2+y_cart.^2); % Calculating the eccentricity (radius) of each point in the meshgrid relative to the center of the 2D image

fixation_space_aperture = ones(p.fixation_space_px)*255; % The aperture is fully transparent to start (the entire sqaure fixation space is let through)
fixation_space_aperture(fixation_space_eccen_radii >= (p.fixation_space_aperture_px/2)) = 0; % Have the aperture block the fixation space texture outside of the specified radius
fixation_space_aperture = imgaussfilt(fixation_space_aperture, 1); % Apply a guassian blur to the aperture to smooth its edge
fixation_space(:,:,2) = fixation_space_aperture; % Apply the aperture layer to the fixation space texture (the gray square)

created_fixation_space = Screen('MakeTexture', window, fixation_space);
