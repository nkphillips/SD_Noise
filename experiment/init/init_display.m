% Initialize display parameters
% Ideal distance: 1 cm equals 1 visual degree at 57 cm

screens = Screen('Screens'); % Grab the available screens

%% Set screen parameters

if p.which_setup == 0 % Macbook

    p.display_setup = 'Macbook';

    w.use_screen = min(screens); % If there are two or more displays, 'max' should grab the most external display.
    w.refresh_rate = 60; % display refresh rate in Hertz, Hz
    w.view_distance = 57; % in centimeters, cm
    w.screen_width = 30; % cm
    w.screen_width_px = 1512; % in pixels, px
    w.screen_height_px = 982; % px

elseif p.which_setup == 1 % 3329C

    p.display_setup = '3329C_ASUS';

    w.use_screen = 0;
    w.gamma_correct = 0;
    w.refresh_rate = 85; % display refresh rate in Hertz, Hz
    w.view_distance = 40;  % cm
    w.screen_width = 40.7; %  cm
    w.screen_width_px = 2560; % px
    w.screen_height_px = 1440; % px

elseif p.which_setup == 2 % S32D850

    p.display_setup = 'S32D850';

    w.use_screen = min(screens); % If there are two or more displays, 'max' should grab the most external display.
    w.refresh_rate = 60; % display refresh rate in Hertz, Hz
    w.view_distance = 57; % in centimeters, cm
    w.screen_width = 70.8; % cm
    w.screen_width_px = 2560; % in pixels, px
    w.screen_height_px = 1440; % px


end

if p.half_screen

    switch p.which_setup
        case 0 % Macbook
            w.view_distance = 57; % in centimeters, cm
            w.screen_width = 42; % cm
            w.screen_width_px = w.screen_width_px/2; % in pixels, px
            w.screen_height_px = w.screen_height_px/2; % px
        case 1 % 3329C
            w.view_distance = 40; % in centimeters, cm
            w.screen_width = 40.7; %  cm
            w.screen_width_px = 1280; % px
            w.screen_height_px = 720; % px
        case 2 % S32D850
            w.view_distance = 57; % in centimeters, cm
            w.screen_width = 70.8; % cm
            w.screen_width_px = w.screen_width_px/2; % in pixels, px
            w.screen_height_px = w.screen_height_px/2; % px
    end

end

%% Calculate visual angle of the display, pixels per degree of visual angle, size of a pixel, and Nyquist frequency

screen_length = w.screen_width;
screen_length_px = w.screen_width_px;

w.visual_angle = 2 * atan2d(screen_length/2,  w.view_distance); % Visual angle of the whole screen in degrees
w.ppd = floor(screen_length_px/w.visual_angle); % Pixels per degree of visual angle
w.px_size = screen_length/screen_length_px; % size of pixel in cm
w.f_Nyquist = 1/(2*w.px_size);

%% Define the background color

w.gray = 127; % halfway point between 0–254 (255 elements total)
w.bg_color = w.gray;

%% Define other colors

w.white = [255 255 255];
w.black = [0 0 0];
w.red = [255 0 0];
w.green = [0 255 0];
w.blue =  [0 0 255];
