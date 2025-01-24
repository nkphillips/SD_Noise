%%% Initialize display parameters
% Ideal distance: 1 cm equals 1 visual degree at 57 cm

%{

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

%% Grab the available screens

screens = Screen('Screens'); 

%% Set screen parameters

p.display_setup = 'Macbook';

w.use_screen = min(screens); % If there are two or more displays, 'max' should grab the most external display.
w.refresh_rate = 60; % display refresh rate in Hertz, Hz
w.view_distance = 57; % in centimeters, cm
w.screen_width = 30; % cm
w.screen_width_px = 1512; % in pixels, px
w.screen_height_px = 982; % px

%% Calculate visual angle of the display, pixels per degree of visual angle, size of a pixel, and Nyquist frequency

screen_length = w.screen_width;
screen_length_px = w.screen_width_px;

w.visual_angle = 2 * atan2d(screen_length/2,  w.view_distance); % Visual angle of the whole screen in degrees
w.ppd = floor(screen_length_px/w.visual_angle); % Pixels per degree of visual angle
w.px_size = screen_length/screen_length_px; % size of pixel in cm
w.f_Nyquist = 1/(2*w.px_size);  

%% Define the background color

gray = 127; % halfway point between 0–254 (255 elements total)
p.bg_color = gray;

%% Define other colors

white = [255 255 255];
black = [0 0 0];
red = [255 0 0];
green = [0 255 0];
blue =  [0 0 255];
