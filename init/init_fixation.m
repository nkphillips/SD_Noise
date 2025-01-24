%%% Initialize fixation

%{

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

%% Fixation space size 

p.fixation_space_deg = 1; % in visual degrees
p.fixation_space_px = round(p.fixation_space_deg * w.ppd); % in pixels
if ~mod(p.fixation_space_px, 2), p.fixation_space_px = p.fixation_space_px-1;end % force the fixation space to an odd # of pixels so that there's a middle

%% Fixation aperture
% The aperture is how much of the fixation space texture will actually go through (the apparent fixation space size).

p.fixation_space_aperture_deg = 0.32; % diameter in visual degreen
p.fixation_space_aperture_px = p.fixation_space_aperture_deg * w.ppd; % px
if ~mod(p.fixation_space_aperture_px,2), p.fixation_space_aperture_px=p.fixation_space_aperture_px-1;end % force to an odd # of pixels

%% Fixation dot

p.fixation_dot_deg = 0.15; % inner fixation dot size in visual degrees (default = 0.15°)
p.fixation_dot_px = round(p.fixation_dot_deg * w.ppd); % fixation size in pixels
if mod(p.fixation_dot_px, 2), p.fixation_dot_px = p.fixation_dot_px -1; end % force to an even # of pixels

%% Set colors for the cue at fixation

p.fixation_dot_color = black; %  Color of the fixation dot during the task blocks
p.rest_fixation_color = white; % Color of the fixation during rest periods