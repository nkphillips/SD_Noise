% Initialize fixation

%% Fixation space size 

p.fixation_space_deg = 1; % in visual degrees
p.fixation_space_px = round(p.fixation_space_deg * w.ppd); % in pixels
if ~mod(p.fixation_space_px, 2), p.fixation_space_px = p.fixation_space_px-1;end % force the fixation space to an odd # of pixels so that there's a middle

%% Fixation aperture
% The aperture is how much of the fixation space texture will actually go through (the apparent fixation space size).

p.fixation_space_aperture_deg = 0.8; % diameter in visual degrees
p.fixation_space_aperture_px = p.fixation_space_aperture_deg * w.ppd; % px
if ~mod(p.fixation_space_aperture_px,2), p.fixation_space_aperture_px=p.fixation_space_aperture_px-1;end % force to an odd # of pixels

% The fixation space is a small circular area around fixation that will be drawn over stimuli to help with fixating
fixation_space_bgcolor = w.gray;
fixation_space = (ones(p.fixation_space_px) * fixation_space_bgcolor); % Create a gray square

% Create an aperture for the fixation space texture
[x_cart, y_cart] = meshgrid(-(p.fixation_space_px/2):(p.fixation_space_px/2) -1, -(p.fixation_space_px/2):(p.fixation_space_px/2)-1); % create coordinates for texture: -x to x, -y to y;
fixation_space_eccen_radii = sqrt(x_cart.^2+y_cart.^2); % Calculating the eccentricity (radius) of each point in the meshgrid relative to the center of the 2D image

fixation_space_aperture = ones(p.fixation_space_px)*255; % The aperture is fully transparent to start (the entire sqaure fixation space is let through)
fixation_space_aperture(fixation_space_eccen_radii >= (p.fixation_space_aperture_px/2)) = 0; % Have the aperture block the fixation space texture outside of the specified radius
fixation_space_aperture = imgaussfilt(fixation_space_aperture, w.ppd*0.1); % Apply a guassian blur to the aperture to smooth its edge
fixation_space(:,:,2) = fixation_space_aperture; % Apply the aperture layer to the fixation space texture (the gray square)

%% Fixation dot

p.fixation_dot_deg = 0.15; % inner fixation dot size in visual degrees (default = 0.15°)
p.fixation_dot_px = round(p.fixation_dot_deg * w.ppd); % fixation size in pixels
if mod(p.fixation_dot_px, 2), p.fixation_dot_px = p.fixation_dot_px -1; end % force to an even # of pixels

%% Set colors for fixation

p.fixation_dot_color = w.black; %  Color of the fixation dot during the task blocks
p.rest_fixation_color = w.white; % Color of the fixation during rest periods