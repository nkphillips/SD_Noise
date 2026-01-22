%%% Create circular aperture

%{
========================================================================
Make an aperture that is effectively a circular mask for your stimuli.
Everything outside of this circular mask remains gray when drawn.
========================================================================

A. Create a circular mask
    
    1. Create a cartesian grid equal to the size of your stimulus using meshgrid(). 
       Every x,y coordinate can be thought of as a pixel coordinate.

    2. Calculate the radius for every pixel location in the grid. 

    3. Initialize a matrix of 1s, size equal to the aperture layer (which itself is already be equal to the size of the stimulus behind the aperture). 
       Values inside the desired radius are then changed to 0. This is your circular mask, before being multiplied by an alpha value!

    4. Smooth the edges of the mask with a gaussian filter (recommended).

B. Apply the desired alpha value to the mask.

   0: completely transparent (the texture containing the alpha layer is completely invisible)
   255: completely opaque (the texture containing the alpha layer is completely visible) 

   Values inside the circular mask were already forced to 0, remaining completely transparent after this operation. 

C. Store the mask in the second layer of the gray patch.

   The end result is a gray patch with an alpha layer that is a transparent circlular mask, allowing the stimulus behind the gray patch to be visible only within the circular mask.
   (i.e., the gray patch drawn over the stimulus is visible only outside a given radius)

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

%% Create aperture mask

% Create cartesian coordinates for texture: -x to x, -y to y
[x,y] = meshgrid(-stimuli.aperture_diameter_px/2:stimuli.aperture_diameter_px/2-1,-stimuli.aperture_diameter_px/2:stimuli.aperture_diameter_px/2-1);

% Calculate eccentricity (e.g., radius) of each pixel coordinate in the cartesian grid, from the center of the grid
r = sqrt((x).^2 + (y).^2); % solving for r in: x^2 + y^2 = r^2

% Create a matrix of 1s
mask = ones(stimuli.aperture_diameter_px);

% Set values inside desired eccentricity to be fully transparent
mask(r <= stimuli.mask_radius_px) = 0;

% Create gaussian filter and apply to the aperture to smooth its edges %
gauss_filter_sd = 0.1 * w.ppd; % smooth with respect to a standard deviation, sigma, that is a fraction of the ppd
mask = imgaussfilt(mask, gauss_filter_sd); % apply gaussian filter to the mask

%% Store aperture with alpha layer

stimuli.aperture(:,:,1) = ones(size(mask)) * 127; % gray background
stimuli.aperture(:,:,2) = mask * stimuli.aperture_alpha; % circular transparent mask

%% Visualize aperture
%{

figure

imshow(stimuli.aperture)

%}

%% Make aperture

% Aperture/mask
aperture_made = Screen('MakeTexture', window, stimuli.aperture);
