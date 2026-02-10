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

function aperture = create_circular_aperture(width, height, radius)

%% Create aperture mask

% Create cartesian coordinates for texture: -x to x, -y to y
[x,y] = meshgrid(-width/2:width/2-1,-height/2:height/2-1);

% Calculate eccentricity (e.g., radius) of each pixel coordinate in the cartesian grid, from the center of the grid
r = sqrt((x).^2 + (y).^2); % solving for r in: x^2 + y^2 = r^2

% Create a matrix of 1s
aperture = ones(height,width);

% Set values inside desired eccentricity to be fully transparent
aperture(r <= radius) = 0;

end