%%% create_wedge_aperture

%{

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

function wedge = create_wedge(ang, width, height, ecc)

%% Create coordinate space for image 
% -x to x, -y to y

[x_cart, y_cart] = meshgrid(-(width/2):(width/2)-1, -(height/2):(height/2)-1); 

%% Get the polar angle of every x,y coordinate

[theta, radii] = cart2pol(x_cart, y_cart);

angles = rad2deg(theta);

%% Eliminate coordinates that are outside the wedge

% Start with an aperture that lets everything through (in this case the gray background)
wedge = ones(height, width); 

% Set x,y coordinates with polar angles outside the desired angle to 0 (this will knock out the gray and allow the texture underneath through)
wedge(abs(angles) >= ang) = 0;

% Set coordinates outside the radius of the wedge to 1
wedge(radii > ecc) = 1;

end