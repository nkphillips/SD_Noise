%%% open_window
% Opens and preps window for stimulus presentation

% Written by Luis D. Ramirez
% UCSD
% lur003@ucsd.edu

%% Prepare Psychtoolbox for Imaging 

% PsychImaging('PrepareConfiguration');
% PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');

%% Open Window

[w.window, p.screen_rect_sz_px] = PsychImaging('OpenWindow', w.use_screen, w.bg_color, [0 0 w.screen_width_px w.screen_height_px]);

commandwindow; 

if ~p.half_screen
    HideCursor;
end

%% Enable alpha blending

Screen('BlendFunction', w.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%% Load CLUT

w.DefaultCLUT = Screen('ReadNormalizedGammaTable', w.window);

if p.which_setup == 1 && w.gamma_correct % 3329D
    
    LoadCLUT
    
elseif p.which_setup == 2 % If @ the scanner
    
    Screen('LoadNormalizedGammaTable', w.window, linspace(0,1, 256)'*ones(1,3));
    
end

%% Define center coordinates

w.centerX = p.screen_rect_sz_px(3)/2; 
w.centerY = p.screen_rect_sz_px(4)/2; 

%% Text settings

Screen('TextStyle', w.window, 1); % 0=normal, 1=bold, 2=italic
Screen('TextSize', w.window, 18);

