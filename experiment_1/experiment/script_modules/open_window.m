%%% open_window
% Opens and preps window for stimulus presentation

% Written by Luis D. Ramirez
% UCSD
% lur003@ucsd.edu

%% Open Window

w.window = PsychImaging('OpenWindow', w.use_screen, w.bg_color, [0 0 w.screen_width_px w.screen_height_px]);

if ~p.half_screen
    HideCursor;
    commandwindow;
end

%% Enable alpha blending

Screen('BlendFunction', w.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%% Load CLUT

w.DefaultCLUT = Screen('ReadNormalizedGammaTable', w.window);

if any(p.which_setup == [1 2 3]) && w.gamma_correct

    load([dirs.monitor_cal_dir '/corrected_gamma_table_' p.display_setup '.mat'])
    w.CorrectedCLUT = repmat(corrected_gamma.chosen_table,1,3) * 255;
    Screen('LoadCLUT', w.window, w.CorrectedCLUT);

end

%% Define center coordinates

w.centerX = w.screen_width_px/2;
w.centerY = w.screen_height_px/2;

%% Text settings

Screen('TextStyle', w.window, 1); % 0=normal, 1=bold, 2=italic
Screen('TextSize', w.window, 18);

