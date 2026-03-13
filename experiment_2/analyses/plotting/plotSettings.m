% plotSettings
%   Returns a structure of plot settings

function ps = plotSettings()

%% Font settings

ps.font_type = 'Helvetica';
ps.axes_label_font_size = 14;
ps.axes_tick_font_size = 13;
ps.tick_length = 0.020;
ps.line_width = 1;

%% Color settings

ps.colors.red = [204 0 0]/255;
ps.colors.green = [0 153 0]/255;
ps.colors.blue = [0 76 152]/255;
ps.colors.purple = [102 51 204]/255;
ps.colors.yellow = [0.9, 0.8, 0];
ps.colors.orange = [0.8, 0.5, 0.1];
ps.colors.black = [0 0 0];
ps.colors.white = [1 1 1];
ps.colors.gray = ps.colors.white/2;

end