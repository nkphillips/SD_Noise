%%% init_lines
% Using Screen('DrawLines')
% "xy" is a two-row vector containing the x and y coordinates of the line
% segments: Pairs of consecutive columns define (x,y) positions of the starts and
% ends of line segments. All positions are relative to "center" (default center is
% [0 0]).

cue_length = p.fixation_space_aperture_px;
cue_thickness = w.ppd * 0.1;
cue_color = w.black;
cue_horizontal_offset = p.fixation_dot_px + w.ppd/2;

% Define coordinates relative to [0, 0]
left_cue_start = w.centerX - (cue_length + cue_horizontal_offset);
left_cue_end = w.centerX - cue_horizontal_offset;

right_cue_start = w.centerX + cue_horizontal_offset;
right_cue_end = w.centerX + cue_horizontal_offset + cue_length;

y_position = w.centerY;  % Same y-coordinate for both lines to keep them horizontal

% Define the xy matrix for Screen('DrawLines')
cues.xy = [left_cue_start, left_cue_end, right_cue_start, right_cue_end;
      y_position, y_position, y_position, y_position];

% Clamp line widths to range supported by graphics hardware:
[minsmooth,maxsmooth] = Screen('DrawLines', w.window);
cue_thickness = min(max(cue_thickness, minsmooth), maxsmooth);