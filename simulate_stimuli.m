%{

Simulates orientation bandpass-filtered noise adjustment via left/right keys

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

%% Clean and prepare workspace

clc;
clear all; %#ok<CLALL> 
close all

commandwindow; % force cursor to command window
Priority(1); % Set MATLAB/Psychtoolbox to "high" priority level

% Grab date
t.the_date = datestr(now, 'yyyymmdd'); % Grab today's date
t.the_time = datestr(now,'HHMM'); % Grab current time

% Generate unique seed for random number generator (rng)
t.my_rng_seed = sum(100*clock);
rng(t.my_rng_seed);

%% Set directories

script_dir = pwd;
init_dir = 'init'; addpath(init_dir);
functions_dir = 'functions'; addpath(functions_dir);
stim_dir = 'stimuli'; addpath(stim_dir);

%% Toggles

% Sync Test
skip_sync_test = 1;

if skip_sync_test
    Screen('Preference', 'SkipSyncTests', 1); % Set to 1 if running on macOS
else
    Screen('Preference', 'SkipSyncTests', 0);
end

init_device_input % defines the device and keypresses to listen to
init_display % creates 'w' struct with fields related to the display setup
init_fixation % set fixation parameters
init_stimuli % creates stimuli struct

%% Open window

open_window

%% Make circular aperture

create_circular_aperture

%% Create orientation bandpass filtered noise

make_filtered_noise

%% Create grating

make_grating

%% Initialize response params

% Set initial probe orientation (this will be random)
curr_orientation_indx = randi(p.num_orientations);
curr_orientation = p.orientations(curr_orientation_indx);

% Fixation color is black by default, but will turn red if orientation is at the limit
fixation_color = black; 

%% Choose stimulus to draw

stim = 'grating';

if strcmp(stim,'grating')

    stim_made = grating_made;

elseif strcmp(stim,'noise')

    stim_made =  noise_made;

end

%% Step into response loop

while 1

    % Draw probe
    Screen('DrawTexture', window, stim_made(curr_orientation_indx), [], stimulus_patch);

    % Draw aperture
    Screen('DrawTexture', window, aperture_made,[], stimulus_patch);

    % Draw fixation
    Screen('DrawTexture', window, created_fixation_space, [], fixation_space_patch); % Fixation circle
    Screen('FillOval', window, fixation_color, fixation_dot_patch); % Fixation dot

    Screen('Flip', window);

    % Check for response
    [key_pressed, ~, first_press] = KbCheck(device_number);
    which_press = find(first_press);

    if key_pressed

        [~, curr_orientation_indx] = find(p.orientations == curr_orientation);

        if which_press(1) == keypress_numbers(1) % Left arrow, decrease orientation

            new_orientation_indx  = curr_orientation_indx - 1;

            % Check if orientation index is within lower bound
            if new_orientation_indx  <= 1
                new_orientation_indx  = 1;
                fixation_color = red;
            else
                fixation_color = black;
            end

            curr_orientation = p.orientations(new_orientation_indx );

            disp(['orientation = ' num2str(curr_orientation) ' cpd'])

        elseif which_press(1) == keypress_numbers(2) % Right arrow, increase orientation

            new_orientation_indx  = curr_orientation_indx + 1;

            % Check if orientation index is within upper bound
             if new_orientation_indx  > length(p.orientations)
                new_orientation_indx  = length(p.orientations);
                fixation_color = red;
             else 
                fixation_color = black;
            end

            curr_orientation = p.orientations(new_orientation_indx );

            disp(['orientation = ' num2str(curr_orientation) ' cpd'])

        elseif which_press(1) == keypress_numbers(3)

            if strcmp(stim,'grating')

                stim = 'noise';
                stim_made = noise_made;

            elseif strcmp(stim,'noise')
                
                stim = 'grating';
                stim_made = grating_made;

            end

        elseif which_press(1) == keypress_numbers(end)
            
            sca;
            disp('User exited.')
            break;

        end

        % Check for response
        [key_pressed, ~, first_press] = KbCheck(device_number);
        which_press = find(first_press);

    end

end