%%% Initialize device input
% Enable unified mode of KbName, so KbName accepts identical key names on
% all operating systems:

%{

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

%%

KbName('UnifyKeyNames');

%% Define device and relevant keys

device_number = 0;
[Kb_indices, product_names, ~] = GetKeyboardIndices;

device_string = 'Apple Internal Keyboard / Trackpad'; % Macbook
%     device_string =  'USB-HID Keyboard'; % external keyboard

% Define keypress numbers
keypress_numbers = [KbName('LeftArrow') KbName('RightArrow') KbName('DownArrow') KbName('Escape')];


%% Scan for device number

for i = 1:length(product_names)
    if strcmp(product_names{i}, device_string)
        device_number = Kb_indices(i);
        break;
    end
end

%% Error if no device is found

if device_number == 0
    error('No device by that name was detected');
end