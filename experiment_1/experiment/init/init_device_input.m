% Initialize device input

% Enable unified mode of KbName, so KbName accepts identical key names on
% all operating systems:
KbName('UnifyKeyNames');

%% Define device and relevant keys

p.device_number = 0;
[Kb_indices, product_names, ~] = GetKeyboardIndices;

switch p.which_setup
    case 0 % If using Macbook
        p.device_string = 'Apple Internal Keyboard / Trackpad';
    case 1 % If using 3329B
        p.device_string = 'Dell Dell USB Keyboard';
    case 2 % If using 3329C
        p.device_string = 'LiteOn Lenovo Traditional USB Keyboard';
    case 3 % If using 3329D
        p.device_string = 'Dell Dell USB Keyboard';
    case 4 % If using S32D850
        p.device_string =  'USB-HID Keyboard';
end


% Define keypress numbers
p.keypress_numbers = [KbName('1!') KbName('2@')];

% Assign trigger key
p.trigger_key = KbName('Space');

%% Scan for device number

for i = 1:length(product_names)
    if strcmp(product_names{i}, p.device_string)
        p.device_number = Kb_indices(i);
        break;
    end
end

%% Error if no device is found

if p.device_number == 0
    error('No device by that name was detected');
end

%% Turn on keyboard input

KbQueueCreate(p.device_number);
KbQueueStart(p.device_number);
