%%% simulate_trials

clear all; close all; clc;

%% Set directories

script_dir = pwd;
functions_dir = '../functions'; addpath(functions_dir);

%% Define parameters

% Define contrasts
stimuli.contrast = [0.9 0.5 0.25]; % high contrast, medium, low; unit: % Michelson contrast

% Define orientation bandpass filter widths
stimuli.bp_filter_width = [0.1 5 20]; % low noise, medium, high ; unit: °

% Check that the number of levels between stimulus contrast and bp filter widths match
if length(stimuli.bp_filter_width) == length(stimuli.contrast)
    p.num_levels = length(stimuli.contrast);
else
    disp('Condition levels do not match in length!');
end

% Define orientations
stimuli.orientation_min = 0;
stimuli.orientation_max = 179;

%% Define number and sequence of events

p.cond_names = {'contrast', 'filter'};
p.num_conds = numel(p.cond_names);

p.num_blocks = 6;
while mod(p.num_blocks, p.num_conds) ~= 0, p.num_blocks = input(['Error! Number of blocks must be a multiple of ' num2str(p.num_conds) ': ']); end

p.num_blocks_per_cond = p.num_blocks / p.num_conds;

p.block_order = repmat(1:p.num_conds, 1, p.num_blocks_per_cond);
p.block_order = Shuffle(p.block_order);

first_block = nan(1, p.num_conds);
for cond = 1:p.num_conds
    first_block(cond) = find(p.block_order == cond, 1, 'first');
end

p.num_trials_per_unique_cond = 3000; % default = 35
while mod(p.num_trials_per_unique_cond, p.num_levels) ~= 0, p.num_trials_per_unique_cond = input(['Error! Number of trials must be a multiple of ' num2str(p.num_levels) ': ']); end

p.num_trials_per_block = (p.num_levels * p.num_trials_per_unique_cond)/p.num_blocks_per_cond;
p.trial_events = nan(p.num_trials_per_block, 3, p.num_blocks); % num_trials x [test_orientation, probe_orientation, cond_lvl] x num_blocks 
p.correct_response = nan(p.num_trials_per_block, p.num_blocks);

%% Generate level order, orientations, correct response

for n_block = 1:p.num_blocks
   
    curr_cond = p.block_order(n_block);

    if p.block_order(n_block) == 1
        
        level_order = BalanceFactors(p.num_trials_per_block/p.num_levels, 1, 1:length(stimuli.contrast));
        
    elseif p.block_order(n_block) == 2
        
        level_order = BalanceFactors(p.num_trials_per_block/p.num_levels, 1, 1:length(stimuli.bp_filter_width));
        
    end

    % Sample Test orientations
    test_orientation = round(stimuli.orientation_min + (stimuli.orientation_max - stimuli.orientation_min) .* rand(p.num_trials_per_block, 1));

    % Pre-allocate the probe orientations
    probe_orientation = calc_probe_orientation(test_orientation, 5);
    corrected_probe_orientation = correct_orientation(probe_orientation);

    % Storing trial events
    p.trial_events(:,:,n_block) = [test_orientation, corrected_probe_orientation, level_order];
    test_orientation_col = 1;
    probe_orientation_col = 2;
    level_order_col = 3;

    % Storing correct response
    cclockwise_trials = double(probe_orientation < test_orientation);
    cclockwise_trials(cclockwise_trials == 0) = 2;
    p.correct_response(:,n_block) = cclockwise_trials;

end

p.num_trials = p.num_trials_per_block * p.num_blocks;

%% Check condition distribution

% "scoreboard" for the # of level pairs
pair_count = zeros(p.num_levels, p.num_levels, p.num_conds);

% store the indices of the current trial if the current and previous trial
% match the current level pair
% note: you must use a cell array
pair_indices = cell(p.num_levels, p.num_levels, p.num_blocks_per_cond, p.num_conds);

delta_theta = cell(p.num_levels, p.num_levels, p.num_conds);

[subplotX, subplotY] = get_subplot_dimensions(p.num_levels * p.num_levels);

for cond = 1:p.num_conds
    
    curr_blocks = find(p.block_order == cond);
    
    curr_lvls = squeeze(p.trial_events(:,end,curr_blocks));
    curr_probe_orientations = squeeze(p.trial_events(:,2,curr_blocks));
    
    % Count # of unique trial pairs (ignore 1st trial)
    plot_counter = 1;
    for a = 1:p.num_levels % curr
        for b = 1:p.num_levels % prev

            % use a and b to identify the current trial pair (trial_n and trial_{n-1})
            % count the number of those pairs across each block
            for n_block = 1:length(curr_blocks)
                for n_trial = 2:size(curr_lvls,1)
                    
                    prev_trial_lvl = curr_lvls(n_trial-1, n_block);
                    curr_trial_lvl = curr_lvls(n_trial, n_block);
                    
                    if prev_trial_lvl == a  && curr_trial_lvl == b
                        
                        pair_count(b, a, cond) = pair_count(b, a, cond) + 1;
                        
                        % store the index of the current trial
                        pair_indices{b, a, n_block, cond} = [pair_indices{b, a, n_block, cond}, n_trial];
                   
                        % Calculate delta theta and store it into delta_theta
                        delta_theta{b, a, cond} = [delta_theta{b, a, cond}, ...
                            curr_probe_orientations(n_trial-1, n_block) - curr_probe_orientations(n_trial, n_block)];   

                    end
                end
            end  


        end
    end
    
end

% for cond = 1:p.num_conds
%     figure('Name',  p.cond_names{cond}, 'Color', [1 1 1])
%     plot_counter = 1;
%     for a = 1:p.num_levels
%         for b = 1:p.num_levels

%             % Plot histogram
%             subplot(subplotX, subplotY, plot_counter);
%             histogram(delta_theta{b, a, cond});
%             box off;
%             xlim([-180 180]);
%             % ylim([0 max(pair_count(:))])    

%             plot_counter = plot_counter + 1;

%         end
%     end
% end

% test for diff in number of pairs between condition

% x = pair_count(:,:,1);
% y = pair_count(:,:,2);
% 
% [h, p_value] = ttest(x(:), y(:));
% 
% figure, histogram(x(:)), hold on;  histogram(y(:))

%% Simulate subject 

% Define parameters
amplitude = 4.5;
sigma = 10;
width = 10;
noise = 1;

% Pre-allocate output
y = cell(p.num_levels, p.num_levels, p.num_conds);

for cond = 1:p.num_conds
    figure('Name',  p.cond_names{cond}, 'Color', [1 1 1])
    plot_counter = 1;
    for a = 1:p.num_levels
        for b = 1:p.num_levels

            % Generate values
            params = [amplitude, width];
            y{b, a, cond} = gaussian_prime(params, delta_theta{b, a, cond}) + noise .* randn(1,length(delta_theta{b, a, cond}));

            % Plot 
            subplot(subplotX, subplotY, plot_counter);
            scatter(delta_theta{b, a, cond}, y{b, a, cond}, 20, 'Marker','o' , 'MarkerFaceColor',[0 0 0], 'MarkerEdgeColor',[1 1 1]);
            xlim([-180 180]);

            hold on
            line([min(xlim), max(xlim)], [0 0],'LineStyle','-','Color', [0 0 0])
            box off;
            xlabel('\Delta\theta')
            ylabel('Bias (°)')
            ylim([-10 10])    
            set(gca, 'TickDir','out')

            plot_counter = plot_counter + 1;
            
        end
    end
end

