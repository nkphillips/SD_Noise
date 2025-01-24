%%% gen_unique_seq
% Given a matrix size, seq_size
% Populate each row with elements from values
% These values do not repeat within a given window_size

function seq = gen_unique_seq(seq_size, values, window_size)

    % Initialize an empty array to hold the sequence
    seq = nan(seq_size);
    
    for i = 1:size(seq,1)
        for j = 1:size(seq,2)
    
            % Start with all candidate values
            candidate_values = values;
    
            % Determine the range of past values to avoid within window_size
            if j > 1
    
                % Set the starting index for previous values to avoid
                start_idx = max(1, j - window_size);
    
                % Gather all prior values within the window to remove
                recent_values = seq(i, start_idx:j-1);
    
                % Remove these values from candidates
                candidate_values(ismember(candidate_values, recent_values)) = [];
    
            end
    
            % Randomly select a candidate from remaining values
            seq(i,j) = datasample(candidate_values,1);
    
        end
    end

end
