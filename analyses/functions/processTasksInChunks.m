function [all_results, task_indices] = processTasksInChunks(task_list, num_chunks, use_parallel, processing_function, toggles, varargin)
    % processTasksInChunks - Process tasks in chunks with optional parallelization
    %
    % This function handles the entire workflow of chunking tasks, processing them
    % (either in parallel or sequentially), and compiling the results.
    %
    % Inputs:
    %   task_list - cell array of tasks to process
    %   num_chunks - number of chunks to create
    %   use_parallel - boolean, whether to use parfor (true) or for (false)
    %   processing_function - function handle to process each task
    %   varargin - additional arguments to pass to processing_function
    %
    % Outputs:
    %   all_results - cell array of results from all tasks
    %   task_indices - matrix of task indices for reconstruction
    
    num_tasks = length(task_list);
    
    % Pre-allocate output arrays
    all_results = cell(num_tasks, 1);
    
    % Calculate chunk size and number of chunks
    chunk_size = ceil(num_tasks / num_chunks);
    actual_num_chunks = ceil(num_tasks / chunk_size);
    
    % Display chunking information
    if use_parallel && toggles.disp_on
        disp(['    - Chunking ' num2str(num_tasks) ' tasks into ' num2str(actual_num_chunks) ' chunks (~' num2str(chunk_size) ' tasks per chunk)']);
    end
    
    % Process chunks
    if use_parallel 
        if toggles.disp_on, disp('    - Starting parallel processing...'); end
        parfor i_chunk = 1:actual_num_chunks
            chunk_start = (i_chunk - 1) * chunk_size + 1;
            chunk_end = min(i_chunk * chunk_size, num_tasks);
            chunk_tasks = task_list(chunk_start:chunk_end);
            
            % Pre-allocate temporary arrays for this chunk
            chunk_results = cell(chunk_end - chunk_start + 1, 1);
            
            % Process each task in the chunk
            for i_task = 1:length(chunk_tasks)
                
                % Call the processing function with the task and additional arguments
                if isempty(varargin)
                    result = processing_function(chunk_tasks{i_task});
                else
                    result = processing_function(chunk_tasks{i_task}, varargin{:});
                end
                
                % Store result in temporary array
                chunk_results{i_task} = result;
            end
            
            % Return chunk results as a struct
            chunk_output(i_chunk).results = chunk_results;
            chunk_output(i_chunk).start_idx = chunk_start;
            chunk_output(i_chunk).end_idx = chunk_end;
        end
        
        % Combine results from all chunks
        for ic = 1:actual_num_chunks
            all_results(chunk_output(ic).start_idx:chunk_output(ic).end_idx) = chunk_output(ic).results;
        end
        
    else
        if toggles.disp_on, disp('    - Starting sequential processing...'); end
        % Sequential processing
        for i_chunk = 1:actual_num_chunks
            chunk_start = (i_chunk - 1) * chunk_size + 1;
            chunk_end = min(i_chunk * chunk_size, num_tasks);
            chunk_tasks = task_list(chunk_start:chunk_end);
            
            % Process each task in the chunk
            for i_task = 1:length(chunk_tasks)
                
                % Call the processing function with the task and additional arguments
                if isempty(varargin)
                    result = processing_function(chunk_tasks{i_task});
                else
                    result = processing_function(chunk_tasks{i_task}, varargin{:});
                end
                
                % Store result directly
                all_results{chunk_start + i_task - 1} = result;
            end
        end
    end
    
    % Create task indices for reconstruction (if needed)
    if nargout > 1
        task_indices = (1:num_tasks)';
    end
    
    % (Removed nested progress update function)
    
end 