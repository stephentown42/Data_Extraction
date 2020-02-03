function put_MCS_data_in_tanks
%
% This function searches the Multichannel experimenter and Data Manager
% directories for data where h5 conversion process is complete. These data
% are then moved to the corresponding block on the relevant tank.
%
% Stephen Town 22nd March 2019

% Options
min_time_delta = 15; % Seconds

% Define paths
dirs.MCS_experimenter = 'C:\Users\steph\Documents\Multi Channel Systems\Multi Channel Experimenter';
dirs.MCS_data_manager = 'E:\MCS_DataManager_h5';
dirs.tanks = 'E:\UCL_Behaving';

% Build Block table
fprintf('Building block table\n')
block_table = build_block_table(dirs.tanks);

% List h5 files to consider
[n_h5, h5_files] = nDir( dirs.MCS_data_manager, '*.h5');

% Escape if nothing to be done
if n_h5 == 0
    fprintf('\tNo files to be moved\n')
    return
end

% Otherwise, get the time stamps of each file
h5_files = struct2table(h5_files);

% For each h5 file
for i = 1 : n_h5
   
    % Strip date from file name
    h5_datetime = datetime( h5_files.name{i}(1:19), 'InputFormat', 'yyyy-MM-dd''T''HH-mm-ss');
    h5_datenum = datenum( h5_datetime);
    
    % Get time difference 
    time_delta = abs(block_table.datenum - h5_datenum); % In days
    time_delta = time_delta .* 86400;   % In seconds
    
    % Find blocks below threshold
    [min_t, m_idx] = min(time_delta);
    
    % Warn if minimum time difference is too big
    if min_t > min_time_delta
        fprintf('Could not find block in time window for %s\n', h5_files.name{i})
        continue
    end
    
    % Otherwise, set destination path
    dirs.target = fullfile( dirs.tanks, block_table.Ferret{m_idx}, block_table.Block{m_idx});
    fprintf('Block found: %s\n', dirs.target)
    
    % Find primary data
    [n_pd, orig_files] = nDir( dirs.MCS_experimenter, strrep( h5_files.name{i}, '.h5', '*'));
    
    % For each primary file
    for j = 1 : n_pd
         
       % Move file
        my_src = fullfile( orig_files(j).folder, orig_files(j).name);
        my_dest = fullfile( dirs.target, orig_files(j).name);
        
        fprintf('\tMoving: %s\n', orig_files(j).name) % Update user        
        [status, msg, msgID] = movefile(my_src, my_dest);
        
        % Check for errors
        if status ~= 1
           fprintf('WARNING: CHECK FILE TRANSFER\nFrom: %s\nTo: %s\n', my_src, my_dest) 
        end           
    end
    
    % Move exported data
    my_src = fullfile( h5_files.folder{i}, h5_files.name{i});
    my_dest = fullfile( dirs.target, h5_files.name{i});
    
    fprintf('Moving: %s\n', h5_files.name{i}) % Update user
    [status, msg, msgID] = movefile(my_src, my_dest);
    
    % Check for errors
    if status ~= 1
        fprintf('WARNING: CHECK FILE TRANSFER\nFrom: %s\nTo: %s\n', my_src, my_dest)
    end
end
