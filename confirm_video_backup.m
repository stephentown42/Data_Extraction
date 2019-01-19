function confirm_video_backup
% 
% This function runs through a directory in which video files are initially
% recorded, and then copied to the relevant tank using the
% arrangeFilesForBackup.m function.
%
% However, it's possible that this backup may go wrong or not be done at
% all, and so before deleting critical data, we want to go through a
% confirmation step in which files are highlighted for deletion by
% confirming they exist in a tank and in the directory that is being
% cleared. 
%
% ST: 05 Dec 2018


% Define paths
dirs.to_clear = 'D:\Python_Videos';
dirs.tanks = 'D:\UCL_Behaving';
dirs.container = fullfile(dirs.to_clear, 'duplicates');

% List tanks
tanks = dir( fullfile(dirs.tanks, 'F*'));

% Create log
logName = datestr(now,'YYYY_mm_dd_backup_log_HH_MM.txt');
fid = fopen( fullfile(dirs.container, logName), 'wt+');
    
% For each tank
for i = 1 : numel(tanks)

    % Log
    fprintf(fid,'%s\n', tanks(i).name);
    
    % Extend paths
    dirs.test_tank = fullfile(dirs.tanks, tanks(i).name);

    % List blocks
    blocks = dir( fullfile( dirs.test_tank, 'Block*'));
    
    % For each block
    for j = 1 : numel(blocks)
    
        % Log
        fprintf(fid,'\t%s\n', blocks(j).name);
        
        % Extend path
        dirs.block = fullfile( dirs.test_tank, blocks(j).name);
            
        % Run comparison function
        main(dirs, '.txt', fid)
        main(dirs, '.avi', fid)
        main(dirs, '.mp4', fid)
    end    
end

% Close log
fclose(fid);


function main(dirs, ext, fid)
        
% List files of available types
[n_files, files] = nDir( dirs.block, ['*' ext]);

% Return if nothing to consider
if n_files == 0, return; end

% For each file discovered
for i = 1 : n_files
    
    % See if it exists in the folder under consideration for deletion
    test_path = fullfile( dirs.to_clear, files(i).name); 
    
    % If found
    if exist(test_path, 'file')
        
        % Report
        fprintf(fid, '\t\t%s\t', files(i).name);
        
        % Move file into deletion container
        dest_path = fullfile(dirs.container, files(i).name);
        [move_ok, ~, ~] = movefile( test_path, dest_path);
        fprintf(fid, '%d\n', move_ok);
    end
end