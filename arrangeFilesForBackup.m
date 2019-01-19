function arrangeFilesForBackup(flag)

try

% Define paths 
dirs.dropbox  = 'C:\Users\Dumbo\Dropbox\Data';
dirs.behavior = 'D:\Behavior';
dirs.webcam   = 'D:\Python_Videos';
dirs.TDT      = 'D:\UCL_Behaving';
dirs.RV2      = '\\192.168.0.102\data\recordings';

% Define flag if undefined
if ~exist('flag','var'), flag = 0; end

% Define ferrets of interest
ferrets = dir( fullfile(dirs.TDT, 'F1*'));

% Note today's date
currentTime = now + flag;
str = [datestr(currentTime, 'dd_mm_yyyy') '*'];

% List webcam videos
webcam.str   = [datestr(currentTime, 'yyyy-mm-dd') '*'];
webcam.files = dir( fullfile( dirs.webcam, sprintf('%s.avi', webcam.str)));
% webcam.files = dir( fullfile( dirs.webcam, sprintf('%s.mp4', webcam.str)));

for i = 1 : numel(webcam.files)
    
    Y = str2num(webcam.files(i).name(1:4));
    M = str2num(webcam.files(i).name(6:7));
    D = str2num(webcam.files(i).name(9:10));
    H = str2num(webcam.files(i).name(18:19));    
    MN = str2num(webcam.files(i).name(21:22));
    S = str2num(webcam.files(i).name(24:25));
    
   webcam.date(i) = datenum(Y,M,D,H,MN,S); 
end

if exist('webcam_copy','var')
    webcam.copy  = webcam_copy == 1;
else
    webcam.copy  = false;
end

% For each ferret
for i = 1 : numel(ferrets)
    
    % Extend paths
    dirs.ferret.behavior = fullfile( dirs.behavior, ferrets(i).name);
    dirs.ferret.dropbox  = fullfile( dirs.dropbox, ferrets(i).name);
    dirs.ferret.TDT = fullfile( dirs.TDT, ferrets(i).name);
    dirs.ferret.RV2 = fullfile( dirs.RV2, ferrets(i).name);
    
    % List behavioral files for today
    files = dir(fullfile( dirs.ferret.behavior, str));

    % List blocks
%     blocks = dir( fullfile( dirs.ferret.TDT, 'Block*'));
    
    % For each file
    for j = 1 : numel(files)            
        
        % Get files 
        block.idx = strfind(files(j).name,'Block');
        block.str = files(j).name(block.idx:end-4);
        
        % Extend paths
        dirs.block.TDT = fullfile( dirs.ferret.TDT, block.str);
        dirs.block.RV2 = fullfile( dirs.ferret.RV2, block.str);
        
        % Copy RV2 file
%         aviFile = dir( fullfile(dirs.block.RV2, '*.avi'));
%         srcFile = fullfile(dirs.block.RV2, aviFile.name);
%         tarFile = fullfile(dirs.block.TDT, aviFile.name);                
%         copyfile_myVersion(srcFile, tarFile)
        
        % Copy files to dropbox 
        srcFile = fullfile(dirs.ferret.behavior, files(j).name);
        tarFile = fullfile(dirs.ferret.dropbox,  files(j).name);        
        copyfile_myVersion(srcFile, tarFile)     
        
        % Identify webcam file
        % Get time difference between webcam files and behavioral files
        webcam.delta = abs(webcam.date - getBehavioralFileDate(files(j).name));
        [webcam.discrep, webcam.idx] = min(webcam.delta);
        webcam.target = webcam.files(webcam.idx).name;
        
        % Warn if min difference in file time is large (> 5 mins)
        if webcam.discrep > seconds2days(300)
            warning('Large difference in webcam timing')            
        end
        
        % Move webcam file
        srcFile = fullfile(dirs.webcam,    webcam.target);
        tarFile = fullfile(dirs.block.TDT, webcam.target);
        copyfile_myVersion(srcFile, tarFile)    
        
        
        % Move webcam matlab file
        webcam.target = strrep(webcam.target,'.avi','.txt');
        srcFile = fullfile(dirs.webcam,    webcam.target);
        tarFile = fullfile(dirs.block.TDT, webcam.target);
        copyfile_myVersion(srcFile, tarFile)   
    end
end

catch err
    err
    keyboard
end

function copyfile_myVersion(srcFile, tarFile)
% Copy file if it doesn't already exist
if ~exist(srcFile,'file')
    fprintf('%s does not exist\n', srcFile)
    return
end

if ~exist(tarFile,'file')
    copyfile( srcFile, tarFile,'f')
end