function get_Extraction_Metadata

% Get list of tanks
tankDir = 'C:\Data\UCL_Behaving';
ferrets = dir( fullfile(tankDir, 'F*'));
[S, T]  = deal([]);        % Preassign
        
% For each ferret
for i = 1 : numel(ferrets)
             
    % run main function
    [Si, Ti] = main( tankDir, ferrets(i).name);
    
    % Add to cell arrays
    S = [S; Si];
    T = [T; Ti];    
end

% Get rid of short blocks (indicated by empty cells in ferret column)
S = S(cellfun(@ischar, S(:,1)),:);
T = T(cellfun(@ischar, T(:,1)),:);


% Define headers for saving
headers  = {'Ferret','Block','Date','Start_Time','Duration',...
            'Wireless_Original','Wireless_Extracted','RV2_Video',...
            'Webcam_Video','TDT_mat','Wireless_mat','SpikeFiles','TimingFiles'};
        
TDTheaders = {'Ferret','Block','TDT Stores','','','','','','','','','',''};           


% Write to file
savePath = fullfile( tankDir, sprintf('Frontiers_Metadata_%s.xlsx', datestr(now,'dd_mmm_yyyy')));
xlswrite( savePath, [headers; S],'Metadata')
xlswrite( savePath, [TDTheaders; T],'TDT Extraction')



function [S, T] = main(tankDir, ferret)

% Define paths
spikeDir = fullfile('D:\Frontiers Data Analysis\Spikes', ferret);
timeDir  = 'D:\Frontiers Data Analysis\Timing\Digital';

% Connect to TDT    
TTfig = figure('visible','off');
TT    = actxcontrol('TTank.X'); 
tank  = fullfile(tankDir, ferret);

TT.ConnectServer('Local','Me'); 
TT.OpenTank(tank , 'R' );   % R-Read, W-Write, C-Control, M-Monitor    

% get block list
files  = dir( fullfile( tank, 'Block_J*'));

% Preassign
nFiles = length(files);
[S, T] = deal(cell(nFiles,13));


% for each block
for i = 1 : nFiles
           
    % define variables
    [S{i,1}, T{i,1}] = deal(ferret);
    [S{i,2}, T{i,2}] = deal(files(i).name);       % Block
    blockDir = fullfile(tank, S{i,2});
    saveName = fullfile( blockDir, sprintf('%s_extracted.mat', S{i,2}));
    
    % check for existing file
    if exist(saveName,'file')
        fprintf('\t%s exists - skipping\n', saveName); continue
    end
    
    % open block
    if ~TT.SelectBlock( S{i,2}),
        fprintf('Could not open %s\n',S{i,2}); 
        continue
    end        
    
    % Get block start and time
    startT = TT.CurBlockStartTime;
    stopT  = TT.CurBlockStopTime;    
    S{i,3} = TT.FancyTime(startT ,'D/O/Y'); % Start time
    S{i,4} = TT.FancyTime(startT ,'H:M:S.U'); % Start hour
    S{i,5} = stopT - startT;
    
    % Skip if shorter than two minutes
    if S{i,5} < 200, 
        fprintf('Short block: %s %s\n', S{i,1}, S{i,2}) % ferret, block
        
        [S{i,1}, T{i,1}] = deal(0); 
        continue 
    end
                
    % Check for other extracted files
    S{i,6}  = dirOrNone( fullfile( blockDir, '*.html'));            % MCS_original
    S{i,7}  = nDir( fullfile( blockDir, '*.h5'));              % MCS_extracted
    S{i,8}  = dirOrNone( fullfile( blockDir, '*Vid0*'));            % RV2_video
    S{i,9}  = dirOrNone( fullfile( blockDir, '*Block*.avi'));       % Matlab_video
    S{i,10} = dirOrNone( fullfile( blockDir, '*TDT_data*'));        % TDT_extracted
    S{i,11} = nDir( fullfile( blockDir, '*McsRecording*.mat')); % wireless_Mat            
    S{i,12} = nDir( fullfile( spikeDir, S{i,2}, '*.mat'));        % Get number of spike files
    S{i,13} = nDir( fullfile( timeDir, sprintf('%s_%s_*',ferret, S{i,2})));  % Get timing files
    
    % Get fields in TDT data
    tdtVars = matwho(fullfile(blockDir, S{i,10}));        
    
    for j = 1 : numel(tdtVars)  % Assign to cell array
       T{i,j+2} = tdtVars{j}; 
    end     
    
    % Reformat data
    S{i,2}  = strrep(S{i,2},'Block_','');
    S{i,10} = nDir( fullfile( blockDir, '*TDT_data*'));        % TDT_extracted
end

close(TTfig)


function n = nDir( filePath)

n = dir( filePath);
n = numel(n);


function str = dirOrNone( filePath)

str = dir(filePath);

if isempty(str)
    str = 'none';
else
    str = str(1).name;
end
