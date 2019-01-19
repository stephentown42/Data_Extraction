function verifyTrialOnset


% List subjects
rootDir = 'C:\Data\UCL_Behaving';
ferrets = dir( fullfile(rootDir,'F*'));

% For each ferret
for i = 1 : length(ferrets)
    
    % List blocks
    ferDir = fullfile( rootDir, ferrets(i).name);
    blocks = dir( fullfile( ferDir, 'Block*'));
    
    % For each block
    for j = 1 : numel(blocks)
               
        % Run main function
        main(ferrets(i).name, blocks(j).name)
    end
end




function main(ferret, block)

if nargin == 0
    
    ferret = 'F1506_Phoenix';
    block  = 'Block_J2-49';
end

% Define paths 
txtDir = fullfile('C:\Data\Behavior', ferret);
tdtDir = fullfile('C:\Data\UCL_Behaving', ferret, block);
saveDir = fullfile('D:\Frontiers Data Analysis\Timing\StimOnsets', ferret);
saveBeh = fullfile('D:\Frontiers Data Analysis\Behavior',ferret);

if ~isdir(saveDir), mkdir(saveDir); end
if ~isdir(saveBeh), mkdir(saveBeh); end

% Check for existing files
exFiles = dir( fullfile(saveDir, sprintf('*%s.*', block)));

if ~isempty(exFiles)
    fprintf('%s %s already processed\n', ferret, block)
    return
end


% Find text file
txtFile = dir( fullfile(txtDir, sprintf('*%s.txt', block)));

if isempty(txtFile)
    fprintf('Could not find text file for %s %s\n', ferret, block)
    return
end

txtFile = txtFile(1).name;

% Get text file start times
B = importdata( fullfile(txtDir, txtFile));

if ~isfield(B,'data')    
    fprintf('No trials for %s %s\n', ferret, block)
    return
end

% LEDs      = B.data(:,strcmp(B.colheaders,'LED Location'));
% Sounds    = B.data(:,strcmp(B.colheaders,'LED Location'));

% Format header text
headers = B.colheaders;
removeQ = @(x) strrep(x,'?','');
removeS = @(x) strrep(x,' ','');
headers = cellfun(removeQ, headers,'un',0);
headers = cellfun(removeS, headers,'un',0);

% Convert to table
T = array2table(B.data,'VariableNames',headers);

T.StartTime = T.StartTime;

nTrials   = numel(T.StartTime);

% Load TDT data
load( fullfile( tdtDir, sprintf('%s_TDT_data.mat',block)),'DOut')

% Create TDT time vector as samples
txtSamps = T.StartTime .* DOut.fs;

% Get TDT traces at start times
window      = [-1.5 2]; 
windowSamps = window.* DOut.fs;
windowTime  = linspace(window(1), window(2), 1+diff(windowSamps));
traces      = cell(nTrials, 1);

for i = 1  : nTrials
    
    windowIdx = round(txtSamps(i)+windowSamps);
    traces{i} = DOut.data(windowIdx(1):windowIdx(2));    
end

% Convert to constant length
nMin = min( cellfun(@numel, traces));

for i = 1 : nTrials, 
   traces{i} = traces{i}(1:nMin); 
end

traces = cell2mat(traces);
trialIdx = 1:nTrials;

% Draw stimulus delivery
f = figure('name', txtFile);
subplot(1,2,1)
imagesc(windowTime, trialIdx, traces)
drawnow



% Get last positive event to identify manual plays
[stimStart, stimEnd] = deal( nan(nTrials, 1));

for i = 1 : nTrials
    
   stimOn = find(traces(i,:) == max(traces(i,:)));
   
   stimStart(i) = windowTime(min(stimOn));
   stimEnd(i)   = windowTime(max(stimOn));
end


% Base trial mask on t-test
p    = normcdf(stimStart, mean(stimStart), std(stimStart));
isOk = abs(p-0.5) < 0.1;

subplot(122)
imagesc(windowTime, trialIdx(isOk), traces(isOk,:))
drawnow

% Apply mask
T.CorrectedStartTime = T.StartTime + stimStart;
T.CorrectedStartTime(~isOk) = nan;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Add computed variables to data
%
% Reconstruct dominant modality
%   Dominance assessed in two steps: (1) for each modality separately (2) across modalities
domMod(:,1) = (T.LEDLocation == T.TargetSpout) .* 1;    % Visual = 1, not visual = 0
domMod(:,2) = (T.SpeakerLocation == T.TargetSpout).*2;  % Not auditory = 0, Auditory = 2
T.domMod    = sum(domMod,2)-1;                          % Bring together across modalities
        
% Identify trial classes
trialClass = T.Modality + T.domMod.*10;       % 0 = V; 2 = VA; 11 = A; 12 = AV
T.V_Trial  = trialClass == 0;   % Visual only
T.VA_Trial = trialClass == 2;   % Audiovisual - target visual
T.A_Trial  = trialClass == 11;  % Auditory only
T.AV_Trial = trialClass == 12;  % Audiovisual - target auditory

% Get spout angles
angles = [-30 -60 -90 -120 -150 -180 150 120 90 60 30 0]; % Angles on the clockface
T.TargetAngle = transpose( angles( T.TargetSpout));

% Get catch spout 
%   defined as the location (2,10 or 12 o'clock) where neither visual or
%   auditory stimulus is presented.
T.CatchSpout = T.SpeakerLocation + T.LEDLocation;   % The sum of locations uniquely defines where the catch spout is
T.CatchSpout(T.CatchSpout == 14) = 10;              % Stimuli presented at 2 & 12 o'clock (2+12 = 14) - catch spout must be 10
T.CatchSpout(T.CatchSpout == 22) = 2;               % Stimuli presented at 10 & 12 o'clock (10+12 = 22) - catch spout must be 2
% Stimuli presented at 2 and 10 o'clock (10+2 = 12) - catch spout is already 12



% Save
save(fullfile(saveBeh, strrep(txtFile,'.txt','.mat')),'traces','T')
saveas(f, fullfile(saveDir,strrep(txtFile,'.txt','.fig')))
close(f)

