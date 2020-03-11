function verifyTrialOnset_keepCorrectionTrials
%
% Stephen Town
% 06 March 2017
%
% Inputs:
% - List of tanks (ferrets) with associated file trees for:
%   - behavioural files (text files generated from GoFerret)
%   - output directory (for saving data in a .mat and .fig) combination
%   - stim onset (unknown - can probably be deleted)
%
% Formatting requirements: Assumes that TDT data has been converted into
% .mat format 
%
% Objectives:
%   - Correct start times to the actual beginning of a stimulus
%   presentation sequence
%   - Add fields to behavioral data containing trial masks for particular
%   conditions used in analysis for the frontiers project
%
% Background:
% For reasons that have never been solved, poor center spout holding can
% lead to odd sequences of stimulus presentation that I have never been
% able to debug. This means that some of the start times that would have
% been reliable in other projects are not always bulletproof and must be
% checked. 
%
% Additionally in the frontiers project, stimuli are presented as sequences
% of multiple stimuli, which the animal may respond to before the end of
% the sequence. It is thus useful to know what stimulus was presented and
% when it was terminated 


% List subjects
rootDir = 'F:\UCL_Behaving';
ferrets = dir( fullfile(rootDir,'F*'));

% For each ferret
for i = 9 : length(ferrets)
    
    % List blocks
    ferret = ferrets(i).name;
    ferDir = fullfile( rootDir, ferret);
    blocks = dir( fullfile( ferDir, 'Block_J*'));
    
    % Specify paths
    dirs.behavSrc  = fullfile('F:\Data\Behavior', ferret);
    dirs.stimOnset = fullfile('E:\Frontiers Data Analysis\Timing\StimOnsets', ferret);
    dirs.behavSave = fullfile('E:\Frontiers Data Analysis\Behavior_All',ferret);
       
    % Create directories 
    if ~isdir(dirs.stimOnset), mkdir(dirs.stimOnset); end
    if ~isdir(dirs.behavSave), mkdir(dirs.behavSave); end
    
    % For each block
    for j = 1 : numel(blocks)
               
        % Run main function
        main(ferrets(i).name, blocks(j).name, dirs)
    end
end



function main(ferret, block, dirs)

try
    % Define paths 
    tdtDir  = fullfile('F:\UCL_Behaving', ferret, block);

    % Check for existing files
    exFiles = dir( fullfile(dirs.behavSave, sprintf('*%s.*', block)));

    if ~isempty(exFiles)
        fprintf('%s %s already processed\n', ferret, block)
        return
    end

    % Find text file
    txtFile = dir( fullfile(dirs.behavSrc, sprintf('*%s.txt', block)));

    if isempty(txtFile)
        fprintf('Could not find text file for %s %s\n', ferret, block)
        return
    end

    txtFile = txtFile(1).name;

    % Get text file start times
    B = importdata( fullfile(dirs.behavSrc, txtFile));

    if ~isfield(B,'data')    
        fprintf('No trials for %s %s\n', ferret, block)
        return
    end

    % Format header text
    headers = B.colheaders;
    removeQ = @(x) strrep(x,'?','');
    removeS = @(x) strrep(x,' ','');
    headers = cellfun(removeQ, headers,'un',0);
    headers = cellfun(removeS, headers,'un',0);

    % Convert to table
    T = array2table(B.data,'VariableNames',headers);

    % Load TDT data
    load( fullfile( tdtDir, sprintf('%s_TDT_data.mat',block)),'DOut')

    % Draw stimulus delivery
    f = figure('name', txtFile);
    sp = dealSubplots(1,2);

    % Get TDT traces at start times
    [traces, t] = getTraces( DOut, T.StartTime);

    [stimStart, Duration] = drawStimOnsets(sp(1), t, traces);

    % Correct for offsets
    T.CorrectedStartTime = T.StartTime + stimStart; 

    % Remove trials with failed duration (error of > 2 ms)
%     tError = abs(Duration - 1.70);
%     T(tError > 0.002,:) = [];

    % Recalculate start times
    [new_traces, t] = getTraces( DOut, T.CorrectedStartTime);

    % Escape if no trials left
    if isempty(T), 
        set(sp(2),'color','r'); 
        close(f); return; 
    end
    
    % Draw
    [stimStart, stimEnd] = drawStimOnsets(sp(2), t, new_traces);


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
    
    % Response modality
    T.A_Response = T.SpeakerLocation == T.Response;
    T.V_Response = T.LEDLocation == T.Response;    

    % Save
    saveFile = fullfile(dirs.behavSave, strrep(txtFile,'.txt','_all.mat'));
    save(saveFile,'traces','T')
    saveas(f, strrep(saveFile,'.mat','.fig'))
    close(f)

catch err
    err
    keyboard
end


function [traces, t] = getTraces( DOut, startTime)
%
% DOut is the record from the device
% startTime is the time recorded in the text file for stimulus onset
%
% traces is a rectangular matrix of device records
 

% Create TDT time vector as samples
startSamps = startTime .* DOut.fs;

% Analysis window
window = [-1.5 2]; 
windowSamps = window.* DOut.fs;
nSamps = 1+diff(windowSamps); 
t = linspace(window(1), window(2), nSamps);

nTrials     = numel(startSamps);
traces      = cell(nTrials, 1);

for i = 1  : nTrials
    
    windowIdx = round(startSamps(i)+windowSamps);
    traces{i} = DOut.data(windowIdx(1):windowIdx(2));    
end

% Convert to constant length
nMin = min( cellfun(@numel, traces));

for i = 1 : nTrials, 
   traces{i} = traces{i}(1:nMin); 
end

traces = cell2mat(traces);


function [stimStart, Duration] = drawStimOnsets(ax, t, traces)

% Draw PSTH
nTrials = size(traces,1);
imagesc(t, 1:nTrials, traces, 'parent', ax)

% Get last positive event to identify manual plays
[stimStart, stimEnd] = deal( nan(nTrials, 1));

for i = 1 : nTrials
    
   stimOn = find(traces(i,:) == nanmax(traces(i,:)));
   
   stimStart(i) = t(nanmin(stimOn));
   stimEnd(i)   = t(nanmax(stimOn));      
end

Duration = stimEnd - stimStart;

% Draw edges
plot3(stimStart, 1:i, repmat(3,i,1), 'w', 'parent', ax)
plot3(stimEnd, 1:i, repmat(3,i,1), 'y', 'parent', ax)

axis(ax,'tight')
