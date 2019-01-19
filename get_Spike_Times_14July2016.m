function get_Spike_Times
%
% Loads snippets from extracted broadband trace. There are n snippets,
% where n is the number of trials. Each snippet is m samples long, where m
% equals the snippet duration x sampling rate.
%
% Plots spectra for each channel together on one channel
% Plots spectrograms for each channel on separate axes in a second figure

% Define paths
timeDir = 'D:\Frontiers Data Analysis\Timing\Digital';
saveDir = 'D:\Frontiers Data Analysis\Spikes';

% List subjects
rootDir = 'C:\Data\UCL_Behaving';
ferrets = dir( fullfile(rootDir,'F*'));

% Anon functions
isPFC = @(x) ~isempty( strfind(x,'K082'));
isAC  = @(x) ~isempty( strfind(x,'K080'));

% For each ferret
for i = 2 : length(ferrets)
    
    % List blocks
    ferDir = fullfile( rootDir, ferrets(i).name);
    blocks = dir( fullfile( ferDir, 'Block_J2-4*'));
    
    % Make save directory 
    saveDir_i = fullfile( saveDir, ferrets(i).name);
    if ~isdir(saveDir_i), mkdir(saveDir_i); end
    
    % For each block
    for j = 1 : numel(blocks)
        
        % List recording data files
        blockDir = fullfile( ferDir, blocks(j).name);       
        recFiles = dir( fullfile( blockDir, '*McsRecording*.mat'));
        
        % Make save directory
        saveDir_j = fullfile( saveDir_i, blocks(j).name);
        if ~isdir(saveDir_j), 
            mkdir(saveDir_j); 
%         else
%             fprintf('%s already complete - skipping\n', blocks(j).name); continue
        end
                
        % For each recording file
        for k = 1 : numel(recFiles)
            
            % Infer location from file name 
            if isPFC(recFiles(k).name), str = 'PFC'; end
            if isAC(recFiles(k).name),  str = 'AC'; end
            
            % Check if timing file exists            
            timeName = sprintf('%s_%s_%s.mat', ferrets(i).name, blocks(j).name, str); 
            timePath = fullfile(timeDir, timeName);
            
            if ~exist( timePath,'file')
                fprintf('%s does not exist - skipping\n', timeName); continue
            else
                fprintf('%s in progress\n', timeName)
                load(timePath,'P')
            end
            
            % Load data
            recPath = fullfile( blockDir, recFiles(k).name);        
            load( recPath, 'Filter_1')
            
            % For each channel                        
            for chan = 1 : 16             
                                
                spikeFile = sprintf('%s_Chan_%02d.mat', str, chan);   
                spikePath = fullfile( saveDir_j, spikeFile);
                
                if exist(spikePath,'file')
                    fprintf('%s %s already complete - skipping\n', recFiles(k).name, spikeFile); continue
                end
                
                % Get spike times
                [t, wv] = getSpikeTimes(P, Filter_1.ChannelData(:,chan));
                
                % Save times             
                save( spikePath, 't','wv');
            end                                                    
        end
    end
end


function [t, wv] = getSpikeTimes(b, x)

% This is taken from getMClustEvents_AlignedInterpolated with the threshold
% parameters adjusted for the different electrodes.
%
% x = filtered voltage trace from which to extract spikes
% b = regression coefficients describing offset and sampling rate
% difference between devices

% Parameters
wInt = 1;
interpFactor = 4;

interpInt = wInt / interpFactor;    % Method dependent (i.e. if you interpolate or not)
window = -15 : wInt : 16;
interpWind = -15 : interpInt  : 16;

nW = numel(window)+1;               % These are regardless of method (interpolated or not)
alignmentZero = find(window == 0);

% Preassign
[t, wv] = deal([]);

% Format
if iscolumn(x), x = transpose(x); end
x = single(x);

% Get upper (ub) and lower (lb) bounds
lb = -4 * std(x);
ub = -10 * std(x);

% Identify thrshold crossings
lcIdx = find(x < lb);
ucIdx = find(x < ub);

% Remove events exceeding the upper threshold                    
lcIdx = setdiff(lcIdx,ucIdx);                                   %#ok<*FNDSB>

% Move to next trial if no events were found
if isempty(lcIdx); return; end

% Identify crossing points in samples
crossThreshold = lcIdx([0 diff(lcIdx)]~=1);

% Remove events where window cannot fit
crossThreshold(crossThreshold < nW) = [];
crossThreshold(crossThreshold > (length(x)-nW)) = [];

% Make row vector
if iscolumn(crossThreshold),
    crossThreshold = crossThreshold';
end

% Get interim waveforms
wvIdx  = bsxfun(@plus, crossThreshold', window);
wv     = x(wvIdx);

% Move to next trial if no waveforms are valid
if isempty(wv); return; end

% Interpolate waveforms
wv = spline(window, wv, interpWind);

% Align events
[~, peakIdx] = min(wv,[],2); 
peakIdx = round(peakIdx / interpFactor);     % Return interpolated peakIdx to original sample rate
alignmentShift = peakIdx' - alignmentZero;
alignedCrossings = crossThreshold + alignmentShift;

% Reset events where window cannot fit (i.e. don't
% throw away, just include without alignment)
alignedCrossings(alignedCrossings < nW) = crossThreshold(alignedCrossings < nW);                     
alignedCrossings(alignedCrossings > (length(x)-nW)) = crossThreshold(alignedCrossings > (length(x)-nW));

% Make row vector
if iscolumn(alignedCrossings),
    alignedCrossings = alignedCrossings';
end

% Get waveforms
wvIdx = bsxfun(@plus, alignedCrossings', window); % But sample aligned waveforms
wv    = x(wvIdx);
               
% Calculate time in TDT from wireless MCS samples
t = (crossThreshold - b(2)) ./ b(1); % y = a + bx => x = (y-a)/b

