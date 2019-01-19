function getSessionLFPs_wavelet
%
% Loads snippets from extracted broadband trace. There are n snippets,
% where n is the number of trials. Each snippet is m samples long, where m
% equals the snippet duration x sampling rate.
%
% Plots spectra for each channel together on one channel
% Plots spectrograms for each channel on separate axes in a second figure
%
% Requires
%   - Chronux library
%   - Electrode data in Matlab format (McsRecording*.mat)
%   - Behavioral files with corrected onset times
%   - Synchronization file

% Add paths
addpath( genpath( 'C:\Users\Dumbo\Documents\MATLAB\lib\chronux_2_11')) 
saveDir = 'D:\Frontiers Data Analysis\Session LFPs';
syncDir = 'D:\Frontiers Data Analysis\Timing\Digital';

% List subjects
rootDir = 'C:\Data\UCL_Behaving';
ferrets = dir( fullfile(rootDir,'F*'));

% For each ferret
for i = 2 : length(ferrets)
    
    % Extend directories
    ferDir     = fullfile( rootDir, ferrets(i).name);
    ferSyncDir = fullfile( syncDir, ferrets(i).name);
    saveFerDir = fullfile(saveDir, ferrets(i).name);
    
    % List blocks
    blocks = dir( fullfile( ferDir, 'Block*'));
    
    % Create save block
    if ~isdir( saveFerDir), mkdir( saveFerDir); end
    
    % For each block
    for j = 1 : numel(blocks)
        
        % Get sync data
        [fS, stimOnsets] = getMCS_sampleRate( ferSyncDir, blocks(j).name);                
                
        % List recording data files
        blockDir  = fullfile( ferDir, blocks(j).name);
        blockSaveDir = fullfile( saveFerDir, blocks(j).name);
        recFiles  = dir( fullfile( blockDir, '*McsRecording*.mat'));
        
        if ~isdir(blockSaveDir)
            mkdir(blockSaveDir)
        end
        
        % Load behavior 
        T = loadBehavior( blockDir);
        
        % Identify trial onsets (as opposed to the broad range of events in
        % the digital out)
                
        % For each recording file
        for k = 1 : numel(recFiles)
            
            % Note headstage (from which we can cross ref location later)
            str = recFiles(k).name(end-7:end-4);
            
            % Check if file exists            
            saveName = sprintf('%s_%s_%s.fig', ferrets(i).name, blocks(j).name, str); 
            savePath = fullfile(blockSaveDir, saveName);
            
            if exist( savePath,'file')
                fprintf('%s exists - skipping\n', saveName); continue
            else
                fprintf('%s in progress\n', saveName)
            end
            
            % Load data
            recPath = fullfile( blockDir, recFiles(k).name);        
            load( recPath, 'Raw_Electrode_1')         

            % For every chanel
            nChans = size(Raw_Electrode_1.ChannelData,2);
            
            for chan = 1 : nChans
                
                % User feedback
                fprintf('\tRunning channel %02d\n', chan) 

                % Run main function
                figs = formatData( Raw_Electrode_1.ChannelData(:,chan), stimOnsets, fS);             
                
                % Save figure
                saveChanPath = strrep(savePath,'.fig',sprintf('_%02d.fig', chan));
                saveas(figs(1), strrep(saveChanPath,'.fig','_ITC_Spectra.fig'))
                saveas(figs(2), strrep(saveChanPath,'.fig','_ERP.fig'))
                close(figs)
            end
        end
    end
end



function T = loadBehavior( blockDir)

    behavFile = dir( fullfile( blockDir, '*level*.mat'));

    % Sanity check that we don't have multiple behavioral files
    if numel(behavFile) > 1, keyboard; end

    % Load behavior
    load( fullfile( blockDir, behavFile(1).name),'T')


function [fS, onsets]  = getMCS_sampleRate( myDir, myBlock)
        
    syncFiles  = dir( fullfile( myDir, sprintf('*%s*.mat', myBlock)));
    nSyncFiles = numel(syncFiles);

    fS = nan(nSyncFiles,1);
    onsets = cell(nSyncFiles,1);

    for s = 1 : nSyncFiles

        S = load( fullfile( myDir, syncFiles(s).name));
        fS(s) = S.fS.MCS;
        onsets{s} = S.onsets.MCS;
    end

    onsets = cell2mat(onsets);
    onsets = mean(onsets,2);
    fS = mean(fS);

    
function f = formatData( data, trig, fS)

% Options
window = [-1 1];

% Define window index
window_samps = round(window .* fS); 
window_samps = window_samps(1) : window_samps(2);
window_time  = linspace(window(1),window(2),numel(window_samps));

trigSamps = bsxfun(@plus, trig, window_samps);

trigSamps(any(trigSamps < 1, 2),:) = [];            % Remove windows before stim
trigSamps(any(trigSamps > numel(data), 2),:) = [];  % Remove windows after stim

trigData = single(data(trigSamps));

% Get ERP
erp_mean = mean(trigData,1);
erp_std  = std(trigData,[],1);
erp_n    = size(trigData,1);
erp_sErr = erp_std./ erp_n;


% Get ITC and LFP
[mySpect, ITC, f(1)] = get_Spectra_and_ITC(trigData, fS, window_time);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw
f(2) = figure; %'Name',recFiles(k).name);
sp = axes('nextPlot','add');
plotSE( window_time, erp_mean, erp_sErr, sp(1));
plotYLine(0,sp(1));
xlabel('Time (s)')
ylabel('Voltage')

