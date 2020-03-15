function verifyTrialOnset(varargin)
%
% Stephen Town
% 11th March 2020
%
%
% Inputs: (Optional: one of...)
% - Config file with list of tanks (ferrets) and where to find them
% - File path of:
%   - Tank for analysis of all blocks (e.g. "G:\UCL_Behaving\F1605_Snorlax")
%   - Block for specific analysis of that file (e.g. "G:\UCL_Behaving\F1605_Snorlax\Block_J2-10")
%
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


% Comments:
%   - Removed frontiers metadata generation (this is now done when summarizing
%     behavioral files.
%   - Control graphics to make optional
%   - Direct output files as required
%     Move output options to config file
%     Add response time verification (turn file into verifyTiming.m)
%   - label number of stimuli (inc. zero when stimuli are absent (e.g. L51
%   snoralx)

% Load configuration file and include TDT library in path
config_file = 'directory_config.json';
config = jsondecode(fileread( config_file));
addpath( genpath( config.tdt_lib))

set(0,'DefaultFigureColormap',eval('gray'));      % Make everything black and white (does this reduce png size?)
                   
%%% Parse inputs
if nargin == 0
    file_path = config_file;
else
    file_path = varargin{1};
end

% If using a configuration file
[tank, block, ext] = fileparts(file_path);

if contains(ext,'.json')                
    block_table = get_blocks_to_analyze( config);
end

% If targetting a specific tank / block
if isfolder( file_path)

    config.n_ferrets = 1;
    
    % If not directed to a specific block, assume that we're given a tank 
    % from which to look at all blocks
    if ~contains( block, 'Block')
        
        [config.tank_dir, config.ferrets] = fileparts( varargin{1});
        config.ferrets = {config.ferrets};
        block_table = get_blocks_to_analyze( config);
    else
        block_table = table( tank, block, 'variableNames',{'tank','name'});
    end       
end

%%% Run main function on each block
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')

for i = 1 : size( block_table, 1)
            
    main( block_table(i,:), config)
end

warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')


function main( block, config)

try
    
    % Skip if this file has already been processed
    if outputs_exist(block, config), return; end       

    % Get text file start times
    T = readtable( block.behav_file{1});

    if isempty(T)    
        fprintf('No trials for %s %s\n', block.ferret{1}, block.name{1})
        return
    end    

    % Load TDT data
    DOut = load_tdt_data( block, 'DOut', T.StartTime, config);
    Sens = load_tdt_data( block, 'Sens', T.RespTime, config);
    
    if size(Sens.data, 1) > 1   % Collapse across separate channels (response spouts)
        Sens.data = any(Sens.data);     
    end
    
    % Get signals surrounding events
    [startTraces, stimTvec] = getTraces( DOut, T.StartTime, config);    
    [respTraces, respTvec] = getTraces( Sens, T.RespTime, config);
    
    % Align traces (repeat to be sure)
    T.StartTimeCorrected = align_onset(startTraces, stimTvec, T.StartTime);
    [corrected, ~] = getTraces( DOut, T.StartTimeCorrected, config);
    T.StartTimeCorrected = align_onset(corrected, stimTvec, T.StartTimeCorrected);    
    
    % Apply correction and report actual stimulus number    
    [corrected, ~] = getTraces( DOut, T.StartTimeCorrected, config);
    T = correct_nStim( corrected, stimTvec, T);
            
    % Visually confirm correction
    if config.output.show_image
                
        fig = show_traces_as_im( stimTvec, startTraces, corrected, block);
        
        draw_trace(respTraces, respTvec, subplot(133)); 
        xlabel('Time w.r.t. Response (s)', 'FontSize', 8) 
        title('(Original)', 'FontSize', 8) 
        
        add_markers( T.StimMissing, [0.63 0.11], 'Missing')
        add_markers( T.Correct, [0.91 0.11], 'Reward')      
        
        % Save file if required
        if config.output.save_image
            
            save_path = get_output_path( block, config, 'png');           
            myPrint( save_path, 'png', 150)
            close(fig)
        end
    end
    
    % Save output
    if config.output.save_data           
        save_path = get_output_path( block, config, 'csv');                
        writetable( T, save_path, 'delimiter', ',') 
        
        save_path = get_output_path( block, config, 'mat');                
        save( save_path, 'stimTvec', 'corrected') 
    end

catch err
    err
    keyboard
end


function block_table = get_blocks_to_analyze( config)

% For each ferret
config.n_ferrets = numel( config.ferrets);
block_table = [];

for i = 1 : config.n_ferrets
    
    % List blocks    
    tank_path = fullfile( config.tank_dir, config.ferrets{i});
    blocks = dir( fullfile( tank_path, 'Block_J*'));
    
    % Add tank and ferret name to table and concatenate    
    if ~isempty(blocks)        
        blocks = struct2table( blocks);
        blocks.tank = repmat({tank_path}, size(blocks, 1), 1);
        blocks.ferret = repmat(config.ferrets(i), size(blocks, 1), 1);
        block_table = [block_table; blocks];
    end      
end

% Add path to behavioral file
if isfield(config, 'behav_dir')
   
    % For each block
    n_blocks = size(block_table, 1);
    block_table.behav_file = cell( n_blocks, 1);
    row_to_remove = [];
    
    for i = 1 : n_blocks
        
        % Search for behavioral file
        temp_path = fullfile( config.behav_dir, block_table.ferret{i});
        
        files = dir( fullfile( temp_path, ['* ' block_table.name{i} '.txt']));
        
        if numel(files) == 1
            block_table.behav_file{i} = fullfile( temp_path, files.name);
        else
           row_to_remove = [row_to_remove; i]; 
            
           warning('Found %d files for %s, %s', numel(files),...
               block_table.ferret{i}, block_table.name{i}) 
        end        
    end    
    
    % Remove blocks without behavioural data
    block_table( row_to_remove, :) = [];
end


function should_skip = outputs_exist(block, config)
% 
% Returns true if either the figure or data file associated with this block
% already exists in the save directory
          
should_skip = false;

if config.output.save_data        % Behavioural data 
    should_skip = report_existing_file( block, config, 'csv');    
end

if config.output.save_image        % Image figure   
    should_skip = report_existing_file( block, config, 'png');
end

if config.output.save_image        % Matlab mask    
    should_skip = report_existing_file( block, config, 'mat');
end


function should_skip = report_existing_file(block, config, file_type)

file_path = get_output_path( block, config, file_type);
should_skip = exist( file_path, 'file');

if should_skip
    fprintf('%s %s - already processed\n', block.ferret{1}, block.name{1})
end


function save_path = get_output_path( block, config, fileType)
   
% Get info from original file
[~, original_file] = fileparts(block.behav_file{1});
sep_idx = strfind( original_file, ' ');

level_str = original_file( sep_idx(1)+1 : sep_idx(2)-1);
day_str = original_file(1:sep_idx(1)-1);
time_str = original_file(sep_idx(2)+1:sep_idx(3)-1);

% Revise date format for name        
dn = datenum( [day_str ',' time_str],'dd_mm_yyyy,HH_MM_SS.FFF');
ds = datestr(dn, 'yyyy-mm-ddTHH-MM-SS');

save_name = [ds '_' block.name{1} '_' level_str '.' fileType];      

% Ensure save path is available
eval( sprintf('save_path = config.save_dir.%s;', fileType))
save_path = fullfile( save_path, block.ferret{1});
        
if ~isfolder( save_path)
    mkdir( save_path)
end

save_path = fullfile( save_path, save_name);


function S = load_tdt_data(block, store, ev_times, config)
%
% Input: 
%   - block: row from block table containing path to data
%   - start_times: events of interest for this block
%   - config: settings for time window for all analyses
%
% Ouput:
%   - DOut: struct containing time-varying signal and sample rate

% Remove events that didn't happen (stamped with -1 for no response)
ev_times( ev_times < 0) = []; 

% Compute the time window of the session to load in (makes code run faster)
t_end = max(ev_times) + max(config.window) + 1;
t_start = min(ev_times) + min(config.window) - 3;
t_start = max([0 t_start]); % Disallow negative values

% Load data
block_path = fullfile(block.tank{1}, block.name{1});

S = TDTbin2mat( block_path, 'STORE', store, 'T1', t_start, 'T2', t_end);
eval( sprintf('S = S.streams.%s;', store))
S.duration = numel(S.data) ./ S.fs;
S.endTime = S.startTime + S.duration;


function [traces, tvec] = getTraces( DOut, ev_times, config)
%
% Input:
%   - DOut: record of stimulus presentations made from the TDT
%   - startTime: recorded time of sound onset in text file 
%
% Output:
%   - traces: rectangular matrix of device records
 
% Check start times are compatible
rm_rows = ev_times - config.window(1) < 0;

if any(rm_rows)        
    keyboard
end

% Remove events that didn't happen (stamped with -1 for no response)
ev_times( ev_times < 0) = []; 

% Define the samples that make up each window and take from DOut 
windowSamps = round(config.window.* DOut.fs);
windowSamps = windowSamps(1) : windowSamps(2);

ev_times = ev_times - DOut.startTime;   % Offset for data not loaded
startSamps = round(ev_times .* DOut.fs);


windowIdx = bsxfun(@plus, startSamps, windowSamps);

if max(windowIdx(:)) > numel(DOut.data) % Edge case where block terminated quickly after last response
    warning('Padding end of signal with zeros')
    n_samples_missing = max(windowIdx(:)) - numel(DOut.data);
    padding = zeros(1, n_samples_missing);
    DOut.data = [DOut.data padding];
end

% Bug fix for edge case where we get exact timing match with zero delay
if sum(windowIdx(:) == 0) == 1 && windowIdx(1) == 0
    windowIdx(1) = 1;
end

% Get data and corresponding time vector
traces = DOut.data(windowIdx) > 0;      % Return as binary
tvec = windowSamps ./ DOut.fs;          % Add time vector


function new_times = align_onset(traces, tvec, original_times)
%
% Input:
%   traces: alignment data
%  
%
% Output:
%  

[~, idx] = max(traces, [], 2);  % Find first positive value
offset = tvec(end, idx);    % Get time value 
new_times = original_times + transpose(offset);    % Apply offset




function T = correct_nStim( traces, tvec, T)
    
% Ask if there was anything in the first second 
t_idx = tvec >= 0 & tvec < 1;
T.StimMissing = all(traces(:, t_idx) == 0, 2);
    
%%% ADD CODE to detect excactly how many stimuli there were

    
    % Remove trials with failed duration (error of > 2 ms)
    
    % Draw stimulus delivery    
%     tError = abs(Duration - 1.70);
%     T(tError > 0.002,:) = [];

function fig = show_traces_as_im( tvec, traces, corrected, block)
        
fig = figure('position',[2680 560 870 420]);

draw_trace(traces, tvec, subplot(131)); 
[~, original_file] = fileparts(block.behav_file{1});
title( strrep(original_file,'_',' '),'FontWeight','Normal','FontSize',8)
xlabel('Time w.r.t. Start (s)')

draw_trace(corrected, tvec, subplot(132)); 
title('Corrected','FontWeight','Normal','FontSize',8)
xlabel('Time w.r.t. Start (s)')
        

function draw_trace(traces, tvec, ax)
%
% Draw heatmap of binary signal
axes(ax)
imagesc(tvec, 1:size(traces, 1), traces)
ylabel('Trial')
axis('tight')
set(gca,'fontsize',8)


function add_markers( trial_val, pos, str)
        
axes('position',[pos  0.01 0.815])      

y = find(trial_val);
scatter( zeros(size(y)), y, '.k')

set(gca,'color','none',...
    'xcolor','none',...
    'ycolor','none',...
    'xlim',[-1 1],...
    'ylim',[0 numel(trial_val)]+0.5,...
    'ydir','reverse')

title(str,'FontAngle','italic','FOntSize',7,'FontWeight','normal')
