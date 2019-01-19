function compare_digital_out_times
%
% Loads snippets from extracted broadband trace. There are n snippets,
% where n is the number of trials. Each snippet is m samples long, where m
% equals the snippet duration x sampling rate.
%
% Plots spectra for each channel together on one channel
% Plots spectrograms for each channel on separate axes in a second figure

% Add paths
addpath( genpath( 'C:\Users\Dumbo\Documents\MATLAB\lib\chronux_2_11')) 
saveDir = 'D:\Frontiers Data Analysis\Timing\Responses';

% List subjects
rootDir = 'C:\Data\UCL_Behaving';
ferrets = dir( fullfile(rootDir,'F*'));

% Anon functions
isPFC = @(x) ~isempty( strfind(x,'K081'));
isAC  = @(x) ~isempty( strfind(x,'K080'));

% For each ferret
for i = 1 : length(ferrets)
    
    % List blocks
    ferDir = fullfile( rootDir, ferrets(i).name);
    blocks = dir( fullfile( ferDir, 'Block*'));
    
    % For each block
    for j = 1 : numel(blocks)
        
        % Define files
        blockDir = fullfile( ferDir, blocks(j).name);
        saveName = fullfile( saveDir, sprintf('%s_%s', ferrets(i).name, blocks(j).name));
        
        if exist([saveName '.mat'],'file'), 
            fprintf('%s already extracted\n', blocks(j).name)
            continue
        else
            fprintf('%s in progress\n', blocks(j).name)
        end
        
        % Get TDT file        
        TDT_file = sprintf('%s_TDT_data.mat', blocks(j).name);
        
        load( fullfile( blockDir, TDT_file),'Sens','Valv') % DOut is struct: data = voltage trace, fs = sample rate
        
        % Downsample and plot
        plotSens = downsample( transpose(Sens.data),10);
        plotValv = downsample( transpose(Valv.data),10);
        
        valv.N = size(Valv.data,2);
        sens.N = size(Sens.data,2);
        
        valvTime = downsample(linspace(1,valv.N, valv.N) ./ Valv.fs, 10);
        sensTime = downsample(linspace(1,sens.N, sens.N) ./ Valv.fs, 10);
        
        f = figure('name', sprintf('%s_%s', ferrets(i).name, blocks(j).name));
        subplot(211)
        plot(sensTime, plotSens); title('Sensors')
        subplot(212)
        plot(valvTime, plotValv); title('Valves')
        xlabel('Time (s)')
        legend
        
        % Get onsets
        for k = 1 : size(Sens.data,1)            
            onsets.Sens{k} = detectOnsets(Sens.data(k,:)); % In samples
            onsets.Sens{k} = onsets.Sens{k} ./ Sens.fs;   % Convert samples to time
        end
                
        for k = 1 : size(Valv.data,1)            
            onsets.Valv{k} = detectOnsets(Valv.data(k,:)); % In samples
            onsets.Valv{k} = onsets.Valv{k} ./ Valv.fs;   % Convert samples to time
        end
                                            
        % Save figure
        save([saveName '.mat'],'onsets')
        saveas(f, [saveName '.fig'])
        close(f)
    end
end


function samples = detectOnsets(x)

% Normalize to one
x = x ./ max(x);

% Escape if no onsets
if all( isnan(x)), samples = []; return; end

% Check if data is binary (not sure why it wouldn't be but check anyway)
ux = unique(x);
ux = transpose(ux(:));

if ~ismember(ux,[0 1],'rows')
   warning('Data is not binary?!')
   keyboard
end

% Onsets defined as change from 0 to 1
samples = find(diff(x) == 1);





