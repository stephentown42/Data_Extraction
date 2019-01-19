function compare_digital_out_times
%
% Loads snippets from extracted broadband trace. There are n snippets,
% where n is the number of trials. Each snippet is m samples long, where m
% equals the snippet duration x sampling rate.
%
% Plots spectra for each channel together on one channel
% Plots spectrograms for each channel on separate axes in a second figure

try

% Add paths
addpath( genpath( 'C:\Users\Dumbo\Documents\MATLAB\lib\chronux_2_11')) 
saveDir = 'D:\Frontiers Data Analysis\Timing\Digital';

% List subjects
rootDir = 'C:\Data\UCL_Behaving';
ferrets = dir( fullfile(rootDir,'F*'));

% Anon functions
isPFC = @(x) ~isempty( strfind(x,'K082'));
isAC2 = @(x) ~isempty( strfind(x,'K046'));
isAC  = @(x) ~isempty( strfind(x,'K080'));

% For each ferret
for i = 1 : length(ferrets)
    
    % List blocks
    ferDir = fullfile( rootDir, ferrets(i).name);
    blocks = dir( fullfile( ferDir, 'Block*'));
    
    % For each block
    for j = 1 : numel(blocks)
        
        % Check if extraction has been completed
        saveName   = sprintf('%s_%s_*', ferrets(i).name, blocks(j).name); 
        savedFiles = dir(fullfile(saveDir, saveName));
            
        if numel(savedFiles) == 4
            fprintf('%s complete - skipping\n', saveName); continue
        else
            fprintf('%s in progress\n', saveName)
        end
        
        % Get TDT file        
        blockDir = fullfile( ferDir, blocks(j).name);
        TDT_file = sprintf('%s_TDT_data.mat', blocks(j).name);
        
        load( fullfile( blockDir, TDT_file),'DOut') % DOut is struct: data = voltage trace, fs = sample rate
        
        plotTDT = downsample(DOut.data, 1000);  % Downsample for quick plotting
                
        % Detect TDT onsets
        onsets.TDT = detectOnsets(DOut.data);        
        onsets.TDT_time = onsets.TDT ./ DOut.fs;
        
        % List recording data files
        recFiles = dir( fullfile( blockDir, '*McsRecording*.mat'));
        
        % For each recording file
        for k = 1 : numel(recFiles)
            
            % Reset
            axColor = 'w';
            
            % Infer location from file name 
            if isPFC(recFiles(k).name), str = 'PFC'; end
            if isAC(recFiles(k).name),  str = 'AC'; end
            if isAC2(recFiles(k).name),  str = 'AC'; end
            
            % Check if file exists            
            saveName = sprintf('%s_%s_%s.fig', ferrets(i).name, blocks(j).name, str); 
            savePath = fullfile(saveDir, saveName);
            
            if exist( savePath,'file')
                fprintf('%s exists - skipping\n', saveName); continue
            else
                fprintf('%s in progress\n', saveName)
            end
            
            % Load data
            recPath = fullfile( blockDir, recFiles(k).name);        
            load( recPath, 'Digital_1')
            
            % Convert to single date type
            Digital_1 = single(Digital_1.ChannelData);
            plotMCS   = downsample(Digital_1,100);
            
            % Create visualization
            f = figure('Name',recFiles(k).name,...
                       'units','normalized',...
                       'position',[0.01 0.05 0.98 0.87]);            
            
            subplot(2,4,[1 2 3])
            plot(plotMCS)
            xlabel('Sample')
            ylabel('Voltage')
            title('Wireless (downsampled x 100)')
            axis tight
            
            subplot(2,4,[5 6 7])
            plot(plotTDT)
            xlabel('Sample')
            ylabel('Voltage')
            title('TDT (downsampled x 1,000)')
            axis tight            
            
            % Detect wireless onsets
            onsets.MCS = detectOnsets(Digital_1);
            
            % If lots more onsets in wireless sequence, likely that the wireless
            % system was left recording        
            if numel(onsets.TDT)+10 < numel(onsets.MCS)                
                onsets.MCS = onsets.MCS(1:numel(onsets.TDT)); 
                axColor = [1 1 0];
            end
            
            % Single value shift usually indicates TDT lost first onset
            if numel(onsets.TDT) ~= numel(onsets.MCS)            
                onsets.MCS = onsets.MCS(2:end); %keyboard
            end 
            
            % Less clear cases
            if numel(onsets.TDT) ~= numel(onsets.MCS)
                                
                diffN = numel(onsets.MCS) - numel(onsets.TDT);
                
                warning('Wireless and TDT data do not match for this block (%d)', diffN)
                
                if diffN < 0 || diffN > 3
                    keyboard
                    continue                                        
                end
                                           
                for ni = 1 : diffN + 1                   
                    
                    WL  = numel(onsets.TDT)-1; % Window length
                    Xni = onsets.TDT;
                    Yni = onsets.MCS(ni:ni+WL);
                    
                    [Pni, Sni] = polyfit(Xni(:)', Yni(:)', 1);
                    pred = polyval(Pni, Xni(:)');   % Prediction based on fit
                    res = pred(:) - Yni;    % Residuals
                    totalRes(ni) = sum(abs(res));
                end
                
                idx = find(totalRes == min(totalRes));  % Find minimal sum of residuals
                onsets.MCS = onsets.MCS(idx:idx+WL);    
                axColor = [1 0.5 0];
            end
                        
            % Fit equation to get sample rate and offset
            P = polyfit(onsets.TDT_time(:), onsets.MCS(:), 1);
            
            % Save Data 
            matSaveFile = strrep( savePath,'.fig','.mat');
            save( matSaveFile, 'onsets','P')
                        
            
                        
            subplot(2,4,4)
            hold on
            scatter(onsets.TDT_time, onsets.MCS)
            plot( onsets.TDT_time, P(2)+P(1).*onsets.TDT_time,'-r')  
            set(getAxes,'color',axColor)
                                                            
            % Save figure            
            saveas(f, savePath)
            close(f)                        
        end
    end
end

catch err
    err
    keyboard
end


function samples = detectOnsets(x)

% Normalize to one
x = x ./ max(x);

% Check if data is binary (not sure why it wouldn't be but check anyway)
ux = unique(x);
ux = transpose(ux(:));

if ~ismember(ux,[0 1],'rows')
   warning('Data is not binary?!')
   keyboard
end

% Onsets defined as change from 0 to 1
samples = find(diff(x) == 1);





