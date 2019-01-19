function h5_to_Mat_batch

% Extract all data from the wireless recordings obtained with the Phoenix
% Experimental Setup on Jumbo
%

% List ferrets
% rootDir = 'C:\Data\UCL_Behaving\';
rootDir  = 'C:\Users\Dumbo\Documents\Multi Channel DataManager';
% savePath = 'D:\h5_to_Mat';
% ferrets = dir( fullfile( rootDir, 'F*'));

% % For each ferret
% for i = 1 : length(ferrets)
% 
%     % Get block list
%     tank   = fullfile( rootDir, ferrets(i).name);
%     blocks = dir( fullfile(tank, 'Block_J2-*'));

%     % For each block
%     for j = 1 : numel(blocks)
% 
%         % Request h5 file
%         blockPath = fullfile(tank, blocks(j).name);
%         h5file = dir( fullfile(blockPath, '*.h5'));
        h5file = dir( fullfile(rootDir, '2017*.h5'));

%         % Skip if no file
%         if isempty(h5file)
%             continue
%         end

        % Otherwise run main function
        for k = 1 : numel(h5file)

%             main( blockPath, h5file(k).name)
            main(rootDir, h5file(k).name)
        end
%     end
% end



% function main(loadPath, loadName)
function main(loadPath, loadName)

try
    
    %Define data location
    loadFile  = fullfile(loadPath, loadName);
    headstage = str2num(loadName(end-5:end-3));

    % Define save file name    
    saveName = strrep(loadName,'.h5','.mat');
    saveFile = fullfile(loadPath, saveName);   
    
      
    % Skip if file already exists
    if exist(saveFile,'file')
        fprintf('%s already exists - skipping\n', saveName); return
    else
        fprintf('%s - extracting\n', saveName)
    end
    
    % Check file structure
    info = h5info(loadFile);
%     h5disp(loadFile);
%     
%     % Append metadata
%     M.date = info.Groups(1).Attributes(6).Value; % Day
%     M.nSamps = info.Groups.Groups.Attributes(4).Value; % Duration in samples
        
    % Define file structure
    trunk    = '/Data/Recording_0/';
    streams  = {'Analog','Event'}; %,'Segment'};% ,'TimeStamp'};
    
    fieldStr.Analog    = {'ChannelData','ChannelDataTimeStamps','InfoChannel'};
    fieldStr.Event     = {'EventEntity_0','InfoEvent'};
    fieldStr.Segment   = {'InfoSegment','SegmentData_*','SegmentData_ts_*','SourceInfoChannel'};
    fieldStr.TimeStamp = {'InfoTimeStamp','TimeStampEntity_*'};
    
    nStreams = [4         1       1         1];
                     
    label.Analog    = {'Analog_1','Filter_1','Raw_Electrode_1','Digital'};
    label.Event     = {'Digital_Events_1'};
%     label.Segment   = {'Spike_Data_1'};
%     label.TimeStamp = {'Spike_Timestamps_1'};
            
    % For each stream type
    for i = 1 : numel(streams)
        
%         warning('Processing analog streams only')

        % Skip segment data when looking at files from PFC (they're only in
        % the files from 32 channel headstages)
        if i == 3 && ~isempty(headstage)
            if headstage == 67 || headstage == 82
                continue
            end
        end
        
        % Define branch
        iBranch = [trunk streams{i} 'Stream/'];
                 
        % Get field names for this type of data
        eval( sprintf('jNames = fieldStr.%s;', streams{i}))

        % For each stream within stream type
        for k = 1 : nStreams(i)
           
            % Get output field
            eval( sprintf('outputField = label.%s{k};', streams{i}))
            
            % For each field
            for j = 1 : numel(jNames)
                                
                                
                % Update user about extraction progress
                fprintf('\tExtracting %s...\n', outputField);
                
                % If all channels stored together
                if isempty( strcmp(jNames{j},'*'))
                    
                    % Define next branch
                    kBranch = [iBranch sprintf('Stream_%d/', k-1) jNames{j}];     
                    
                    % Read data to output
                    eval(sprintf('MCS.%s.%s = h5read(loadFile, kBranch);', outputField, jNames{j}))  
                    
                % If channels have separate fields (i.e. spike data)
                else                    
                    % For each channel
                    for chan = 1 : 16
                                                
                        % Define next branch                
                        chanName = strrep(jNames{j},'*',num2str(chan-1));
                        kBranch  = [iBranch sprintf('Stream_%d/', k-1) chanName];  
                                                
                        % Read data to output
                        eval(sprintf('MCS.%s.%s = h5read(loadFile, kBranch);', outputField, chanName)) 
                    end                    
                end                
            end
        end
    end
    
    % Save data
    save( saveFile, '-struct','MCS','-v7.3')    
    
catch err
    err
    keyboard
end

