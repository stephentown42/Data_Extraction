function verifyTrialOnset_CF(ferret)


if nargin == 0    
    ferret = 'F1807*';
end

% Specify paths
dirs.root = 'D:\UCL_Behaving';
dirs.behav = 'D:\Behavior';


% List subjects
ferrets = dir( fullfile(dirs.root, ferret));

% For each ferret
for i = 1 : length(ferrets)
    
    % Extend path
    dirs.tank = fullfile( dirs.root, ferrets(i).name);
    dirs.behav_ferret = fullfile( dirs.behav, ferrets(i).name);    
    
    % List blocks
    blocks = dir( fullfile( dirs.tank, 'Block*'));
    
    % Filter for blocks today
    block_datenum = cat(1, blocks.datenum);
    idx = now - block_datenum < 1/2;
    blocks = blocks(idx);
    
    % For each block
    for j = 1 : numel(blocks)
               
        % Extend path
        dirs.block = fullfile( dirs.tank, blocks(j).name);
        
        % Load TDT data
        DOut = TDT2mat(dirs.tank, blocks(j).name, 'STORE','DOut','Verbose',0);
        Sens = TDT2mat(dirs.tank, blocks(j).name, 'STORE','Sens','Verbose',0);
        
        % Crop strucutre
        DOut = DOut.streams.DOut;
        Sens = Sens.streams.Sens;
        
        % Draw distribution
        pxl_edges = 0 : 5 : 255;
        pxl_inst = histc(DOut.data(4,:), pxl_edges);
        
        figure
        hold on
        bar(pxl_edges, pxl_inst, 'histc')
        xlabel('Pixel Value')
        ylabel('Instances (n)')
        
        
        % Create time vectors
        DOut.tvec = [1 : size(DOut.data, 2)] ./ DOut.fs;
        Sens.tvec = [1 : size(Sens.data, 2)] ./ Sens.fs;
        
        % Extract center trace from sensor record
        centerSpout = Sens.data == 64;
        
        % Create figure
        figureST( sprintf('%s: %s', ferrets(i).name, blocks(j).name));
        hold on
        plot(Sens.tvec, centerSpout, 'b')
        plot(DOut.tvec, DOut.data(1,:), 'r') % Stimulus play        
        xlabel('Time (s)')
        
        % Find behavioral file
        behavior_file = dir( fullfile( dirs.behav_ferret, ['*' blocks(j).name '.txt']));
        
        if numel(behavior_file) ~= 1           
            warning('Behavioral file number not equal to 1'); keyboard;
        end
        
        % Import behavioral file
        B = importdata( fullfile( dirs.behav_ferret, behavior_file(1).name));
        
        % Mark time of trial onset
        set(gca,'ylim',[0 2.5])
        
        startTimes = B.data(:, strcmp(B.colheaders, 'StartTime'));
        
        for start_time_idx = 1 : numel(startTimes)            
           plotYLine(startTimes(start_time_idx));
        end
        
    end
end




