function get_TDT_data

% Paths
rootDir = 'C:\Data\UCL_Behaving';
ferrets = dir( fullfile(rootDir,'F*'));

% Run program for each ferret
for i = 1 : length(ferrets)
    
    % Define paths
    tank = fullfile( rootDir, ferrets(i).name); 
        
    % List blocks
    blocks = dir( fullfile(tank,'Block_J*'));    
        
    % For each file
    for j = 1 : length(blocks)
        
        % Check save file 
        saveName = sprintf('%s_TDT_data.mat',blocks(j).name);
        savePath = fullfile(tank, blocks(j).name, saveName);
        
        if exist(savePath,'file')
            fprintf('%s already extracted - skipping\n', saveName)
            continue 
        else
            fprintf('Extracting %s...\n', saveName)
        end
        
        % Extract
        try
            S = TDT2mat(tank, blocks(j).name);  

            % Format to include metadata
            S.streams.metadata = S.info;
            S.streams.tick = S.epocs.Tick;
            S = S.streams;

            % Save
            save( savePath,'-struct','S'); 
        catch err
            err
            keyboard
        end
    end
end

