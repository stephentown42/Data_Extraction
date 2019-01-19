function reportTodaysPerformance(flag)

try
    
    % Define paths
    dirs.behavior = 'D:\Behavior';
    
    % Define flag if undefined
    if ~exist('flag','var'), flag = 0; end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % String mode: Report performance for last week tested of specific ferret
    if ischar(flag)
        
        % Extend paths
        ferrets.name = flag;
        dirs.ferret.behavior = fullfile( dirs.behavior, flag);
        
        % Find offset for day since last monday (use day number, 1 = sunday, 7 = saturday)
        this_day = weekday(now);
        offset_to_monday = 2 - this_day;
        if offset_to_monday >= 0, offset_to_monday = offset_to_monday - 7; end % Wrap for Sunday / current monday
        
        % Define array of date numbers for desired days
        currentTime = now + offset_to_monday + [0:4];
        
        % Check if not tested last week        
        str = [datestr(currentTime(1), 'dd_mm_yyyy') '*'];                        
        nFiles = nDir(fullfile( dirs.ferret.behavior, str)); 
        
        while nFiles == 0 && abs(now-currentTime(1)) < 31
            
            currentTime = currentTime - 7;
            
            str = [datestr(currentTime(1), 'dd_mm_yyyy') '*'];
            nFiles = nDir(fullfile( dirs.ferret.behavior, str));                        
        end
        
        
        % For each day of the week
        for i = 1 : 5
            
            % Note today's date
            str = [datestr(currentTime(i), 'dd_mm_yyyy') '*'];                        
            
            % Run main function
            main(str, ferrets, dirs)     
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Numeric mode: Report performance on last n (flag) days for all animals
        % tested on that day
        
    else
        
        % Define ferrets of interest
        ferrets = dir( fullfile(dirs.behavior, 'F1*'));
        
        % Note today's date
        currentTime = now + flag;
        str = [datestr(currentTime, 'dd_mm_yyyy') '*'];
        
        
        % For each ferret
        for i = 1 : numel(ferrets)
            
            % Extend paths
            dirs.ferret.behavior = fullfile( dirs.behavior, ferrets(i).name);            
            
            % Run main function
            main(str, ferrets(i), dirs)            
        end
    end
    
catch err
    err
    keyboard
end


function main(str, ferret, dirs)

% List behavioral files for today
files = dir(fullfile( dirs.ferret.behavior, str));


% For each file
for j = 1 : numel(files)
    
    
    % Report subject and time
    fprintf('%s\t', datestr(files(j).datenum,'HH:MM'))
    fprintf('%s\t', ferret.name)
    
    % Load data
    srcFile = fullfile(dirs.ferret.behavior, files(j).name);
    B = importdata(srcFile);
    
    % Catch formatting problems
    if ~isstruct(B), fprintf('No Data\n'); continue; end
    
    % Get performance variable
    Correct  = B.data(:,strcmp(B.colheaders,'Correct'));
    corTrial = B.data(:,strcmp(B.colheaders,'CorrectionTrial'));
    Response = B.data(:,strcmp(B.colheaders,'Response'));
    Spkr_loc = B.data(:,strcmp(B.colheaders,'Speaker Location'));
    
    % Get number of trials
    nTrials = numel(Correct);
    
    % Remove correction trials
    Correct = Correct(corTrial == 0);
    Response = Response(corTrial == 0);
    Spkr_loc = Spkr_loc(corTrial == 0);
    
    % Remove probe trials
    idx = Spkr_loc == 6 | Spkr_loc == 12;
    
    if any(idx)
        Correct = Correct(idx);
        Response = Response(idx);
    end
    
    % Get performance
    nCorrect = sum(Correct);
    nTest    = numel(Correct);
    pCorrect = mean(Correct) * 100;
    myBias   = mean(Response) - 6;
    
    % Report
    fprintf('%d Trials\t', nTrials)
    fprintf('%.1f %% Correct\t', pCorrect)
    fprintf('(%d / %d)\t', nCorrect, nTest)
    fprintf('Bias: %.3f\n', myBias)
    
end