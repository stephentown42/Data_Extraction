function contrastStartTimes_And_Dout


% Demo arguments
ferret  = 'F1703_Grainger';
tank    = fullfile('D:\UCL_Behaving\', ferret);
block   = 'Block_J5-120';
bPath   = fullfile('D:\Behavior', ferret);
bFile   = '24_04_2018 level53_Granger 08_13_42.910 Block_J5-120.txt';

% Load behavioral file
B = importdata( fullfile( bPath, bFile));

% Convert behavioral data to table
rmSpace = @(x) strrep(x,' ','');
myHeaders = cellfun( rmSpace, B.colheaders,'un',0);
B = array2table(B.data,'VariableNames',myHeaders);
nTrials = size(B,1);

% Load tdt data
T = TDT2mat(tank, block);

% Create time vector for TDT
S = T.streams.DOut;
tvec = 1 : numel(S.data);
tvec = tvec ./ S.fs;

% Compare 
figure
hold on
plot(tvec, S.data(1,:))
plot(B.StartTime, zeros(1,nTrials)+3,'dk','MarkerFaceColor','k')


% Get window
windowTime = [-1 1];
windowSamps = round(windowTime .* S.fs);
windowIdx = windowSamps(1):windowSamps(2);
B.startSamp = round(B.StartTime.*S.fs); 
windowIdx = bsxfun(@plus, B.startSamp, windowIdx);
startWindow = S.data(windowIdx);

% Draw window
window_tvec = linspace(windowTime(1), windowTime(2), numel(windowSamps));
figure
hold on
imagesc(window_tvec, B.Trial, startWindow)

% Add graphical representation of key events
duration = 0.25;
isi = 0.05;
set(plotYLine(0),'color','w','LineWidth',1) % Start time
set(plotYLine(-isi),'color','w','LineWidth',1) % ISI
set(plotYLine(-isi-duration),'color','w','LineWidth',1) % Stimulus + ISI

% Add trial labels of interest
CorrectionTrials = B.Trial(B.CorrectionTrial == 1);
scatter(repmat(-0.15, size(CorrectionTrials)), CorrectionTrials,'.y')

CorrectTrials = B.Trial(B.Correct == 1);
scatter(repmat(-0.2, size(CorrectTrials)), CorrectTrials,'xg')

CorrectTrials = B(B.Correct == 1,:);



% Get stimulus drive in key window
keyWindow_time = [-duration 0] - isi;
keyWindow_Samps = round(keyWindow_time .* S.fs);
keyWindowIdx = keyWindow_Samps (1):keyWindow_Samps (2);
keyWindowIdx = bsxfun(@plus, B.startSamp, keyWindowIdx);
keyWindow = S.data(keyWindowIdx);
stimDrive = sum(keyWindow,2);

% Filter for ok data
K = B(stimDrive > 0,:);
K = K(K.CorrectionTrial == 0,:); % Remove correction trials





