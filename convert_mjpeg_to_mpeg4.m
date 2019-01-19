function convert_mjpeg_to_mpeg4(filePath, fileName)
%
% Note: From the demo quality, this isn't really very good. Better to do
% the encoding correctly during recording or live with the excessively
% large file sizes in the past.
%
% ST: 06 April 2018


% Demo inputs
filePath = 'D:\Matlab_Videos';
fileName = '2018-04-06_Track_08-32-42.avi';

% Open old object
oldObj = VideoReader( fullfile( filePath, fileName));

% Figure out number of frames
nFrames = ceil(oldObj.Duration * oldObj.FrameRate);
frameCount = 0; % Initialize counter

% Change name
newFile = strrep(fileName,'.avi','_mpg4.avi');

% Create new object
newObj =  VideoWriter( fullfile(filePath, newFile),'MPEG-4');
newObj.FrameRate = 10;
open(newObj)

% Set up waitbar
str = 'Converting to MPEG4';
h = waitbar(0, str);

% For each frame
while hasFrame(oldObj)
    
    % Update user
    frameCount = frameCount + 1;
    myProgress = frameCount / nFrames;
    waitbar(myProgress, h, sprintf('%s: %.2f%%', str, myProgress*100))
    
    % Read frame and write to new video    
    writeVideo(newObj, readFrame(oldObj));
end

% Close objects
 close(oldObj) 
 close(newObj)
 close(h)