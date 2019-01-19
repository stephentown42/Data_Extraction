function reduceVideoQuality

% Options
targetQuality   = 50;
frameRate       = 10;

% List files
dirs.root = 'D:\Matlab_Videos';
files = dir( fullfile( dirs.root, '*.avi'));

% For each file
for i = 1 : numel( files)
    
   % Load file
   vidPath = fullfile( dirs.root, files(i).name);
   vidObj = vision.VideoFileReader(vidPath);
    
   % Create new file
   newFile = strrep(files(i).name,'Track','Edit');
   newObj  = vision.VideoFileWriter( fullfile( dirs.root, newFile));
   
   newObj.FrameRate = frameRate;
   newObj.VideoCompressor = 'MJPEG Compressor';
   newObj.Quality = 50;
   
end