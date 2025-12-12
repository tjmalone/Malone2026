%% stack2tif
% converts all specified tif stacks to average projection tifs and saves in
% new location.


%% Identify folders

stackF = 'ref/refStackMCFinal*.tif';
saveF = 'mcStacks';
baseF = 'F:\Mice\MotionCorrected\ID_210413\';
subF = 'motionCorrected5times';
subN = 2;
idxs = 0;

files = findSubF(subF,subN,baseF,idxs);
filesN = length(files);


%% Generate projection

for f = 1:filesN
    d = dir([files{f} stackF]);
    imStack = loadtiff(fullfile(d(1).folder,d(1).name));
    
    imAvg = mean(imStack,3);
    
    fUnq = erase(files{f},baseF);
    
    saveastiff(int16(imAvg),fullfile(saveF,fUnq,[saveF '.tif']));
end