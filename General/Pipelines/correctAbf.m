%% correctAbf
% Corrects abf file. Loads and calculates abf file from md voltage
% recording file.


clear all; clc

% Establish inputs
mouseID = 'TM230304-2\';
baseFolder = 'D:\Mice\';

% Select mouse raw data folder
cd([baseFolder mouseID])

top = pwd;

% Find all days
d = dir();
isub = [d(:).isdir]; %# returns logical vector
nameDays = {d(isub).name}';
nameDays(ismember(nameDays,{'.','..'})) = [];

% Find all locs
allLocs = {};

for ff = 1:length(nameDays)
    cd(nameDays{ff})

    d2 = dir('loc*');
    isub = [d2(:).isdir]; %# returns logical vector
    subLocs = {d2(isub).name}';

    for gg = 1:length(subLocs)
        allLocs{end+1} = [nameDays{ff} '\' subLocs{gg}];
    end

    cd(top)
end


%% Run correction

for ff = 1:length(allLocs)
    cd([baseFolder mouseID allLocs{ff}])

    if ~exist([pwd '\claimed.mat'],'file')
        continue
    end

    % folder where pipeline code is saved
    copyFolder = 'D:\Matlab\suite2p\';

    % copy self to current directory
    if ~exist([pwd '\correctAbf.m'],'file')
        copyfile([copyFolder 'correctAbf.m'],'correctAbf.m');
    end

    % identify folder with original tiffs
    f='TSeries*';
    d = dir(f);
    folderPath = d(1).name;

    % generate Z folder names
    outFolders = {[folderPath '_1'],[folderPath '_2']};

    dwnSample = 2;      % downsampling factor


    %% Generate abf

    % cycle through Z planes
    for ii = 1:2
        if ~exist([pwd '\' outFolders{ii}],'dir')
            continue
        end

        cd(outFolders{ii})

        %% Generate abf

        % nFrameDetected is the number of volatge trace time points
        load('md.mat')

        abf = struct();

        % update md with parameters
        fRate = mean(diff(md(:,1)))/1000;   % frame rate (sec)
        save('md.mat','md','fRate','dwnSample','-v7.3')

        nFrameDetected = size(md,1);

        % nOriginalFrames is the number of imaging frames after motion correction
        load('motionCorrected5times/nOriginalFrames.mat');

        if nOriginalFrames>nFrameDetected
            nOriginalFrames=nFrameDetected;
        end

        abf.imageIndex = (1:1:nOriginalFrames)';
        abf.t = fRate*(abf.imageIndex-1);

        abf.y = zeros(nOriginalFrames,1);
        for jj = 1:nOriginalFrames
            abf.y(jj) = mean(md(jj:jj+dwnSample-1,3)*100);
        end


        %% delete old analysis

        cd('suite2p')

        d2 = dir();

        isDir = cat(1,d2.isdir);
        d3 = {d2.name}';
        d3 = d3(~isDir);

        n = find(strcmp(d3,'Fall.mat'));
        if ~isempty(n)
            d3{n} = [];
        end

        for k=1:length(d3)
            delete(d3{k})
        end


        %% Save updated abf

        save('abf.mat','abf');

        abfFake = abf;
        save('abfFake.mat','abfFake')

        cd ..\..
    end


    disp([baseFolder mouseID allLocs{ff}])
end

cd(baseFolder)

