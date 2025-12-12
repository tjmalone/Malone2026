%% countCellsPerArea
% Counts number of cells
% Calculates nissl cells per area

clear; close all; clc

%%%% Input
versionName = 'FinalVersion';
foldBase = 'Z:\SNM\labMembers\KC\AD_Project\Histology\reelCalNeurodegen';
%%%%%

useT = 1;
useIdx = 1; %calbindin channel
types = {'reelin','calbindin'};
cNum = length(types);
images = 'MEC_brightCal.tif';
degenRegion = 'degen_brightCal_v4.roi';

% area per pixel (um^2/pixel)
APPum = (850.19/1024)^2; %Confocal Microscope
APP = (APPum)*(1/1000)^2; % in mm^2

%% find folder paths to mice for analysis
foldNeuroden = 'Z:\SNM\labMembers\KC\AD_Project\Histology\tau_RC\neurodegen';
cd(foldNeuroden)
load('gcampMiceFolds.mat')
fNum = length(folds);

% load cell numbers
cd([foldBase '\' versionName '\remBorderCells'])
load('removeBorderCells.mat')

% create save fold
cd([foldBase '\' versionName])
mkdir('countCellsPerArea')
cd('countCellsPerArea')
saveFold = pwd;

%% Loop through all mice

nRoi = zeros(fNum,cNum);
area = zeros(fNum,cNum);
roisPerArea = zeros(fNum,cNum);
cellFiles = cell(fNum,cNum);

% for c = 1 % reelin only
for c = 1:cNum

    % for f = 12
    for f = 1:fNum % loops through mice

        cd(folds{f})

        % Load Image
        for i = 1 % loop through 1 planes
            curIm = tiffreadVolume(images);
            
            % normalize image
            curIm = double(curIm);
            curIm = curIm(:,:,useIdx(useT));
             
            % Extract ROIs
            nRoi(f,c) = size(cellsInsideBorder{f,c},1);
    
            % Extract "degen_1st", the section of layer 2 used for the
            % analysis
            [rows,columns] = size(curIm,1:2);
             reg = extractImageJROI(rows,columns,degenRegion);       
    
            % calculate area
            rProps = regionprops(reg(:,:,1),'Area','PixelIdxList');
            area(f,c) = rProps.Area*APP; %% should be in mm^2
        end    
    end
end

roisPerArea = nRoi./area;

cd(saveFold)
    
save(['countCellsPerArea.mat'],'roisPerArea','nRoi','area','folds','types');

% folder where code is saved
baseFolder = 'D:\MATLAB\Code\AnalysisCode\Analysis\HistologyCode\ReelCalNeurodegen\withRemoveBorderCells';

% copy self to current directory
copyfile([baseFolder '\countCellsPerArea_remBord.m'],'countCellsPerArea_remBord.m');