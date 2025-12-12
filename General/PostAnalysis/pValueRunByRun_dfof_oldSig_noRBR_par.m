function pValueRunByRun_dfof_oldSig_noRBR_par(trackEnd,binWidth)
%% pValueRunByRun_sig_oldSig_noRBR_par
%
% Define grid cells using calcium transients (dfof). Uses
% parallelization for implementation. Run-by-run consistency must be
% previously calculated.
%
% Inputs:
%   trackEnd - length of track (in cm)
%   binWidth - track bin size (in cm)
%


%% Load parameters

% load dfof_sig
name = 'dfof_*';
fn = dir(name);
load(fn(1).name,'dfof');
load(fn(2).name,'dfof_sig');

% load interpolated dfof_sig
load('RunByRun_dfof\dfofM.mat','dfofM')
N = length(dfofM);
useRunNumbers = 1:1:size(dfofM{1},1);

% load abfFake
load('abfFake.mat','abfFake');
abf = abfFake;

% calculate speed threshold
[speedThreshold] = speedThreshold1D(1,length(abfFake.t),abfFake,0);

% set parameters
Pvalue1=0.8;
Pvalue2=0.8;
Pvalue3=0.25;

trackStart=0;
nCellStart = 1;
nCellEnd = N;
imageStartNumber = 1;
imageEndNumber = length(abfFake.t);
abfStartNumber = imageStartNumber;
abfEndNumber = imageEndNumber;


%% Generate output structures

mkdir('gridAnalysis_dfof');
mkdir('gridAnalysis_dfof\grid');
mkdir('gridAnalysis_dfof\nongrid');

% save parameter information
params = cell(1,N);
for n=1:N
    params{n}.useRunNumbers=useRunNumbers;
    params{n}.nCellStart=nCellStart;
    params{n}.nCellEnd=nCellEnd;
    params{n}.imageStartNumber=imageStartNumber;
    params{n}.imageEndNumber=imageEndNumber;
    params{n}.abfStartNumber=abfStartNumber;
    params{n}.abfEndNumber=abfEndNumber;
    params{n}.trackStart=trackStart;
    params{n}.trackEnd=trackEnd;
    params{n}.binWidth=binWidth;
    params{n}.speedThreshold=speedThreshold;
    params{n}.Pvalue1=Pvalue1;
    params{n}.Pvalue2=Pvalue2;
    params{n}.Pvalue3=Pvalue3;
end

% set field names and data types
fieldNames = {'dfofUse','dfofaveragesmooth','dfofaveragesmoothFields',...
    'Pvalues','inFieldBins','outFieldBins','transitions','percentageBins',...
    'fieldStartEnds','fieldCenters','fieldSpacings','fieldWidths',...
    'meanInField','meanOutField','ratioInOut','minSpacing','maxWidth','indices'};
fieldInit = {[],[],[],...
    [],{},{},[],[],...
    {},{},{},{},...
    [],[],[],[],[],[]};

% initialize output structs
grid = struct();
nongrid = struct();
allCells = struct();

% initialize struct fields
for ff = 1:length(fieldNames)
    grid.(fieldNames{ff}) = fieldInit{ff};
    nongrid.(fieldNames{ff}) = fieldInit{ff};
    allCells.(fieldNames{ff}) = fieldInit{ff};
end


%% Run grid cell classifier

disp('Running grid cell classifier');

isGrid = zeros(1,length(nCellStart:nCellEnd));
fieldsAll = cell(1,length(nCellStart:nCellEnd));

% initialize parallel pool if none exists
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool;
end

% save command window text with classification reasons
diary gridAnalysis_dfof\output.out

% run grid cell tests in parallel
parfor n=1:length(nCellStart:nCellEnd)
    cellN = nCellStart+n-1;
    
    [answer,fields]=isGridCellKYModified2_includeSig(cellN,...
        imageStartNumber,imageEndNumber,abfStartNumber,abfEndNumber,...
        trackStart,trackEnd,binWidth,speedThreshold,...
        Pvalue1,Pvalue2,Pvalue3,dfof,dfof_sig,abf);
    
    fields.minSpacing = min(fields.fieldSpacings);
    fields.maxWidth = max(fields.fieldWidths);
    
    isGrid(n) = answer;
    fieldsAll{n} = fields;
    
    if answer
        savedir ='gridAnalysis_dfof\grid';
    else
        savedir ='gridAnalysis_dfof\nongrid';
    end
    
    %save the figure
    filename1 = sprintf('%s_%d.fig', 'cell', cellN);
    savefig(fullfile(savedir,filename1))
%     filename2 = sprintf('%s_%d.tif', 'cell', cellN);
%     saveas(gcf,fullfile(savedir,filename2));
end

%stop saving command
diary off


%% Store field information

% store field information
for n=1:length(nCellStart:nCellEnd)
    fields = fieldsAll{n};
    
    for ff = 1:length(fieldNames)-1
        if iscell(grid.(fieldNames{ff}))
            if isGrid(n)
                grid.(fieldNames{ff}){end+1} = fields.(fieldNames{ff});
            else
                nongrid.(fieldNames{ff}){end+1} = fields.(fieldNames{ff});
            end
            allCells.(fieldNames{ff}){n} = fields.(fieldNames{ff});
        else
            if isGrid(n)
                grid.(fieldNames{ff})(:,end+1) = fields.(fieldNames{ff});
            else
                nongrid.(fieldNames{ff})(:,end+1) = fields.(fieldNames{ff});
            end
            allCells.(fieldNames{ff})(:,n) = fields.(fieldNames{ff});
        end
    end
end

%find indices for grid and non grid cells
grid.indices = (find(isGrid)+nCellStart-1)';
nongrid.indices = (find(~isGrid)+nCellStart-1)';
allCells.indices = (nCellStart:nCellEnd)';

%save results:
save('gridAnalysis_dfof\grid.mat','grid');
save('gridAnalysis_dfof\nongrid.mat','nongrid');
save('gridAnalysis_dfof\params.mat','params');
save('gridAnalysis_dfof\allCells.mat','allCells');


%% fix grid fields
clear allCellsCorrected
cd('gridAnalysis_dfof');
plotstuff=0;
trackLength=trackEnd-trackStart;
correctPValue(allCells,grid,nongrid,params,binWidth,trackLength,dfof,abf,plotstuff);
cd ..

end

