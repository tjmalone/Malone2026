function pValueRunByRun_sig_subset(trackEnd,binWidth,speedThreshold,dfof,dfof_sig,abf,useRuns)
%% pValueRunByRun_subset
%
% Identify fields using a subset of runs. Uses parallelization for
% implementation. dfofM_sig  must be previously calculated.
%
% Inputs:
%   trackEnd - length of track (in cm)
%   binWidth - track bin size (in cm)
%


%% Load parameters

N = size(dfof,2);

% set parameters
Pvalue1=0.8;
Pvalue2=0.8;
Pvalue3=0.25;

trackStart = 0;
nCellStart = 1;
nCellEnd = N;
imageStartNumber = 1;
imageEndNumber = length(abf.t);
abfStartNumber = imageStartNumber;
abfEndNumber = imageEndNumber;


%% Generate output structures

% save parameter information
params = struct();
params.useRunNumbers=useRuns;
params.nCellStart=nCellStart;
params.nCellEnd=nCellEnd;
params.imageStartNumber=imageStartNumber;
params.imageEndNumber=imageEndNumber;
params.abfStartNumber=abfStartNumber;
params.abfEndNumber=abfEndNumber;
params.trackStart=trackStart;
params.trackEnd=trackEnd;
params.binWidth=binWidth;
params.speedThreshold=speedThreshold;
params.Pvalue1=Pvalue1;
params.Pvalue2=Pvalue2;
params.Pvalue3=Pvalue3;


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
allCells = struct();

% initialize struct fields
for ff = 1:length(fieldNames)
    allCells.(fieldNames{ff}) = fieldInit{ff};
end


%% Run grid cell classifier

svName = 'figures';
mkdir(svName)
disp('Running field identification');

fieldsAll = cell(1,length(nCellStart:nCellEnd));

% initialize parallel pool if none exists
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool;
end

% run grid cell tests in parallel
parfor n=1:length(nCellStart:nCellEnd)
    disp(n)
    cellN = nCellStart+n-1;
    
    fields = shuffle1D_dfof_includeSig_subset(cellN,imageStartNumber,imageEndNumber,abfStartNumber,abfEndNumber,...
        trackStart,trackEnd,binWidth,speedThreshold,Pvalue1,Pvalue2,Pvalue3,dfof_sig,dfof_sig,abf);

    fields.minSpacing = min(fields.fieldSpacings);
    fields.maxWidth = max(fields.fieldWidths);
    
    fieldsAll{n} = fields;
    
    %save the figure
    filename1 = sprintf('%s_%d.fig', 'cell', cellN);
    saveas(gcf,fullfile(svName,filename1));
end


%% Store field information

% store field information
for n=1:length(nCellStart:nCellEnd)
    fields = fieldsAll{n};
    
    for ff = 1:length(fieldNames)-1
        if iscell(allCells.(fieldNames{ff}))
            allCells.(fieldNames{ff}){n} = fields.(fieldNames{ff});
        else
            allCells.(fieldNames{ff})(:,n) = fields.(fieldNames{ff});
        end
    end
end

%find indices for grid and non grid cells
allCells.indices = (nCellStart:nCellEnd)';

%save results:
save('params.mat','params');
save('allCells.mat','allCells');


%% fix grid fields
clear allCellsCorrected
trackLength=trackEnd-trackStart;
correctPValue_clean(allCells,binWidth,trackLength);
close

end

