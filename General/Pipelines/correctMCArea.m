function stackCorrected = correctMCArea()
%% correctMCArea


%% Get file info

% list all correctino info files
infoPath = [pwd '\correctionInfo\'];
fnInfo=dir([infoPath 'correctionInfo_*']);
infoNumber=length(fnInfo);
infoAllPaths=cellfun(@(x)([infoPath x]),{fnInfo(1:infoNumber).name},'uniformoutput',0);

% list all motion corrected tif files
matrixPath = [pwd '\'];
fnMatrix =dir([matrixPath '*tif']);
matrixNumber=length(fnMatrix);
matrixAllPaths=cellfun(@(x)([matrixPath x]),{fnMatrix(1:matrixNumber).name},'uniformoutput',0);

if infoNumber~=matrixNumber
    error('Motion correction was unsuccessful')
end

stackNumber = infoNumber;


%% Identify maximal shifts in all directions

options.overwrite = true;
xP = 1;
xN = 2;
yP = 3;
yN = 4;

maxShift = zeros(1,4);

% load each info to determine the shift across all stacks,
for n = 1:stackNumber
    load(infoAllPaths{n});
    
    maxShift(xP) = max([maxShift(xP) xx']);
    maxShift(xN) = min([maxShift(xN) xx']);
    maxShift(yP) = max([maxShift(yP) yy']);
    maxShift(yN) = min([maxShift(yN) yy']);
end

% calculate lost pixels
absShift = ceil(abs(maxShift));


%% Crop matrices

try
    load('stackCorrected.mat')
catch
    stackCorrected = zeros(stackNumber,1);
end

for n = 1:stackNumber
    if stackCorrected(n)==1; continue; end
    
    curStack = loadtiff(matrixAllPaths{n});
    
    newStack = curStack(absShift(yP)+1:end-absShift(yN),absShift(xP)+1:end-absShift(xN),:);
    
    saveastiff(newStack,matrixAllPaths{n},options);

    stackCorrected(n) = 1;
    save('stackCorrected.mat','stackCorrected')
end

