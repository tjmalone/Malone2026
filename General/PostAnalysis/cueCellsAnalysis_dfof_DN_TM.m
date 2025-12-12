function cueCellsAnalysis_dfof_DN_TM(~,binWidth)

%only calcualte left and right cells
% screen cells: at least one in field one out field
disp('screen cells')

load('allROIs.mat');
load('abfFake.mat');
name='dfof*.mat';
file=dir(name);

for n=1:length(file)
    load(file(n).name);
end

load('gridAnalysis_dfof\allCells.mat');

% mkdir('cueAnalysisNew_sig');
%calculate COMs
COMs=[]; %first column is X second is Y

for n=1:size(roi,3)
    thisroi=roi(:,:,n);
    [r,c]=find(thisroi);
    if ~isempty(c)
        X=mean([min(c) max(c)]);
    else
        X=0;
    end
    if ~isempty(r)
    Y=mean([min(r) max(r)]);
    else
        Y=0;
    end
    COMs(n,1)=X;    %#okay
    COMs(n,2)=Y;    %#okay
end

cd('cueAnalysis_dfof');

%screen the cells
% (1) at least one in field 

%has in fields
C1=zeros(size(allCells.indices,1),1);
%has out fields
C2=ones(size(allCells.indices,1),1);
for n=1:size(allCells.indices,1)
    if ~isnan(allCells.inFieldBins{n})
        C1(n)=1;
    end
    
    if ~isnan(allCells.outFieldBins{n})
        C2(n)=1;
    end
        
end

allIdx = [1:1:size(allCells.indices,1)]';   %#okay
C = C1.*C2;
useIdx=allIdx(C==1);

save('useIdx.mat','useIdx');

useDfof=dfof(:,useIdx);
useCOMs=COMs(useIdx,:);
useDfofaverage=dfofaveragesmooth(:,useIdx);
save('useDfof.mat','useDfof');
save('useCOMs.mat','useCOMs');
save('useDfofaverage.mat','useDfofaverage');

% 
% d = dir('temp*.mat');
% folderNames = [];
% templates = [];
% for n = 1:numel(d);
%     iii = strfind(d(n).name, '.');
%     
%     F = load(d(n).name);
%     templates(:, end+1) = F.(d(n).name(1:iii(end)-1));  %#okay
%     
%     folderNames{n} = d(n).name(1:iii(end)-1);           %#okay
% end

% Modified to specifically analyze left and right cues
load('tempL')
load('tempR')
templates = [tempL,tempR,min(tempR+tempL,1)];
folderNames = {'Left','Right','All'};

 
%% loop through cells
%this is the new code that makes shuffles using shuffling cue template
disp('get basic info');
nShuffle = 200;

mkdir('newScoreShuffleTemplate');
cd('newScoreShuffleTemplate');

indicesRange = 1:size(templates,1);
for ii = 1:numel(folderNames)
    disp(['Calculating ' folderNames{ii} ' cue scores']);
    mkdir(folderNames{ii});
    cd(folderNames{ii});
    clear cueCells
    temp = templates(:, ii);
    
    cueCells.useDfofaverage=useDfofaverage;
    cueCells.useCOMs=useCOMs;
    cueCells.useIdx=useIdx;
    cueCells.cueTemplate = temp;
    cueCells.maxLag = 30;

    cellNumber = 1;
    allCells.dfofaveragesmooth = temp;
    allCells.dfofaveragesmoothFields = temp;
    [cueCells.shuffleTemplates, ~] = extractFrateRandom_YGPValueFields(cellNumber, allCells, indicesRange, nShuffle, binWidth);
    
    lagBins = 30; %under 5cm bin, this is 150 shift
    cueCells.shuffleScores=[];
    for n=1:length(useIdx)
        for m=1:nShuffle
            [cueCells.shuffleScores(m,n), ~] = get_cuescoreYGModified_lessLagsNew(useDfofaverage(:,n),cueCells.shuffleTemplates(:,m),lagBins);
        end
    end
    
    %real scores
    cueCells.realScores=[];
    cueCells.realLags=[];
    for n=1:length(useIdx)
        [cueCells.realScores(1,n), cueCells.realLags(1,n)] = get_cuescoreYGModified_lessLagsNew(useDfofaverage(:,n),temp,lagBins);
    end
    
    cueCells.cueCellInUseIdx=[];
    for n=1:length(useIdx);
        a=cueCells.realScores(n);
        if ~isnan(cueCells.realScores(n))
            s=cueCells.shuffleScores(:,n);
            s=s(~isnan(s));
            if a>prctile(s,95);
                cueCells.cueCellInUseIdx(end+1,1)=n;
            end
        end
    end
    
    cueCells.cueCellRealIdx=cueCells.useIdx(cueCells.cueCellInUseIdx);
    cueCells.cueCellScores=cueCells.realScores(cueCells.cueCellInUseIdx);
    cueCells.cueCellLags=cueCells.realLags(cueCells.cueCellInUseIdx);
    cueCells.cueCellDfofAvg=cueCells.useDfofaverage(:,cueCells.cueCellInUseIdx);
    save('cueCells.mat','cueCells')
    
    cd ../

end



