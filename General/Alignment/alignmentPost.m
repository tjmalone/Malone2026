function alignmentPost(useIdx,svID,fMaster)
%% alignmentPost
% Generates an alignment results file based on a set of pre-run alignments.
% Generates the same output variables as the original alignment pipeline,
% but adds an alignment confidence calculation (the percent of FOV pairs
% that showed alignment for a given proposed alignment set.
%
% Input: Inputs can be left empty or skipped to use default behavior.
%   useIdx - FOV indices to process. Default is all FOV
%   svID - name of output data set. Save folder and results file both use
%       this as a suffix. Default is numerical numbering based on current
%       data saved to avoid overwriting
%   fMaster - folder to perform analysis in. Default is current folder
%
% Outputs: Save results_svID.mat, in results_svID folder. Contains general
%   alignment info
%


%% Process inputs

% set default svID directory
if nargin<2 || isempty(svID)
    d = dir('*\results*.mat');
    nD = length(d);
    svID = num2str(nD+1);
end

% set save file and folder name
svName = [svID '.mat'];
svFold = [svID '\'];
mkdir(svFold)

% set run folder
if nargin<3 || isempty(fMaster)
    fMaster = pwd;
end
orFolder = pwd;
cd(fMaster)

% set spatial footprint directories
% load('allFile_names.mat','allFile_names');
d2 = dir('spatial_footprints*.mat');
allFile_names = {d2(:).name};

% set FOV use indices
if nargin<1 || isempty(useIdx)
    useIdx = 1:length(allFile_names);
end
nIdx = length(useIdx);


%% Calculate number of cells identified under each reference

% initialize output structure
results = struct();
results.days={};
results.useIdx = useIdx;

% cycle each pair of FOV
for ff=1:nIdx
    f = useIdx(ff);
    disp(f)
    load(allFile_names{f});
    
    NRefCells=size(data,1);
    
    %save assignment data of ref with all other days
    results.days{ff}.allFOVs={};
    results.days{ff}.allRefOthers(:,1) = 1:NRefCells;
    kk=setdiff(useIdx,f);
    
    for k=1:length(kk)
        disp(kk(k))
        load(allFile_names{kk(k)});
        
        NCells=size(data,1);
        
        % move to correct FOV pair folder
        if kk(k)<f
            foldername=sprintf('%s_%s',num2str(kk(k)),num2str(f));
            cd(fullfile(num2str(kk(k)),foldername));
        else
            foldername=sprintf('%s_%s',num2str(f),num2str(kk(k)));
            cd(fullfile(num2str(f),foldername));
        end
        
        % load alignment pair
        d=dir('cellRegistered*');
        
        % save pairwise alignment info
        if ~isempty(d)
            load(d(1).name);
            cim=cell_registered_struct.cell_to_index_map;
            
            % sort cell indices according to current reference
            if kk(k)<f
                cim=fliplr(cim);
            end
            [~,ii]=sort(cim(:,1));
            cim=cim(ii,:);
            
            % move unaligned cells to end
            iNoZero=find(cim(:,1)~=0,1);
            cim = cat(1,cim(iNoZero:end,:),cim(1:iNoZero-1,:));
            
            results.days{ff}.allFOVs{k}=cim;
            results.days{ff}.allRefOthers(:,k+1)=cim(1:NRefCells,2);
            results.days{ff}.perctCommonCells(k)=sum(cim(1:NRefCells,2)~=0)/NCells;
            results.days{ff}.perctCommonCellsInRef(k)=sum(cim(1:NRefCells,2)~=0)/NRefCells;
        else
            error('No alignment found')
        end
        
        cd ..\..
    end
    
    results.days{ff}.perctCommonCellsMean=mean(results.days{ff}.perctCommonCells);
    results.days{ff}.perctCommonCellsInRefMean=mean(results.days{ff}.perctCommonCellsInRef);
    
    % number of common partners for individual cells in ref data
    % when the number is 1, means there is no common cell, only reference cell existed
    results.days{ff}.NCommonCells = sum(results.days{ff}.allRefOthers~=0,2);
    
    % sort these N Common cells in descending order
    [results.days{ff}.NCommonCellsSorted,i]=sort(results.days{ff}.NCommonCells,'descend');
    
    % sort allRefOther cells according to the order above.
    results.days{ff}.allRefOthersSorted=results.days{ff}.allRefOthers(i,:);
    
    % number of all common cells (including cells in ref days with common cells)
    results.days{ff}.totalNCommonCells=sum(results.days{ff}.NCommonCellsSorted(results.days{ff}.NCommonCellsSorted>1));
end

% save general common cell info
results.threshFrequentEnough=nIdx-2;
results.NFrequentEnough=[];
for f=1:nIdx
    results.totalNCommonCellsAll(f,1)=results.days{f}.totalNCommonCells;
    results.perctCommonCellsMeanAll(f,1)=results.days{f}.perctCommonCellsMean;
    results.perctCommonCellsInRefMeanAll(f,1)=results.days{f}.perctCommonCellsInRefMean;
    results.NFrequentEnough(f,1)= sum(results.days{f}.NCommonCellsSorted>=results.threshFrequentEnough);
end

save(svName,'results');


%% Plot common cell info

% plot number common cell number across days
figure
A=results.totalNCommonCellsAll;
B=results.perctCommonCellsMeanAll;
C=results.perctCommonCellsInRefMeanAll;
D=results.NFrequentEnough;
E = {A,B,C,D};
ylabs = {'Total N common cells','Perct common cells mean',...
    'Perct common cells mean','Number of frequent enough cells'};

subplot(3,2,1);
plot(A,B,'r.')
xlabel('Total N common cells');
ylabel('Perct common cells mean');

for ii = 1:4
    subplot(3,2,ii+1)
    plot(E{ii},'k','Marker','.','MarkerEdgeColor','r')
    xlabel('Days')
    ylabel(ylabs{ii})
end

saveas(gcf,[svFold 'results.fig']);
close

% plot percent common cell across days for all sessions
fldNames = {'perctCommonCells','perctCommonCellsInRef'};
ttlNames = {'perctCommonCellsFovs','perctCommonCellsRef'};
nSPlots = ceil(sqrt(nIdx));
for ii = 1:length(fldNames)
    figure
    for f=1:nIdx
        subplot(nSPlots,nSPlots,f);
        plot(results.days{f}.(fldNames{ii}),'k','Marker','.','MarkerEdgeColor','r');
        ylim([0 1]);
        title(['ref',num2str(f)]);
    end
    saveas(gcf,[svFold ttlNames{ii} '.fig']);
    close
end


%% Identify common cells

%first order all cells in the columns of FOVs(currently the reference is on
%the first column)
for n=1:nIdx
    results.days{n}.allRefOthersSortedInOrder = results.days{n}.allRefOthersSorted(:,[2:n 1 n+1:end]);
end

%extract common cells in all
for m=1:nIdx
    % identify all alignment groups for m common cells in all FOVs
    d=['commonCells_',num2str(m)];
    results.(d)=[];
    for n=1:nIdx
        % find all common cells sets using each FOV as a reference
        i=find(results.days{n}.NCommonCellsSorted==m);
        results.(d) = [results.(d); results.days{n}.allRefOthersSortedInOrder(i,:)];
    end
    
    % save only unique rows
    curCommon = unique(results.(d),'rows');
    results.(d) = curCommon;
    
    % calculate alignment confidence
    e = ['CCconfidence_',num2str(m)];
    ccConf = findConfidence(curCommon,results);
    results.(e) = ccConf;
end

% save output results
save([svFold svName],'results');

% return to original folder
cd(orFolder)

end

