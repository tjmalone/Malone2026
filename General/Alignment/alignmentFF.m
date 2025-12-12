function alignmentFF(folder1,folder2)
%% Define folders for analysis

% select the paths of all FOV to analyze
folders = {folder1,folder2};

save('folders.mat','folders');


%% Identify and process bad rois

load folders.mat folders
nFold = length(folders);

badFile = [];
badRoi = {};

for f = 1:nFold
    % indentify all bad rois
    bROIs = fixRoi(folders{f});
    
    % record files that have bad rois
    if ~isempty(bROIs)
        badFile(end+1) = f;
        badRoi{end+1} = bROIs;
    end
end

save('badRois.mat','badFile','badRoi')


%%  Load files

load folders.mat folders
nFold = length(folders);

name={};
allFile_names={};

for f=1:nFold
    disp(f)
    
    load(fullfile(folders{f},'allROIsNew.mat'));
    
    name{f}=sprintf('%s_%02d.mat','spatial_footprints',f);
    data=permute(roi,[3 2 1]);
    save(name{f},'data');
    allFile_names{f} = fullfile(pwd,name{f});
end

save('allFile_names.mat','allFile_names');

% copy 'commonCellIdentification.m' to local directory
copyF = 'commonCellIdentification.m';
if ~isfile(copyF)
    srcF = which(copyF);
    copyfile(srcF);
end


%% Align all files pairwise

load allFile_names.mat
nFold = length(allFile_names);

for f = 1:nFold-1
    clear results
    disp(f)
    
    % load reference day
    load(allFile_names{f});
    NRefCells=size(data,1);         % reference roi number
    mkdir(num2str(f));              % directory for reference day
    aa = f+1:1:nFold;     % indices for non-ref days
    
    % initialize results struct
    results.refDay = f; % reference day number
    results.allFOVs = {}; % combined assignment data with other days
    
    % [reference index, common cell index for other days]
    results.allRefOthers(:,1) = 1:1:NRefCells;
    results.perctCommonCells = []; % percentage of overlapping cells
    
    for k = 1:length(aa)
        % load current non-ref data
        load(allFile_names{aa(k)});
        NCells = size(data,1);
        
        % make directory for current alignment (ref day - comparison day)
        foldername = sprintf('%s_%s',num2str(f),num2str(aa(k)));
        mkdir(fullfile(num2str(f),foldername));
        useIdx = [f aa(k)]; % indices of comparison days
        
        % set inputs for alignment code
        useFile_names = allFile_names(useIdx);
        results_directory = fullfile(pwd,num2str(f),foldername);
        save([results_directory '/useFile_names.mat'],'useFile_names');
        
        % delete old alignment data
        name='\cellRegistered*.mat';
        d=dir([results_directory name]);
        for ii = 1:length(d)
            delete([results_directory '\' d(ii).name]);
        end
        
        % run alignment and save data
        clear cell_registered_struct;
        
        if isfile([results_directory '\manualAlignment.mat'])
            copyF = 'commonCellIdentificationManual.m';
            if ~isfile(copyF)
                srcF = which(copyF);
                copyfile(srcF);
            end
            run(fullfile(pwd,'commonCellIdentificationManual.m'))
        else
            run(fullfile(pwd,'commonCellIdentification.m'));
        end
        
        close all
    end
    
end


%% Fix bad cells

load ('badRois.mat','badFile','badRoi')
load folders.mat folders
p = pwd;

for i = 1:length(badFile)
    %% Define fix folders
    
    pre = 1:badFile(i)-1;                   % bad file is not reference
    post = badFile(i)+1:length(folders);    % bad file is reference
    
    fixFolders = {cell(size(post)),cell(size(pre))};
    
    str = '\\%d\\%d_%d';
    
    % full path of post folders
    for j = 1:length(post)
        fixFolders{1}{j} = [p sprintf(str,badFile(i),badFile(i),post(j))];
    end
    
    % full path of pre folders
    for j = 1:length(pre)
        fixFolders{2}{j} = [p sprintf(str,pre(j),pre(j),badFile(i))];
    end
    
    
    %% Fix folders
    
    % fix all post and pre folder alignments
    for n = 1:2
        for f=1:length(fixFolders{n})
            % change to alignment directory for current pair
            cd(fixFolders{n}{f});
            
            % identify current alignment file
            name='cellRegistered*.mat';
            d=dir(name);
            load(d(1).name);
            
            A = cell_registered_struct.cell_to_index_map;
            
            % cycle through each bad roi
            for b = 1:length(badRoi{i})
                
                % increase index for all rois after bad roi
                for k=1:size(A,1)
                    if A(k,n)>=badRoi{i}(b)
                        A(k,n)=A(k,n)+1;
                    end
                end
                
                % add back bad roi base on pre/post status
                if n==1
                    A(end+1,:) = [badRoi{i}(b) 0];
                elseif n==2
                    A(end+1,:) = [0 badRoi{i}(b)];
                end
            end
            
            cell_registered_struct.cell_to_index_map = A;
            save('cellRegisteredFix.mat','cell_registered_struct');
            
        end
    end
end

cd(p)


%% Calculate number of cells identified under each reference

load('allFile_names.mat');
load('badRois.mat');
nFold = length(allFile_names);

clear results;
results.days={};

for f=1:nFold
    disp(f)
    load(allFile_names{f});
    
    badF = find(badFile==f);
    if ~isempty(badF)
        NRefCells=size(data,1)+length(badRoi{badF});
    else
        NRefCells=size(data,1);
    end
    
    %save assignment data of ref with all other days
    results.days{f}.allFOVs={};
    results.days{f}.allRefOthers(:,1) = 1:NRefCells;
    aa=setdiff(1:nFold,f);
    
    for k=1:length(aa)
        disp(k)
        load(allFile_names{aa(k)});
        
        badF = find(badFile==aa(k));
        if ~isempty(badF)
            NCells=size(data,1)+length(badRoi{badF});
        else
            NCells=size(data,1);
        end
        
        if aa(k)<f
            foldername=sprintf('%s_%s',num2str(aa(k)),num2str(f));
            cd(fullfile(num2str(aa(k)),foldername));
            d=dir('cellRegistered*');
            
            if ~isempty(d)
                load(d(1).name);
                ff=cell_registered_struct.cell_to_index_map;
                ff=fliplr(ff);
                [~,ii]=sort(ff(:,1));
                ff=ff(ii,:);
                iNoZero=find(ff(:,1));
                iZero=setdiff([1:1:size(ff,1)],iNoZero);
                fNoZero=ff(iNoZero,:);
                fZero=ff(iZero,:);
                ff=cat(1,fNoZero,fZero);
                results.days{f}.allFOVs{k}=ff;
                results.days{f}.allRefOthers(:,k+1)=ff(1:NRefCells,2);
                results.days{f}.perctCommonCells(k)=length(find(results.days{f}.allRefOthers(:,k+1)))/NCells;
                results.days{f}.perctCommonCellsInRef(k)=length(find(results.days{f}.allRefOthers(:,k+1)))/NRefCells;
            else
                
                e1=zeros(NRefCells,1);
                e2=zeros(NCells,1);
                results.days{f}.allFOVs{k}(:,1)= cat(1,[1:1:NRefCells]',e2);
                results.days{f}.allFOVs{k}(:,2)= cat(1,e1,[1:1:NCells]');
                results.days{f}.allRefOthers(:,k+1)=zeros(NRefCells,1);
                results.days{f}.perctCommonCells(k)=0;
                results.days{f}.perctCommonCellsInRef(k)=0;
            end
            
        else
            foldername=sprintf('%s_%s',num2str(f),num2str(aa(k)));
            cd(fullfile(num2str(f),foldername));
            d=dir('cellRegistered*');
            if ~isempty(d)
                load(d(1).name);
                ff=cell_registered_struct.cell_to_index_map;
                %              f=fliplr(f);
                [~,ii]=sort(ff(:,1));
                ff=ff(ii,:);
                iNoZero=find(ff(:,1));
                iZero=setdiff(1:size(ff,1),iNoZero);
                fNoZero=ff(iNoZero,:);
                fZero=ff(iZero,:);
                ff=cat(1,fNoZero,fZero);
                results.days{f}.allFOVs{k}=ff;
                results.days{f}.allRefOthers(:,k+1)=ff(1:NRefCells,2);
                results.days{f}.perctCommonCells(k)=length(find(results.days{f}.allRefOthers(:,k+1)))/NCells;
                results.days{f}.perctCommonCellsInRef(k)=length(find(results.days{f}.allRefOthers(:,k+1)))/NRefCells;
            else
                
                e1=zeros(NRefCells,1);
                e2=zeros(NCells,1);
                results.days{f}.allFOVs{k}(:,1)= cat(1,[1:1:NRefCells]',e2);
                results.days{f}.allFOVs{k}(:,2)= cat(1,e1,[1:1:NCells]');
                results.days{f}.allRefOthers(:,k+1)=zeros(NRefCells,1);
                results.days{f}.perctCommonCells(k)=0;
                results.days{f}.perctCommonCellsInRef(k)=0;
            end
            
        end
        cd ..\..
    end
    
    results.days{f}.perctCommonCellsMean=mean(results.days{f}.perctCommonCells);
    results.days{f}.perctCommonCellsInRefMean=mean(results.days{f}.perctCommonCellsInRef);
    
    results.days{f}.NCommonCells=[];%number of common partners for individual cells in ref data
    for mm=1:NRefCells
        results.days{f}.NCommonCells(mm,1)=length(find(results.days{f}.allRefOthers(mm,:)));%when the number is 1, means there is no common cell, only reference cell existed
    end
    [results.days{f}.NCommonCellsSorted,i]=sort(results.days{f}.NCommonCells,'descend'); %sort these N Common cells in descending order
    results.days{f}.allRefOthersSorted=results.days{f}.allRefOthers(i,:); %sort allRefOther cells according to the order above.
    results.days{f}.totalNCommonCells=sum(results.days{f}.NCommonCellsSorted(find(results.days{f}.NCommonCellsSorted-1))); %number of all common cells (including cells in ref days); Ref cells with no common cells are not counted. Ref cells with common cells are included
end

results.threshFrequentEnough=nFold-2;
results.NFrequentEnough=[];
for f=1:nFold
    results.totalNCommonCellsAll(f,1)=results.days{f}.totalNCommonCells;
    results.perctCommonCellsMeanAll(f,1)=results.days{f}.perctCommonCellsMean;
    results.perctCommonCellsInRefMeanAll(f,1)=results.days{f}.perctCommonCellsInRefMean;
    results.NFrequentEnough(f,1)=length(find(results.days{f}.NCommonCellsSorted>=results.threshFrequentEnough));
end

save('results.mat','results');

figure
A=results.totalNCommonCellsAll;
B=results.perctCommonCellsMeanAll;
C=results.perctCommonCellsInRefMeanAll;
D=results.NFrequentEnough;
subplot(321);
plot(A,B,'r.')
xlabel('Total N common cells');
ylabel('Perct common cells mean');
subplot(322)
plot(A,'k');
hold on
plot(A,'r.');
ylabel('Total N common cells');
xlabel('Days');
subplot(323)
plot(B,'k');
hold on
plot(B,'r.');
ylabel('Perct common cells mean');
xlabel('Days');
subplot(324)
plot(C,'k');
hold on
plot(C,'r.');
ylabel('Perct common cells in ref mean');
xlabel('Days');
saveas(gcf,'results.fig');
subplot(325)
plot(D,'k');
hold on
plot(D,'r.');
ylabel('Number of frequent enough cells');
xlabel('Days');
saveas(gcf,'results.fig');
close


figure
for f=1:length(results.days)
    subplot(4,4,f);
    plot(results.days{f}.perctCommonCells,'k');
    hold on
    plot(results.days{f}.perctCommonCells,'r.');
    ylim([0 1]);
    title(['ref',num2str(f)]);
end

saveas(gcf,'perctCommonCellsFovs.fig');
close

figure
for f=1:length(results.days)
    subplot(4,4,f);
    plot(results.days{f}.perctCommonCellsInRef,'k');
    hold on
    plot(results.days{f}.perctCommonCellsInRef,'r.');
    ylim([0 1]);
    title(['ref',num2str(f)]);
end

saveas(gcf,'perctCommonCellsRef.fig');
close


%% identify common cells
%first order all cells in the columns of FOVs(currently the reference is on
%the first column)

load('results.mat');

for n=1:length(results.days)
    results.days{n}.allRefOthersSortedInOrder=[];
    
    if n==1
        results.days{n}.allRefOthersSortedInOrder=results.days{n}.allRefOthersSorted;
    else
        results.days{n}.allRefOthersSortedInOrder(:,n)=results.days{n}.allRefOthersSorted(:,1);
        idxBefore=[2:1:n];
        idxAfter=[n+1:1:size(results.days{n}.allRefOthersSorted,2)];
        results.days{n}.allRefOthersSortedInOrder(:,1:n-1)=results.days{n}.allRefOthersSorted(:,idxBefore);
        results.days{n}.allRefOthersSortedInOrder(:,n+1:size(results.days{n}.allRefOthersSorted,2))=results.days{n}.allRefOthersSorted(:,idxAfter);
    end
end

%extract common cells in all
for n=1:length(results.days)
    d=['commonCells_',num2str(n)];
    results.(d)=[]; %common number of n cells in all FOVs
    for m=1:length(results.days)
        i=find(results.days{m}.NCommonCellsSorted==n);
        if ~isempty(i)
            results.(d)(end+1:end+length(i),:)=results.days{m}.allRefOthersSortedInOrder(i,:);
        end
    end
    [results.(d)]=unique(results.(d),'rows');
end


save('results.mat','results');
