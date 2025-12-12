function alignmentPCAICA(fMaster,folders)
%% Process inputs

if nargin==0 || isempty(fMaster)
    fMaster = pwd;
end

if nargin<2 || isempty(idxs)
    load('folders.mat','folders');
end

cd(fMaster)
save('folders.mat','folders')


%% Identify and process bad rois

load folders.mat
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

load folders.mat
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
        
        % run alignment and save data
        clear cell_registered_struct;
        
        d=dir([results_directory '\cellRegistered*']);
        if ~isempty(d)
            for ii = 1:length(d)
                delete([results_directory '\' d(ii).name])
            end
        end
        
        if isfile([results_directory '\manualAlignment.mat'])
            copyF = 'commonCellIdentificationManual.m';
            if ~isfile(copyF)
                srcF = which(copyF);
                copyfile(srcF);
            end
            run(fullfile(pwd,'commonCellIdentificationManual.m'))
        else
            try
                run(fullfile(pwd,'commonCellIdentification.m'));
            catch
            end
        end
        
        close all
    end
    
end


%% Fix bad cells

load ('badRois.mat','badFile','badRoi')
load folders.mat
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
