function alignmentByIdx(alignIdx,fMaster,folders)
%% Process inputs

if nargin<=1 || isempty(fMaster)
    fMaster = pwd;
end

if nargin<=2 || isempty(idxs)
    load('folders.mat','folders');
end

if nargin==0 || isempty(alignIdx)
    alignIdx = length(folders);
end

cd(fMaster)
save('folders.mat','folders')


%%  Load files

nFold = length(folders);

name={};
allFile_names={};

for f=1:nFold
    disp(f)
    
    load(fullfile(folders{f},'allROIs.mat'));
    
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
    if f>alignIdx
        continue
    end

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
    
    if f==alignIdx
        itrIdx = 1:length(aa);
    else
        itrIdx = alignIdx-f;
    end

    for k = itrIdx
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
            run(fullfile(pwd,'commonCellIdentification.m'));
        end
        
        close all
    end
    
end
