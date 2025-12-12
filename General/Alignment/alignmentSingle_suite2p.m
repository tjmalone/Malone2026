function alignmentSingle_suite2p(~,idx1,idx2)
%% Align all files pairwise

load allFile_names.mat

for f = idx1
    clear results
    disp(f)
    
    % load reference day
    load(allFile_names{f});
    NRefCells=size(data,1);         % reference roi number
    mkdir(num2str(f));              % directory for reference day
    aa = idx2;     % indices for non-ref days
    
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
            copyF = 'commonCellIdentificationNR.m';
            if ~isfile(copyF)
                srcF = which(copyF);
                copyfile(srcF);
            end

            run(fullfile(pwd,'commonCellIdentificationNR.m'));
        end
        
        close all
    end
    
end
