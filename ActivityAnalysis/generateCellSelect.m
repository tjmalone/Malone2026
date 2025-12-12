%% generateCellSelect
% Generates a structure (cellSelect) containing the cell indices for all
% analysis subsets. The cellSelect structure can be used with an "by cell"
% metric or distribution. By FOV or by cell avergaing should be performed
% after integrating a cellSelect leaf with a collected data metric.
%
% cellSelect:
%   Struct - learning type (learn/recall)
%   Struct - sex (allSex/female/male)
%   Struct - cell type (allType/common/grid/cue/non-common)
%   Struct - morphology type (allMorph/ste/pyr)
%   Struct - genotype (WT/AD)
%   Cell Array - session (FOV x Day)
%   Num. Array - use cells (empty for excleded mice)
%

clear; close all; clc

p1 = 'D:\AD_Project\imagingData';
cd(p1)

% load learning alignments
load('foldersLearning.mat')

% load recall alignments
load('foldersRecall.mat')

% load genotype and sex information
load('groupIDs.mat','groups','sexes')
load('groupIDsRecall.mat','groupsRecall','sexesRecall')

% load learning cell type information
load('data\globalCellTypes.mat','globalCellTypeIdx')
learnTypes = globalCellTypeIdx;

% load learning cell type information
load('data\globalCellTypesRecall.mat','globalCellTypeIdx')
recallTypes = globalCellTypeIdx;

% load morphology information
load('D:\AD_Project\imagingData\data\morph\morphIdxs.mat','morphIdxs')
useMorphs = morphIdxs.commonCells;

% load recall morphology information
load('D:\AD_Project\imagingData\data\morph\morphIdxsRecall.mat','morphIdxs')
useMorphsRecall = morphIdxs;


%% Initizalize structure

% loop key:
% learning type (lt), sex (ss), cell type (ct), morphology type (mt),
% genotype (gg), FOV (ff), day (dd)

learningTypes = {'learn','recall'};
sexIDs = {'allSex','female','male'};
cellTypes = {'allType','common','grid','cue','other','nonCommon','nongrid'};
morphTypes = {'allMorph','ste','pyr'};
genotypes = {'WT','AD'};

% set structure sizes
nLearningTypes = length(learningTypes);
nSexIDs = length(sexIDs);
nCellTypes = length(cellTypes);
nMorphTypes = length(morphTypes);
nGenotypes = length(genotypes);

% initialize cellSelect struct
cellSelect = struct();

% loop through learning types
for lt = learningTypes
    % set current folder and alignment info
    if strcmp(lt,'learn')
        curFolders = foldersLearning;
    elseif strcmp(lt,'recall')
        curFolders = foldersRecall;
    else
        error('Fix learing type list')
    end

    % set FOV and day number
    nFOV = length(curFolders);
    nDays = length(curFolders{1});

    % loop through sexes
    for ss = sexIDs
        % loop through cell types
        for ct = cellTypes
            % loop through morphologies
            for mt = morphTypes
                % skip cell type morphology interactions
                if strcmp(ct,'allType') && ~strcmp(mt,'allMorph') && strcmp(lt,'learn')
                    continue
                end

                % loop through genotypes
                for gg = genotypes
                    cellSelect.(lt{:}).(ss{:}).(ct{:}).(mt{:}).(gg{:}) =...
                        cell(nFOV,nDays);
                end
            end
        end
    end
end


%% Store sub-analysis indices

% loop through learning types
for lt = learningTypes
    % set current folder and alignment info
    if strcmp(lt,'learn')
        curFolders = foldersLearning;
        curAligns = alignsLearning;
        curGlobalTypes = learnTypes;
        curGroups = groups;
        curSexes = sexes;
        curMorphs = useMorphs;
    elseif strcmp(lt,'recall')
        curFolders = foldersRecall;
        curAligns = alignsRecall;
        curGlobalTypes = recallTypes;
        curGroups = groupsRecall;
        curSexes = sexesRecall;
        curMorphs = useMorphsRecall;
    else
        error('Fix learning type list')
    end

    % set FOV and day number
    nFOV = length(curFolders);
    nDays = length(curFolders{1});

    % loop through FOV
    for ff = 1:nFOV
        % loop through days
        for dd = 1:nDays

            % correct for nan days
            if strcmp(lt,'learn')
                % identify true day index
                trueDay = trueDays(ff,dd);

                % skip days past day limits
                if trueDay>nDays
                    continue
                end
            else
                trueDay = dd;
            end

            % get full cell indices
            roi = matfile(fullfile(curFolders{ff}{dd},'allROIs.mat'));
            nCells = size(roi,'roi',3);

            % loop through cell types
            for ct = 1:length(cellTypes)
                clear idxType
                if strcmp(cellTypes{ct},'allType')
                    idxType = 1:nCells;
                elseif strcmp(cellTypes{ct},'common')
                    idxType = curAligns{ff}(:,dd);
                elseif strcmp(cellTypes{ct},'grid')
                    idxType = curGlobalTypes.grid_sig{ff}(:,dd);
                elseif strcmp(cellTypes{ct},'cue')
                    idxType = [curGlobalTypes.cueL_sig{ff}(:,dd);...
                        curGlobalTypes.cueR_sig{ff}(:,dd);...
                        curGlobalTypes.cueAll_sig{ff}(:,dd)];
                elseif strcmp(cellTypes{ct},'other')
                    idxBase = curAligns{ff}(:,dd);
                    idxRemove = [curGlobalTypes.grid_sig{ff}(:,dd);...
                        curGlobalTypes.cueL_sig{ff}(:,dd);...
                        curGlobalTypes.cueR_sig{ff}(:,dd);...
                        curGlobalTypes.cueAll_sig{ff}(:,dd)];
                    idxType = setdiff(idxBase,idxRemove);
                elseif strcmp(cellTypes{ct},'nonCommon')
                    idxType = setdiff(1:nCells,curAligns{ff}(:,dd));
                elseif strcmp(cellTypes{ct},'nongrid')
                    idxType = setdiff(curAligns{ff}(:,dd),curGlobalTypes.grid_sig{ff}(:,dd));
                end

                % loop through morphologies
                for mt = 1:length(morphTypes)
                    % skip cell type morphology interactions
                    if strcmp(cellTypes{ct},'allType') &&...
                            ~strcmp(morphTypes{mt},'allMorph') && strcmp(lt,'learn')
                        continue
                    end

                    if strcmp(morphTypes{mt},'allMorph')
                        idxMorph = 1:nCells;
                    else
                        idxMorph = curMorphs.(morphTypes{mt}){ff,dd};
                    end

                    % determine final cell indices
                    idxFinal = intersect(idxType,idxMorph,'stable');

                    % loop through sexes
                    for ss = 1:length(sexIDs)
                        % skip unmatches sex
                        if ~ismember(ff,curSexes{ss})
                            continue
                        end

                        % loop through genotypes
                        for gg = 1:length(genotypes)
                            % skip unmatched genotype
                            if ~ismember(ff,curGroups{gg})
                                continue
                            end

                            % store final indices
                            cellSelect.(lt{:}).(sexIDs{ss}).(cellTypes{ct}).(morphTypes{mt}).(genotypes{gg}){ff,trueDay}...
                                = idxFinal;

                            % disp([lt{:} '-' sexIDs{ss} '-' cellTypes{ct} '-'...
                            %     morphTypes{mt} '-' genotypes{gg} '-'...
                            %     num2str(ff) '-' num2str(trueDay)])
                        end
                    end
                end
            end
        end
    end
end

cd(p1)
save('data\cellSelect.mat','cellSelect')

