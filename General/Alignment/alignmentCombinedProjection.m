function [projection_sum,SFC_use] = alignmentCombinedProjection(FOV,idxRef)
%% alignmentPostCheck
% Allows manual checking of input alignments. Will run based on FOV in
% current alignment folder.
%
% Inputs: Inputs with default can be empty or skipped
%   FOV - which FOV in the current alignment to analyze. Must be set as 1D
%       numerical array
%   idxRef - reference FOV index. Must be set to an interger in FOV
%
% Outputs:
%   projection_sum - pairwise alignment projections. Can be used to run
%       additional checks without repeating pairwise alignment. Must
%       use same reference FOV.
%


%% Generate aligned projections

nFOV = length(FOV);     % number of FOV

% initialize alignment variables
SFC = cell(2,nFOV-1);
SFC_use = cell(1,nFOV-1);
projection_sum = cell(1,nFOV-1);

% cycle through FOV
for ff = 1:nFOV
    if FOV(ff)==idxRef; continue; end

    % get pairwise alignment directory
    FOV1 = min(idxRef,FOV(ff));
    FOV2 = max(idxRef,FOV(ff));
    foldername = sprintf('%s_%s',num2str(FOV1),num2str(FOV2));
    results_directory = fullfile(pwd,num2str(FOV1),foldername);

    % load pairwise alignment
    d = dir([results_directory '\cellRegistered*mat']);
    load([d(end).folder '\' d(end).name],'cell_registered_struct')

    % load corrected alignment
    curSFC = cell_registered_struct.spatial_footprints_corrected;
    SFC(:,ff) = curSFC;

    % project correct alignment
    alignProjection = cell(1,2);
    alignProjection{1} = squeeze(sum(curSFC{1},1));
    alignProjection{2} = squeeze(sum(curSFC{2},1));
    for ii = 1:2
        % set all projection values to 1
        alignProjection{ii}(alignProjection{ii}>0) = 1;
    end

    % find alignment sum
    projection_sum{ff} = (alignProjection{1}+alignProjection{2})/2;

    % save individual  cell projections
    if idxRef<FOV(ff)
        idxCur = 2;
    else
        idxCur = 1;
    end
    SFC_use{:,ff} = logical(permute(curSFC{idxCur},[3 2 1]));
end


end

