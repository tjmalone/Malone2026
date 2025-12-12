function [useFOV,projection_sum,SFC] = alignmentPostCheck(FOV,cells,idxRef,cycleSave,projection_sum,SFC)
%% alignmentPostCheck
% Allows manual checking of input alignments. Will run based on FOV in
% current alignment folder.
%
% Inputs: Inputs with default can be empty or skipped
%   FOV - which FOV in the current alignment to analyze. Must be set as 1D
%       numerical array
%   cells - which cells to analyze. Must be 2D cell array (use num2cell to
%       convert numerical array). Each row should correspond to 1
%       alignment. Each column should correspond to 1 of the FOV. If a cell
%       contains multiple cell numbers, all will be plotted together.
%       However, order is not currently identified in plot
%   idxRef - reference FOV index. Must be set to an interger in FOV
%   cycleSave - whether to save user input of good/bad alignment. If set to
%       1, cycles through all alignment rows. Asks user input on good
%       alignment. User must enter numerical input, but any value is okay
%       (not just 0 and 1). If set to 0 (~1), all alignment rows will be
%       plotted simultaneously is separate figures. Default is 0.
%   projection_sum - pairwise alignment projections. If input, does not
%       need to repeat slow pairwise alignment step. Previous projection
%       sum must be run with same reference FOV.
%   SFC - spatial footprints corrected. If input, does not need to repeat
%       slow pairwise alignment step. Previous SFC must be run with same
%       reference FOV.
%
% Outputs:
%   useFOV - array of user input alignment quality. If cycleSvae is off,
%       array will be all zeros.
%   projection_sum - pairwise alignment projections. Can be used to run
%       additional checks without repeating alow pairwise alignment. Must
%       use same reference FOV.
%   SFC - spatial footprints correcteding alow pairwise alignment. Must
%       use same reference FOV.
%


%% Load alignment data

nFOV = length(FOV);     % number of FOV
nC = ceil(sqrt(nFOV));  % number of subplots

% find pairwise alignments to reference frame if alignments are not inputs
if nargin<5 || isempty(projection_sum) || isempty(SFC)
    % initialize alignment variables
    SFC = cell(2,nFOV-1);
    projection_sum = cell(1,nFOV-1);

    % cycle through FOV
    for ff = 1:nFOV
        % skip reference FOV
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
    end
end

%% Plot overlay

% set plotting inputs
colors = {'-b',':r'};   % border color/style
lWidths = [1.5,2];      % line widths
dBuff = 35;             % buffer around cell of interest to show

% initialize output array
useFOV = zeros(size(cells,1),1);

% cycle through all cells
for mm = 1:size(cells,1)
    figure;

    % cycle through all FOV
    for ff = 1:nFOV
        % skip reference FOV
        if FOV(ff)==idxRef; continue; end

        idxRefIdx = find(FOV==idxRef);
        % sort indices
        useN = [min(idxRefIdx,ff),max(idxRefIdx,ff)];

        % plot projection for FOV pair
        subplot(nC,nC,ff)
        imshow(projection_sum{ff})

        % set axis and title
        axis('off','equal')
        title(['Post-alignment: ' num2str(ff)],'FontWeight','Bold','fontsize',14)
        hold on

        % plot input cell overlays
        for ii = 1:2

            for jj = 1:length(cells{mm,useN(ii)})
                % skip empty cells
                if cells{mm,useN(ii)}(jj)==0
                    continue
                end

                % find current regions
                curRegion = squeeze(SFC{ii,ff}(cells{mm,useN(ii)}(jj),:,:));
                curRegion(curRegion>0) = 1;

                % calculate borders
                curBorder = bwboundaries(curRegion);

                % plot borders
                for kk = 1:length(curBorder)
                    plot(curBorder{kk}(:,2),curBorder{kk}(:,1),colors{ii},'linewidth',lWidths(ii))
                end

                % calculate bounding box
                curProps = regionprops(curRegion,'BoundingBox');

                % calculate extremes of bounding box
                bbAll = {curProps(:).BoundingBox};
                bbCat = cat(1,bbAll{:});
                bbMax = max(bbCat,[],1);

                % calculate region of interest
                Y = bbMax(2)-dBuff;
                X = bbMax(1)-dBuff;
                yD = bbMax(2)+bbMax(4)+dBuff;
                xD = bbMax(1)+bbMax(3)+dBuff;

                % set axis to region of interest
                axis([X,xD,Y,yD])
            end
        end
    end

    % set figure title/size
    sgtitle(['Cell ' num2str(mm)])
    set(gcf,'WindowState', 'maximized')

    % if cycle save is on, ask for user input
    if nargin>=4 && cycleSave==1
        try
            useFOV(mm) = input(['Cell ' num2str(mm) ' - good/bad (1/0): ']);
        catch
        end
        close all
    end

end

end

