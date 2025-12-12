function [SFC,CLC,FPC,CPC,MCC,AT,OF] = FOV_manual_alignment_new( fovCells, xdiff, ydiff, degRotation, goAhead, resizeFactor_W, resizeFactor_H)

% INPUT
%   fov1Dir: directory contains 'allROIs.mat' for fov 1
%   fov2Dir: directory contains 'allROIs.mat' for fov 2
%   fovCells: cells for use in alignment (from allROIs.mat)
%
%   xdiff = displacement in x-axis
%       calculate by substracting x-position of fov2 in Illustrator from
%       that of fov1
%   ydiff = displacement in y-axis
%       calculate by substracting y-position of fov2 in Illustrator from
%       that of fov1
%
%   degRotation = Rotate degree of fov2 as indicated in Illustrator
%
%   cellOverlapRatio = threshold of area overlapped percentage to be
%       identified as same cell
%
%   goAhead = 1 or 0, 1 to go ahead and align cell ID
%
% OUTPUT
%   plot overlapped fov
%   cellAligned = pair of overlapped cells
%
%   SFC = spatial_footprints_corrected
%   FPC = footprints_projections_corrected
%   CLC = centroid_locations_corrected
%   CPC = centroid_projections_corrected
%   MCC = maximal_cross_correlation
%   AT = alignment_translations
%   OF = overlapping_FOV
%

tic

warning('off')

roi1 = fovCells{1};
roi2 = fovCells{2};

deg = degRotation;
x0 = size(roi2, 1)/2;
y0 = size(roi2, 2)/2;

sz1 = size(roi1,1:2);
sz2 = size(roi2,1:2);
mxSz = max(sz1,sz2);
roi2_aligned = zeros(size(roi2));

% ROTATE and RESIZE
for cc = 1:size(roi2, 3)+1
    if cc<=size(roi2, 3)
        xx = roi2(:, :, cc);
    else
        xx = ones(sz2);
    end
    
    [row0, col0] = find(xx==1);
    
    % RESIZE
    row0 = row0 * resizeFactor_H;   % resize heigth
    col0 = col0 * resizeFactor_W;   % resize width
    
    % ROTATE
    row = cosd(deg)*(row0-x0) - sind(deg)*(col0-y0) + x0 + ydiff;
    row = round(row);
    
    col = sind(deg)*(row0-x0) + cosd(deg)*(col0-y0) + y0 + xdiff;
    col = round(col);
    
    useIdx = row>0 & row<=sz2(1) & col>0 & col<sz2(2);
    indices = sub2ind(sz2, row(useIdx),col(useIdx));
    tempCellROI = zeros(sz2);
    tempCellROI(indices) = 1;
    tempCellROI = imfill(tempCellROI,'holes');
    
    if cc<=size(roi2, 3)
        roi2_aligned(:, :, cc) = tempCellROI;
    else
        % overlapping_FOV
        OF = padarray(tempCellROI,mxSz-sz2,'post')';
    end
end

oldROI = (sum(roi1, 3)>=1)*2;
oldROI = padarray(oldROI,mxSz-sz1,'post');

newROI = (sum(roi2_aligned,3))>=1;
newROI = padarray(newROI,mxSz-sz2,'post');

overlappedRoi = oldROI + newROI;

figure; imagesc(overlappedRoi);
axis('off','equal')
% magenta-green-white
% map = [0 0 0;
%     0 1 0;
%     1 0 1;
%     1 1 1];

% red-green-yellow
map = [0 0 0;
    0 1 0;
    1 0 0;
    1 1 0];
colormap(map);


%% Generate Outputs

if goAhead
    microns_per_pixel=0.8;
    
    SFC_raw = cell(1,2);        % spatial_footprints_corrected
    FPC = cell(1,2);        % footprints_projections_corrected
    
    ROIs = cell(1,2);
    ROIs{1} = roi1;
    ROIs{2} = roi2_aligned;
    
    % spatial_footprints_corrected
    for ii = 1:2
        SFC_raw{ii} = permute((ROIs{ii}),[3 2 1]);  
    end
%     SFC_norm = normalize_spatial_footprints(SFC_raw);
    SFC = adjust_FOV_size(SFC_raw);
    
    % footprints_projections_corrected
    for ii = 1:2
        FPC{ii} = squeeze(sum(SFC{ii},1));      
    end
    % centroid_locations_corrected
    CLC = compute_centroid_locations(SFC,microns_per_pixel);
    
    % centroid_projections_corrected
    CPC = compute_centroids_projections(CLC,SFC);
    
    % maximal_cross_correlation
    full_FOV_correlation = normxcorr2(FPC{1},FPC{2});
    MCC = max(max(full_FOV_correlation));
    
    % alignment_translations
    AT = [xdiff,resizeFactor_W;ydiff,resizeFactor_H;0,degRotation];
    
    % figure;
    % % overlapping_FOV
    % imagesc(OF)
end

toc

end

