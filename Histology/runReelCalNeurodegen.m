%% runReelCalNeurodegen

% 1. Save indices of rois that are 90% percent within the MEC layer 2 roi
run('remBorderCells.m')
    % Output: removeBorderCells.mat
    % path: Z:\SNM\labMembers\KC\AD_Project\Histology\reelCalNeurodegen\FinalVersion\remBorderCells
        % indices of rois that pass the border cutoff

% 1b.Creates new roi files for reelin and calbindin. 
% Allows you to see the final rois that go into the calculations. 
% Not a necessary step.
% saveFilteredRoisReelCal.m	

% 2. count cells per area
run('countCellsPerArea_remBord.m')
    % Output: countCellsPerArea.mat
    % path: Z:\SNM\labMembers\KC\AD_Project\Histology\reelCalNeurodegen\FinalVersion\countCellsPerArea
        % roisPerArea: rois per area for all FOVS

% 3. plot
run('plotCellsPerAreaRC_GcampMice.m')


