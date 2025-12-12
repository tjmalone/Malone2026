%% runTau_RC

% "processOverlay_TRC.m"
    % Need to edit slice directory in the code for each slice 
    % Must run for each slice before proceeding to step 2
        % Creates overlay images
        % Output:
            % Example: "overlay_calbindin-tau_5.fig"
            % path: Z:\SNM\labMembers\KC\AD_Project\histology_KC1.5mm\
            % reelCalTau\Batch6\TM240316-K\OtherFiles
        % Use these images to make the following zip files:
            % Overlay_reelin-tau: Reelin and tau overlap
            % Overlay_calbindin-tau: Calbindin and tau overlap
            % Full_overlay: Overlap of “Overlay_reelin-tau” and 
            % “Overlay_calbindin-tau”
                % Example folder: Z:\SNM\labMembers\KC\AD_Project
                % \histology_KC1.5mm\reelCalTau\Batch6\TM240316-K
        
run('calculateFinalCellCount_Manual.m')
% counts all the cells
% Output
    % cellCountFinal.mat
    % path: Z:\SNM\labMembers\KC\AD_Project\Histology\tau_RC\FinalVersion
        % has the numbers of the following cells for each FOV:
            % reelin
            % calbindin
            % tau
            % reelin-tau overlaps
            % calbindin-tau overlaps
            % reelin-calbindin-tau overlaps

run('plotOverlay_RCByFOV.m')
run('plotRCTBySexByFOV.m')
run('plotRCTBySex_forFigsByFOV.m')
% Makes plots

