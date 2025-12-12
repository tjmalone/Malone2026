%%

clear all; close all; clc

trackEnd = 400;     % length of your track
binWidth = 5;

base = 'D:\4m_optoPaper';
cd(base)

folds = findSubF('pcaica',3,[],0);

%%

% Run Functions:
% postFncs = {'dfofMCorrelationFull',...
%     'pValueRunByRun_dfof_oldSig_noRBR_par'...
%     'pValueRunByRun_sig_noRBR_par'...
%     'spatialInfo_dual'...
%     'speedScoreCalculation'...
%     'cueCellsAnalysis_dfof_DN_TM'...
%     'cueCellsAnalysis_sig_DN_TM'...
%     'popGeometryAnalysis};

runFncs = [4];
claimName = 'claimedS3.mat';

% cue templates must be copied separately beforehand
for ff = 1:length(folds)
    cd(folds{ff})
    
    % check claim
    if exist(claimName,'file')~=0
        cd(base)
        disp(0)
        continue
    else
        claimed = 1;
        save(claimName,'claimed')
        disp(1)
    end
    
    tic
    try
        copyfile('D:\AnalysisCode\PostAnalysis\postAnalysis_TM.m','postAnalysis_TM.m');
        postAnalysis_TM(trackEnd,binWidth,runFncs);
    catch err
        delete(claimName)
        delete('postAnalysis_TM.m')
        disp(err);
    end
    toc

end

cd(base)

