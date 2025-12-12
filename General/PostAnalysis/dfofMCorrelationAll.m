%%

clear all; close all; clc

trackEnd = 600;     % length of your track
binWidth = 5;
binWidthClean = 2.5;

base = 'D:\AD_Project\imagingData\data';
cd(base)

folds = findSubF('suite2p',4,[],0);


%%

for ff = 1:length(folds)
    cd(folds{ff})
    disp(ff)

    try
        % perform standard dfof calculation
        copyfile('D:\AnalysisCode\PostAnalysis\dfofMCorrelationFull.m','postAnalysis\dfofMCorrelationFull.m');
        dfofMCorrelationFull(trackEnd,binWidth);

        % perform dfof calculation with 2.5cm bin
        copyfile('D:\AnalysisCode\PostAnalysis\dfofMCorrelation_clean.m','postAnalysis\dfofMCorrelation_clean.m');
        dfofMCorrelation_clean(trackEnd,binWidthClean);

        % perform spatial RBR consistency calculation
        copyfile('D:\AnalysisCode\PostAnalysis\dfofM_RBRbyLoc.m','postAnalysis\dfofM_RBRbyLoc.m');
        dfofM_RBRbyLoc();
    catch
        disp('Error')
    end
end

cd(base)

