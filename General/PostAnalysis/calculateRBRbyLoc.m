%%

clear all; close all; clc

base = 'D:\AD_Project\imagingData';
cd(base)

folds = findSubF('suite2p',4,[],0);

%%

for ff = 1:length(folds)
    cd(folds{ff})
    disp(ff)
    try
        copyfile('D:\AnalysisCode\PostAnalysis\dfofM_RBRbyLoc.m','postAnalysis\dfofM_RBRbyLoc.m');
        dfofM_RBRbyLoc();
    catch
        disp('Error')
    end
end

cd(base)

