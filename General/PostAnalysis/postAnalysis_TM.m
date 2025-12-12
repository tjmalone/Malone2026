function postAnalysis_TM(trackLength,binWidth,runFncs)

p = pwd();
mkdir('postAnalysis');

errMessages = {};

close all;

fileSource = 'D:\AnalysisCode\PostAnalysis\';
postFncs = {'dfofMCorrelationFull',...
    'pValueRunByRun_dfof_oldSig_noRBR_par'...
    'pValueRunByRun_sig_noRBR_par'...
    'spatialInfo_dual'...
    'speedScoreCalculation_TM'...
    'cueCellsAnalysis_dfof_DN_TM'...
    'cueCellsAnalysis_sig_DN_TM'...
    'popGeometryAnalysis'};

if nargin<3 || isempty(runFncs)
    runFncs = 1:length(postFncs);
end

%%

errMessages = {};


for ii = 1:length(runFncs)
    curFnc = postFncs{runFncs(ii)};
    try
        copyfile([fileSource curFnc '.m'],[curFnc '.m'])
        feval(curFnc,trackLength,binWidth);
        cd(p)
        movefile([curFnc '.m'],['postAnalysis\' curFnc '.m']);
    catch er
        cd(p)
        disp(er)
        delete([curFnc '.m'])
        errMessages{end+1} = [curFnc ' '];
    end
        
    close all;
end

close all;

if ~isempty(errMessages)
   [~] = rmdir('postAnalysis'); % remove folder only if empty
   error(cat(2,errMessages{:}));
end

end

