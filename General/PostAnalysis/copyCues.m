%% copyCues

clear; clc

cueFolder = cell(3,1);
cueFolder{1} = 'D:\AnalysisCode\PostAnalysis\Cues\2m_KZ_passive';
cueFolder{2} = 'D:\AnalysisCode\PostAnalysis\Cues\2m_TM_passive_new';
cueFolder{3} = 'D:\AnalysisCode\PostAnalysis\Cues\2m_YM_env4_passsive';

% Env1 indices
env{1} = [74, 75, 80, 81];

% Env2 indices
env{2} = [77, 78];

% Env 3 indices
env{3} = [83, 84];
% env{4} = [47,48,50,51,53,54];

folds = findSubF('suite2p',4,[],0);
% folds = findSubF('suite2p',4,[],0);

%%
for ii = 1:length(folds)
%     if isfolder([folds{i} 'cueAnalysis_dfof'])
%         continue
%     end
    
    disp(folds{ii})
    
    for jj = 1:length(cueFolder)
        if ismember(ii,env{jj})
            copyfile(cueFolder{jj},[folds{ii} 'cueAnalysis_dfof'])
            copyfile(cueFolder{jj},[folds{ii} 'cueAnalysis_sig'])
            continue
        end
    end
end





