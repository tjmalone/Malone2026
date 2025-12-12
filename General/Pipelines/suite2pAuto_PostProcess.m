%% Use motion corrected data for Suit2p
% run suite2p directly using motion corrected tif stacks. Suite2p will
% output Fall.mat
%

clear; close all; clc

d = dir('*D3\*\loc*');
% d = dir('*D*\*\*\TSeries*_*\');


trackEnd = 200;

folders = cell(length(d),1);
for ii = 1:length(d)
    % folders{ii} = [d(ii).folder '\' d(ii).name '\suite2p'];
    folders{ii} = [d(ii).folder '\' d(ii).name];
end

p = pwd;

%%

for ii = 1:length(folders)
    if ~exist(folders{ii},'dir')
        continue;
    end
    
    cd(folders{ii})

    % check claim
    if exist('claimedS2.mat','file')~=0
        cd(p)
        disp(0)
        continue
    else
        claimed = 1;
        save('claimedS2.mat','claimed')
        disp(1)
    end
    
    % folder where code is saved
    baseFolder = 'D:\AnalysisCode\Pipelines\suite2p';
    
    % copy self to current directory
    copyfile([baseFolder '\suite2p_PostProcess_Interleave.m'],'suite2p_PostProcess_Interleave.m');
    % copyfile([baseFolder '\suite2p_PostProcess.m'],'suite2p_PostProcess.m');
    try
        suite2p_PostProcess_Interleave([],trackEnd,80,2,1)
        % suite2p_PostProcess([],trackEnd)
    catch
        delete('claimedS.mat')
    end
    
end

cd(p)
