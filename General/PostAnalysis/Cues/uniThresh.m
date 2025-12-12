
folders=[];
folders{1}='E:\ID20201118';
folders{2}='E:\ID20210206';
folders{3}='E:\ID20210207';
folders{4}='E:\ID20210208';
folders{5}='E:\ID20210209';
folders{6}='E:\ID20210413';
folders{7}='E:\ID20210519_1';
folders{8}='E:\ID20210519_2';
save('folders.mat','folders');


%% old env
clear all
p=pwd;
load('folders.mat');

allShufflesLDfofOld=[];
allShufflesRDfofOld=[];
allShufflesLSigOld=[];
allShufflesRSigOld=[];

for n=1:length(folders);
    cd(folders{n});
    load('oldEnvPath.mat');
    for idx=1:length(oldEnvPath);
        disp(idx)
        cd(oldEnvPath{idx});
        p2=pwd;
        d=dir;
        %remove non-folder files (isdir property is 0)
        dfolders = d([d(:).isdir]) ;
        % remove '.' and '..'
        dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));
        
        for nn=1:length(dfolders);
            cd([p2 '\' dfolders(nn).name]);
            if exist('pcaica','dir');
                cd('pcaica');
                if exist('allROIs.mat')>0;
                    load('allROIs.mat');
                    N=size(roi,3);
                    if N>1;
                        load('cueAnalysisNew_dfof\newScoreShuffleTemplate\Left\cueCells.mat');
                        A=cueCells.shuffleScores;
                        A=reshape(A,size(A,1)*size(A,2),1);
                        allShufflesLDfofOld(end+1:end+length(A))=A;
                        load('cueAnalysisNew_dfof\newScoreShuffleTemplate\Right\cueCells.mat');
                        A=cueCells.shuffleScores;
                        A=reshape(A,size(A,1)*size(A,2),1);
                        allShufflesRDfofOld(end+1:end+length(A))=A;
                        
                        load('cueAnalysisNew_sig\newScoreShuffleTemplate\Left\cueCells.mat');
                        A=cueCells.shuffleScores;
                        A=reshape(A,size(A,1)*size(A,2),1);
                        allShufflesLSigOld(end+1:end+length(A))=A;
                        load('cueAnalysisNew_sig\newScoreShuffleTemplate\Right\cueCells.mat');
                        A=cueCells.shuffleScores;
                        A=reshape(A,size(A,1)*size(A,2),1);
                        allShufflesRSigOld(end+1:end+length(A))=A;
                    end
                end
            end
        end
    end
end
cd(p);

% new env
p=pwd;
load('folders.mat');

allShufflesLDfofNew=[];
allShufflesRDfofNew=[];
allShufflesLSigNew=[];
allShufflesRSigNew=[];

for n=1:length(folders);
    cd(folders{n});
    load('newEnvPath.mat');
    for idx=1:length(newEnvPath);
        disp(idx)
        cd(newEnvPath{idx});
        p2=pwd;
        d=dir;
        %remove non-folder files (isdir property is 0)
        dfolders = d([d(:).isdir]) ;
        % remove '.' and '..'
        dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));
        
        for nn=1:length(dfolders);
            cd([p2 '\' dfolders(nn).name]);
            if exist('pcaica','dir');
                cd('pcaica');
                if exist('allROIs.mat')>0;
                    load('allROIs.mat');
                    N=size(roi,3);
                    if N>1;
                        load('cueAnalysisNew_dfof\newScoreShuffleTemplate\Left\cueCells.mat');
                        A=cueCells.shuffleScores;
                        A=reshape(A,size(A,1)*size(A,2),1);
                        allShufflesLDfofNew(end+1:end+length(A))=A;
                        load('cueAnalysisNew_dfof\newScoreShuffleTemplate\Right\cueCells.mat');
                        A=cueCells.shuffleScores;
                        A=reshape(A,size(A,1)*size(A,2),1);
                        allShufflesRDfofNew(end+1:end+length(A))=A;
                        
                        load('cueAnalysisNew_sig\newScoreShuffleTemplate\Left\cueCells.mat');
                        A=cueCells.shuffleScores;
                        A=reshape(A,size(A,1)*size(A,2),1);
                        allShufflesLSigNew(end+1:end+length(A))=A;
                        load('cueAnalysisNew_sig\newScoreShuffleTemplate\Right\cueCells.mat');
                        A=cueCells.shuffleScores;
                        A=reshape(A,size(A,1)*size(A,2),1);
                        allShufflesRSigNew(end+1:end+length(A))=A;
                    end
                end
            end
        end
    end
end
cd(p);

thsLDfofNew=prctile(allShufflesLDfofNew,95);
thsRDfofNew=prctile(allShufflesRDfofNew,95);
thsLSigNew=prctile(allShufflesLSigNew,95);
thsRSigNew=prctile(allShufflesRSigNew,95);
%

% thsLDfofNew=prctile([allShufflesLDfofNew allShufflesLDfofOld],95);
% thsRDfofNew=prctile([allShufflesRDfofNew allShufflesRDfofOld],95);
% thsLSigNew=prctile([allShufflesLSigNew allShufflesLSigOld],95);
% thsRSigNew=prctile([allShufflesRSigNew allShufflesRSigOld],95);

save('thsLDfofNew.mat','thsLDfofNew');
save('thsRDfofNew.mat','thsRDfofNew');
save('thsLSigNew.mat','thsLSigNew');
save('thsRSigNew.mat','thsRSigNew');

% thsLDfofOld=thsLDfofNew;
% thsRDfofOld=thsRDfofNew;
% thsLSigOld=thsLSigNew;
% thsRSigOld=thsRSigNew;

thsLDfofOld=prctile(allShufflesLDfofOld,95);
thsRDfofOld=prctile(allShufflesRDfofOld,95);
thsLSigOld=prctile(allShufflesLSigOld,95);
thsRSigOld=prctile(allShufflesRSigOld,95);

save('thsLDfofOld.mat','thsLDfofOld');
save('thsRDfofOld.mat','thsRDfofOld');
save('thsLSigOld.mat','thsLSigOld');
save('thsRSigOld.mat','thsRSigOld');

%%
p=pwd;

per=[];
load('folders.mat');
for n=1:length(folders);
    cd(folders{n});
    load('oldEnvPath.mat');
    
    for idx=1:length(oldEnvPath);
        disp(idx)
        cd(oldEnvPath{idx});
        p2=pwd;
        d=dir;
        %remove non-folder files (isdir property is 0)
        dfolders = d([d(:).isdir]) ;
        % remove '.' and '..'
        dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));
        
        for nn=1:length(dfolders);
            cd([p2 '\' dfolders(nn).name]);
            if exist('pcaica','dir');
                cd('pcaica');
                if exist('allROIs.mat')>0;
                    load('allROIs.mat');
                    N=size(roi,3);
                    if N>1;
                        %dfof left
                        cd([p2 '\' dfolders(nn).name '\pcaica\cueAnalysisNew_dfof\newScoreShuffleTemplate\Left']);
                        load('cueCells.mat');
                        thresh=thsLDfofOld;
                        cueCellsUniThresh=[];
                        cueCellsUniThresh.useDfofaverage=cueCells.useDfofaverage;
                        cueCellsUniThresh.useCOMs=cueCells.useCOMs;
                        cueCellsUniThresh.useIdx=cueCells.useIdx;
                        cueCellsUniThresh.cueTemplate=cueCells.cueTemplate;
                        cueCellsUniThresh.maxLag=cueCells.maxLag;
                        cueCellsUniThresh.threshL=thresh;
                        cueCellsUniThresh.shuffleTemplates=cueCells.shuffleTemplates;
                        cueCellsUniThresh.shuffleScores=cueCells.shuffleScores;
                        cueCellsUniThresh.realScores=cueCells.realScores;
                        cueCellsUniThresh.realLags=cueCells.realLags;
                        cueCellsUniThresh.cueCellInUseIdx=find(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellRealIdx=cueCells.useIdx(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellScores=cueCells.realScores(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellLags=cueCells.realLags(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellDfofAvg=cueCells.useDfofaverage(:,cueCells.realScores>thresh);
                        save('cueCellsUniThresh.mat','cueCellsUniThresh');
                        per(end+1,1)=length(find(cueCells.realScores>thsLDfofOld))/N;
                        
                        %dfof right
                        cd([p2 '\' dfolders(nn).name '\pcaica\cueAnalysisNew_dfof\newScoreShuffleTemplate\Right']);
                        load('cueCells.mat');
                        thresh=thsRDfofOld;
                        cueCellsUniThresh=[];
                        cueCellsUniThresh.useDfofaverage=cueCells.useDfofaverage;
                        cueCellsUniThresh.useCOMs=cueCells.useCOMs;
                        cueCellsUniThresh.useIdx=cueCells.useIdx;
                        cueCellsUniThresh.cueTemplate=cueCells.cueTemplate;
                        cueCellsUniThresh.maxLag=cueCells.maxLag;
                        cueCellsUniThresh.threshL=thresh;
                        cueCellsUniThresh.shuffleTemplates=cueCells.shuffleTemplates;
                        cueCellsUniThresh.shuffleScores=cueCells.shuffleScores;
                        cueCellsUniThresh.realScores=cueCells.realScores;
                        cueCellsUniThresh.realLags=cueCells.realLags;
                        cueCellsUniThresh.cueCellInUseIdx=find(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellRealIdx=cueCells.useIdx(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellScores=cueCells.realScores(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellLags=cueCells.realLags(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellDfofAvg=cueCells.useDfofaverage(:,cueCells.realScores>thresh);
                        save('cueCellsUniThresh.mat','cueCellsUniThresh');
                        per(end,2)=length(find(cueCells.realScores>thsRDfofOld))/N;
                        
                        %sig Left
                        cd([p2 '\' dfolders(nn).name '\pcaica\cueAnalysisNew_sig\newScoreShuffleTemplate\Left']);
                        load('cueCells.mat');
                        thresh=thsLSigOld;
                        cueCellsUniThresh=[];
                        cueCellsUniThresh.useDfofaverage=cueCells.useDfofaverage;
                        cueCellsUniThresh.useCOMs=cueCells.useCOMs;
                        cueCellsUniThresh.useIdx=cueCells.useIdx;
                        cueCellsUniThresh.cueTemplate=cueCells.cueTemplate;
                        cueCellsUniThresh.maxLag=cueCells.maxLag;
                        cueCellsUniThresh.threshL=thresh;
                        cueCellsUniThresh.shuffleTemplates=cueCells.shuffleTemplates;
                        cueCellsUniThresh.shuffleScores=cueCells.shuffleScores;
                        cueCellsUniThresh.realScores=cueCells.realScores;
                        cueCellsUniThresh.realLags=cueCells.realLags;
                        cueCellsUniThresh.cueCellInUseIdx=find(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellRealIdx=cueCells.useIdx(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellScores=cueCells.realScores(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellLags=cueCells.realLags(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellDfofAvg=cueCells.useDfofaverage(:,cueCells.realScores>thresh);
                        save('cueCellsUniThresh.mat','cueCellsUniThresh');
                        per(end,3)=length(find(cueCells.realScores>thsLSigOld))/N;
                        
                        %sig right
                        cd([p2 '\' dfolders(nn).name '\pcaica\cueAnalysisNew_sig\newScoreShuffleTemplate\Right']);
                        load('cueCells.mat');
                        thresh=thsRSigOld;
                        cueCellsUniThresh=[];
                        cueCellsUniThresh.useDfofaverage=cueCells.useDfofaverage;
                        cueCellsUniThresh.useCOMs=cueCells.useCOMs;
                        cueCellsUniThresh.useIdx=cueCells.useIdx;
                        cueCellsUniThresh.cueTemplate=cueCells.cueTemplate;
                        cueCellsUniThresh.maxLag=cueCells.maxLag;
                        cueCellsUniThresh.threshL=thresh;
                        cueCellsUniThresh.shuffleTemplates=cueCells.shuffleTemplates;
                        cueCellsUniThresh.shuffleScores=cueCells.shuffleScores;
                        cueCellsUniThresh.realScores=cueCells.realScores;
                        cueCellsUniThresh.realLags=cueCells.realLags;
                        cueCellsUniThresh.cueCellInUseIdx=find(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellRealIdx=cueCells.useIdx(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellScores=cueCells.realScores(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellLags=cueCells.realLags(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellDfofAvg=cueCells.useDfofaverage(:,cueCells.realScores>thresh);
                        save('cueCellsUniThresh.mat','cueCellsUniThresh');
                        per(end,4)=length(find(cueCells.realScores>thsRSigOld))/N;
                    end
                end
            end
        end
    end
end
cd(p);


%%
p=pwd;

per=[];
load('folders.mat');
for n=1:length(folders);
    cd(folders{n});
    load('newEnvPath.mat');
    
    for idx=1:length(newEnvPath);
        disp(idx)
        cd(newEnvPath{idx});
        p2=pwd;
        d=dir;
        %remove non-folder files (isdir property is 0)
        dfolders = d([d(:).isdir]) ;
        % remove '.' and '..'
        dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));
        
        for nn=1:length(dfolders);
            cd([p2 '\' dfolders(nn).name]);
            if exist('pcaica','dir');
                cd('pcaica');
                if exist('allROIs.mat')>0;
                    load('allROIs.mat');
                    N=size(roi,3);
                    if N>1;
                        %dfof left
                        cd([p2 '\' dfolders(nn).name '\pcaica\cueAnalysisNew_dfof\newScoreShuffleTemplate\Left']);
                        load('cueCells.mat');
                        thresh=thsLDfofNew;
                        cueCellsUniThresh=[];
                        cueCellsUniThresh.useDfofaverage=cueCells.useDfofaverage;
                        cueCellsUniThresh.useCOMs=cueCells.useCOMs;
                        cueCellsUniThresh.useIdx=cueCells.useIdx;
                        cueCellsUniThresh.cueTemplate=cueCells.cueTemplate;
                        cueCellsUniThresh.maxLag=cueCells.maxLag;
                        cueCellsUniThresh.threshL=thresh;
                        cueCellsUniThresh.shuffleTemplates=cueCells.shuffleTemplates;
                        cueCellsUniThresh.shuffleScores=cueCells.shuffleScores;
                        cueCellsUniThresh.realScores=cueCells.realScores;
                        cueCellsUniThresh.realLags=cueCells.realLags;
                        cueCellsUniThresh.cueCellInUseIdx=find(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellRealIdx=cueCells.useIdx(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellScores=cueCells.realScores(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellLags=cueCells.realLags(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellDfofAvg=cueCells.useDfofaverage(:,cueCells.realScores>thresh);
                        save('cueCellsUniThresh.mat','cueCellsUniThresh');
                        per(end+1,1)=length(find(cueCells.realScores>thsLDfofNew))/N;
                        
                        %dfof right
                        cd([p2 '\' dfolders(nn).name '\pcaica\cueAnalysisNew_dfof\newScoreShuffleTemplate\Right']);
                        load('cueCells.mat');
                        thresh=thsRDfofNew;
                        cueCellsUniThresh=[];
                        cueCellsUniThresh.useDfofaverage=cueCells.useDfofaverage;
                        cueCellsUniThresh.useCOMs=cueCells.useCOMs;
                        cueCellsUniThresh.useIdx=cueCells.useIdx;
                        cueCellsUniThresh.cueTemplate=cueCells.cueTemplate;
                        cueCellsUniThresh.maxLag=cueCells.maxLag;
                        cueCellsUniThresh.threshL=thresh;
                        cueCellsUniThresh.shuffleTemplates=cueCells.shuffleTemplates;
                        cueCellsUniThresh.shuffleScores=cueCells.shuffleScores;
                        cueCellsUniThresh.realScores=cueCells.realScores;
                        cueCellsUniThresh.realLags=cueCells.realLags;
                        cueCellsUniThresh.cueCellInUseIdx=find(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellRealIdx=cueCells.useIdx(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellScores=cueCells.realScores(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellLags=cueCells.realLags(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellDfofAvg=cueCells.useDfofaverage(:,cueCells.realScores>thresh);
                        save('cueCellsUniThresh.mat','cueCellsUniThresh');
                        per(end,2)=length(find(cueCells.realScores>thsRDfofNew))/N;
                        
                        %sig Left
                        cd([p2 '\' dfolders(nn).name '\pcaica\cueAnalysisNew_sig\newScoreShuffleTemplate\Left']);
                        load('cueCells.mat');
                        thresh=thsLSigNew;
                        cueCellsUniThresh=[];
                        cueCellsUniThresh.useDfofaverage=cueCells.useDfofaverage;
                        cueCellsUniThresh.useCOMs=cueCells.useCOMs;
                        cueCellsUniThresh.useIdx=cueCells.useIdx;
                        cueCellsUniThresh.cueTemplate=cueCells.cueTemplate;
                        cueCellsUniThresh.maxLag=cueCells.maxLag;
                        cueCellsUniThresh.threshL=thresh;
                        cueCellsUniThresh.shuffleTemplates=cueCells.shuffleTemplates;
                        cueCellsUniThresh.shuffleScores=cueCells.shuffleScores;
                        cueCellsUniThresh.realScores=cueCells.realScores;
                        cueCellsUniThresh.realLags=cueCells.realLags;
                        cueCellsUniThresh.cueCellInUseIdx=find(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellRealIdx=cueCells.useIdx(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellScores=cueCells.realScores(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellLags=cueCells.realLags(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellDfofAvg=cueCells.useDfofaverage(:,cueCells.realScores>thresh);
                        save('cueCellsUniThresh.mat','cueCellsUniThresh');
                        per(end,3)=length(find(cueCells.realScores>thsLSigNew))/N;
                        
                        %sig right
                        cd([p2 '\' dfolders(nn).name '\pcaica\cueAnalysisNew_sig\newScoreShuffleTemplate\Right']);
                        load('cueCells.mat');
                        thresh=thsRSigNew;
                        cueCellsUniThresh=[];
                        cueCellsUniThresh.useDfofaverage=cueCells.useDfofaverage;
                        cueCellsUniThresh.useCOMs=cueCells.useCOMs;
                        cueCellsUniThresh.useIdx=cueCells.useIdx;
                        cueCellsUniThresh.cueTemplate=cueCells.cueTemplate;
                        cueCellsUniThresh.maxLag=cueCells.maxLag;
                        cueCellsUniThresh.threshL=thresh;
                        cueCellsUniThresh.shuffleTemplates=cueCells.shuffleTemplates;
                        cueCellsUniThresh.shuffleScores=cueCells.shuffleScores;
                        cueCellsUniThresh.realScores=cueCells.realScores;
                        cueCellsUniThresh.realLags=cueCells.realLags;
                        cueCellsUniThresh.cueCellInUseIdx=find(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellRealIdx=cueCells.useIdx(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellScores=cueCells.realScores(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellLags=cueCells.realLags(cueCells.realScores>thresh);
                        cueCellsUniThresh.cueCellDfofAvg=cueCells.useDfofaverage(:,cueCells.realScores>thresh);
                        save('cueCellsUniThresh.mat','cueCellsUniThresh');
                        per(end,4)=length(find(cueCells.realScores>thsRSigNew))/N;
                    end
                end
            end
        end
    end
end
cd(p);
