%% alignFix.m
%
% Cycles through FOV alignments to check for proper alignment. Allows user
% to manually align poorly aligned FOV pairs.
%

clear; close all; clc

p = pwd;

load('folders')
n = length(folders);

N = 10;
goodAlign = zeros(N-1,N-1);

screenSet = {[50 50 600 600],[700 50 600 600]};

%%

checkFile = 'Figures\Stage 2 - pre vs post alignment.fig';
checkFile2 = 'Figures\Stage 5 - projcetions - final registration.fig';

% checkFiles = findSubF(checkName,3,p,0);

% useIdx = [2 9 16 18 20];
useIdx = 1:n;

% wheter to only check alignments to a specified FOV
useAlign = 0;
alignIdx = 12;


for ii = 1:length(useIdx)-1
    % skip extra FOV pairs
    if useAlign && ii>alignIdx
        continue
    end

    % set iteration
    if ~useAlign || ii==alignIdx
        itrIdx = ii+1:length(useIdx);
    else
        itrIdx = alignIdx;
    end

    for jj = itrIdx



        curFolder = [num2str(useIdx(ii)) '\' num2str(useIdx(ii)) '_' num2str(useIdx(jj))];
        cd(curFolder)
        
        %%
        
        F = openfig(fullfile(p,curFolder,checkFile));
        
        clear set;
        
        usePlot = 5;
        H = figure(2);
        PLO = copyobj(F.Children(usePlot(1)),H);
        set(PLO,'Position',get(0,'DefaultAxesPosition'))
        set(H,'Position',screenSet{2})
        
        close(F);
                
        F = openfig(fullfile(p,curFolder,checkFile2));
        
        clear set;
        
        usePlot = 1;
        H = figure(3);
        PLO = copyobj(F.Children(usePlot(1)),H);
        set(PLO,'Position',get(0,'DefaultAxesPosition'))
        set(H,'Position',screenSet{1})
        axis('equal','off')
        close(F);
        
%         set(gcf, 'Position',screenSet);
        keep = input([num2str(ii) '-' num2str(jj) '. Good alignment?: ']);
       

        if keep==1
            goodAlign(ii,jj-1) = 1;
        end
        
        close all
        cd(p)
    end
end

cd(p)


%%

alignTran = goodAlign';

alignAll = [];

for ii = 1:N-1
    alignAll = [alignAll; alignTran(ii:end,ii)];
end

alignIdxs = zeros(0,2);
for ii = 1:N-1
    for jj= ii+1:N
        alignIdxs(end+1,:) = [ii,jj];
    end
end

save('alignCheckData.mat','alignAll','alignIdxs')