%%  plotAverageBySex
% plots reel/cal/tau percentage for B3456 comparing Males and Females
% normalized version

clear; close all; clc

p1 = 'Z:\labMembers\KC\AD_Project\Histology\tau_RC\FinalVersion\ByFOV\plotRCTBySex';
cd(p1)

% end up 4 groups
groupNames = {'Female','Male'};
legs = {'ste','pyr'};
xGroups = 1:length(groupNames);

types = {'tauDen','tauNum'};
useFiles = {'plotRCTBySex_tauDen.mat','plotRCTBySex_tauNum.mat'};
ttls = {'Percent of tau+ cells that are stellate or pyramidal',...
    'Percent of stellate or pyramidal cells that are tau+'};

pValReg = {};
pValNorm = {};

%% initialize stats parameters and array
outStats = {};
nUnits = 'FOV';
testNameA = 'two-tailed paired Students t-test';
testPairA = 1;
testNameB = 'two-tailed paired Students t-test';
testPairB = 0;
testMC = 0;
testLimitP = 0;

for ii = 1:length(useFiles)
    %%

    load(useFiles{ii})

    % reshape data
    groups2D = flipud(reshape(groups,2,2)');

    % no comparisons
    figure;
    subplot(1,2,1)
    pShow = [1 2;1 3;2 4; 3 4];
    [h,pValReg{ii}] = barGroup(groupNames,groups2D,'violin',[],pShow,'pair');
    ylabel ('Percent (%)')
    xlabel ('Type')
    title(['regular: p Males= ' num2str(pValReg{ii}(3,4))]);
    legend(h,legs,'Location','best')

    % perform full statistics
    statData = flipud(groups2D)';
    testCat = [ttls{ii} '_regular'];
    outStats(end+1,:) = [testCat ttestEffectSize(statData(1,:),statData(2,:),testNameA,nUnits,testPairA,testMC,testLimitP)];

    % no comparisons
    groupsNorm = cell(size(groups2D));
    for jj = 1:2
        for kk = 1:2
            groupsNorm{jj,kk} = groups2D{jj,kk}/mean(groups2D{jj,1});
        end
    end

    subplot(1,2,2)
    pShow = [1 2;1 3;2 4; 3 4];
    [h,pValNorm{ii}] = barGroup(groupNames,groupsNorm,'violin',[],pShow,'pair');
    ylabel ('Normalized percent')
    xlabel ('Type')
    title(['normalized: p Males= ' num2str(pValNorm{ii}(2,4))]);
    sgtitle(ttls{ii} )
    savefig(['plotRCTBySex_' types{ii}])

    statData = flipud(groupsNorm(:,2));
    testCat = [ttls{ii} 'norm'];
    outStats(end+1,:) = [testCat ttestEffectSize(statData(1,:),statData(2,:),testNameB,nUnits,testPairB,testMC,testLimitP)];

    legend(h,legs,'Location','best')
end

