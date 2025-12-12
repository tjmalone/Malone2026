%% Plot heat map

clear; close all; clc

p1 = ['Z:\labMembers\KC\AD_Project\Histology\tau_RC\B3456\1.5mm\ManualOverlap\'...
    'V4_addedMissedCells\ByFOV\plotRCTBySex'];
cd(p1)

dataTypes = {'perPS','perTau'};
data = cell(1,2);

load('plotRCTBySex_tauNum.mat','groups')
data{1} = flipud(reshape(groups,2,2)');

load('plotRCTBySex_tauDen.mat','groups')
data{2} = flipud(reshape(groups,2,2)');

ttls = {'Percent of stellate or pyramidal cells that are tau+',...
    'Percent of tau+ cells that are stellate or pyramidal'};
normTypes = {'Raw','Norm'};

for ii = 1:length(data)
    curData = data{ii};
    
    for jj = 1:2
        if jj==2
            for kk = 1:2
                curMean = mean(curData{kk,1});
                for ff = 1:2
                    curData{kk,ff} = curData{kk,ff}/curMean;
                end
            end
        end


        morphData = cell(2,1);
        for kk = 1:2
            morphData{kk} = cat(1,curData{:,kk});
        end

        % define data
        curMap = cell(3,3,2);
        for kk = 1:numel(curMap)
            curMap{kk} = nan;
        end
        curMap(1,2:3,:) = cat(2,{0;0},morphData);
        curMap(2:3,2:3,:) = cat(3,{0,0;0,0},curData);

        % make heat map
        h = heatmapGridSplit(curMap,{'Sex','All','Female','Male'},...
            {'Morphology','All','ste','pyr'},-1,{2:3},{2:3},{'eastoutside'});
        sgtitle([ttls{ii} ' ' normTypes{jj}],'FontSize',18)

    end
end
