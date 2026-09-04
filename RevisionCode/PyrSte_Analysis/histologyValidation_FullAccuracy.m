
%% Set inputs

clear; close all; clc

% set base folder
p1 = 'Z:\labMembers\TM\AD_Project\Revision\Histology\revisionHistology';
cd(p1)

% find all validation selections
foldsHV = findSubF('HV',2,[],0);

% define use FOV
useFOV = [3:6 11:12 16];

nFOV = length(useFOV);

% ture parameters
trueTypes = {'pyr','ste'};
trueFields = {'P','S'};
trueColors = {'b','r'};

nT = length(trueTypes);
idents = [-1 1];

% putative parameters
putTypes = {'pyr','und','ste'};
putIdents = [-1 0 1];

nP = length(putTypes);

% buffer settings
buffs = 0:100;
nBuff = length(buffs);
svBuff = 11;

% multiplier to deal with comparison direction
cMult = [1 NaN -1];

% define global threshold
globThresh = 156;

% analysis types
biTypes = {'Individual','Global'};
biField = {'ind','glob'};
nB = length(biField);

% analysis types
anaTypes = {'True Positive Rate','Population Purity','False Positve Rate'};
anaFields = {'true','purity','false'};
nAna = length(anaTypes);


%% Calculate confusion matrices

% initialize accuracy array
dataAcc = struct();
dataAcc.ind.ident = cell(nFOV,1);
dataAcc.glob.ident = cell(nFOV,1);
dataAcc.ind.conf = nan(3,2,nFOV,nBuff);
dataAcc.glob.conf = nan(3,2,nFOV,nBuff);

for ff = 1:nFOV
    %% Load inputs
    curF = useFOV(ff);

    % load current roi data
    load([foldsHV{curF} '\roiData.mat'],'roiData')

    % find threshold value
    load([foldsHV{curF} '\indThresh.mat'],'splitThr')


    %% Calculate accuracy

    % current thresholds
    biThresh = [splitThr,globThresh];

    % loop through true cell types
    for tt = 1:nT
        % define current areas
        curAreas = roiData.area(roiData.identity==idents(tt));

        % loop through comparison types
        for bb = 1:nB
            curThresh = biThresh(bb);

            % initialize save identity
            svIdent = nan(size(curAreas));

            % loop through test cell types
            for pp = 1:nP

                % apply size threshold
                if ismember(pp,[1 3])
                    cM = cMult(pp);
                    curTestType = cM*curAreas < cM*curThresh-buffs;
                else
                    curTestType = curAreas>=curThresh-buffs & curAreas<=curThresh+buffs;
                end

                % store cell number
                dataAcc.(biField{bb}).conf(pp,tt,ff,:) = sum(curTestType,1);
            end
        end
    end


    %% Calculate putative identites

    % define current areas
    curAreasAll = roiData.area;

    % loop through comparison types
    for bb = 1:nB
        curThresh = biThresh(bb);

        % initialize save identity
        svIdent = nan(size(curAreasAll));

        % loop through test cell types
        for pp = 1:nP

            % apply size threshold
            if ismember(pp,[1 3])
                cM = cMult(pp);
                curTestType = cM*curAreasAll < cM*curThresh-buffs(svBuff);
            else
                curTestType = curAreasAll>=curThresh-buffs(svBuff) & curAreasAll<=curThresh+buffs(svBuff);
            end

            % store identity
            svIdent(curTestType==1) = putIdents(pp);
        end

        % store identity
        dataAcc.(biField{bb}).ident{ff} = svIdent;
    end
end


%% Calculate accuracies

for bb = 1:nB
    curConf = dataAcc.(biField{bb}).conf;

    % calculate accuracy rates
    trueP = squeeze(curConf(1,1,:,:)./(curConf(1,1,:,:)+curConf(3,1,:,:)));
    trueS = squeeze(curConf(3,2,:,:)./(curConf(3,2,:,:)+curConf(1,2,:,:)));
    purityP = squeeze(curConf(1,1,:,:)./(curConf(1,1,:,:)+curConf(1,2,:,:)));
    purityS = squeeze(curConf(3,2,:,:)./(curConf(3,2,:,:)+curConf(3,1,:,:)));
    falseP = squeeze(curConf(1,2,:,:)./(curConf(1,2,:,:)+curConf(3,2,:,:)));
    falseN = squeeze(curConf(3,1,:,:)./(curConf(3,1,:,:)+curConf(1,1,:,:)));

    % store rates
    dataAcc.(biField{bb}).P.true = trueP;
    dataAcc.(biField{bb}).S.true = trueS;
    dataAcc.(biField{bb}).P.purity = purityP;
    dataAcc.(biField{bb}).S.purity = purityS;
    dataAcc.(biField{bb}).P.false = falseP;
    dataAcc.(biField{bb}).S.false = falseN;
end

save('Z:\labMembers\TM\AD_Project\Revision\Histology\histologyValidation\histValAccuracy.mat','dataAcc')


%% Plot full accuracy comparison

% define plot settings
yLims = [70 100;80 100;0 30];
xLim = [0 30];

% plot with individual FOV
figure; tiledlayout(nB,nAna)


for bb = 1:nB
    for aa = 1:nAna
        % current panel
        nexttile((bb-1)*nAna+aa); hold on

        for tt = 1:nT
            % get current data
            plotData = dataAcc.(biField{bb}).(trueFields{tt}).(anaFields{aa})*100;

            % plot data
            semshade(plotData,0.3,trueColors{tt},buffs);
        end

        % set plot labels
        xlim(xLim)
        ylim(yLims(aa,:))
        xlabel('Buffer (\mum^2)')
        ylabel([anaTypes{aa} ' (%)'])
        title([biTypes{bb} ' Threshold: ' anaTypes{aa}])
        % legend(trueTypes,'Location','best')
    end
end


