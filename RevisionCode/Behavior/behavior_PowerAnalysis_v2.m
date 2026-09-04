
% clear all variables
clear; close all; clc

% reset directory
cd('D:\AD_Project\Behavior')
p = pwd;

% load processed behavioral data
load('data/allMiceDataB5.mat')

% load genotypes and sexes
load('groupIDs.mat')

% define legend
nGroups = length(groups);
nSexes = 2;

% define plot colors
colors6 = {[0 0 1],[1 0 0];[102 105 255]/255,[255 102 102]/255;[0 0 175]/255,[175 0 0]/255};


%% Plot all distributions

% simulation parameters
alpha = 0.05;
useN = 2:85;
nN = length(useN);

% define fields, titles, and labels
fields = {'perSuccess','perSuccessTry','dprimeGlobal','numStops',...
    'avgVelMovingRaw','predLicking'};

plotFields = [1 3];
useDays = 5:14;

% calculate maximum number of plotting sessions
maxSessions = max(cellfun(@length,{B(:).typeAll}));

% define colormap
curCMap = customcolormap([0 0.1999 0.2 1],[0 0 1; 0 0 1; 0.3 0.3 1; 1 1 1]);
curCMap = customcolormap([0 0.2 1],[0 0 1; 0.5 0.5 1; 1 1 1],10);

% define base square
baseN = {[6 4],[4 7]};
baseColor = {[1 0 0],[1 0 1]};

figure; tiledlayout(length(plotFields),2)

% loop through all fields
for pf = 1:length(plotFields)
    % load current field data
    f = plotFields(pf);

    % initialize data matrices
    curData = zeros(length(B),maxSessions);

    % extract data for current field (pad data with nan)
    for ii = 1:length(B)
        curD = B(ii).(fields{f})(B(ii).typeAll);
        curD(end+1:maxSessions) = nan;
        curData(ii,:) = curD;
    end

    % loop through sexes
    for ss = 1:nSexes
        %% Plot current field

        % current sexes
        curSexes = sexes{ss+1};

        % sort data by group for current sex
        data = cell(1,2);
        for gg = 1:nGroups
            data{gg} = curData(intersect(groups{gg},curSexes),useDays);
        end

        %%

        % run anova
        [~,~,~,~,rm] = anovaRM2W_full_BH(data{1},data{2});
        curMiceN = cellfun(@(x) size(x,1),data);

        % run anova table
        ranovatbl = ranova(rm,'WithinModel','Time');

        % extract effect and error terms
        ss_effect = ranovatbl.SumSq(2);
        ss_error = ranovatbl.SumSq(3);

        df_num = ranovatbl.DF(2);
        df_den = ranovatbl.DF(3);

        % calculate effect size (including eta interpretation)
        % (https://www.sciencedirect.com/science/article/pii/S1747938X11000029)
        partial_eta_sq = ss_effect/(ss_effect+ss_error);

        % claculate Cohen's F
        % (https://www.sciencedirect.com/science/article/pii/S1747938X11000029)
        cohen_f = sqrt(partial_eta_sq/(1-partial_eta_sq));

        % calculate F critical (Matlab innate function)
        F_crit = finv(1-alpha,df_num,df_den);

        powerEstPH = nan(nN,nN);

        % loop through sample sizes
        for ii = 1:nN
            for jj = 1:nN
                % current sample size
                curN = [useN(ii) useN(jj)];
                curNTot = sum(curN);

                % re-calculate DF
                cur_df_den = curNTot-2;

                % calculate F critical (Matlab innate function)
                cur_F_crit = finv(1-alpha,df_num,cur_df_den);

                % non-centrality paramater
                % FAQ How is effect size used in power analysis? UCLA: Statistical Consulting Group.
                % https://stats.oarc.ucla.edu/other/mult-pkg/faq/general/effect-size-power/faqhow-is-effect-size-used-in-power-analysis/
                cur_lambda = curNTot*(cohen_f^2);

                % calculate new power (innate matlab function)
                powerEstPH(ii,jj) = 1-ncfcdf(cur_F_crit,df_num,cur_df_den,cur_lambda);
            end
        end


        %% Plot results
        nexttile(2*(pf-1)+ss); hold on

        imagesc('XData',useN,'YData',useN,'CData',powerEstPH)

        % set axes labels
        xlim([useN(1)-0.5 useN(end)+0.5])
        ylim([useN(1)-0.5 useN(end)+0.5])
        xlabel('N_{WT}')
        ylabel('N_{PS19}')
        axis('square')
        title([fields{f} ': ' sexIDs{ss+1}])
        set(gca,'FontSize',12)

        % set colormap
        colormap(curCMap)
        cb = colorbar;
        cb.Label.String = 'Analysis Power';
        cb.Label.Rotation = -90;
        clim([0 1])

        for ii = 1:2
            bX = baseN{ii}(1)-0.5;
            bY = baseN{ii}(2)-0.5;
            rectangle('Position',[bX,bY,1,1],...
                'EdgeColor',baseColor{ii},'LineWidth',1)
        end
    end
end


