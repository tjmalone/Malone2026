%% quantTauLayersKC_AllTimecourses
% For Timecourses A B C
%%%%%
% Quantifies tau for the timecourse
% Updates from Version 2
% - new threshold calculated by finding the xth percentile of the pixels in
% the young mice (2 and 3 month old)
    % - a threshold for each of the 3 batches
% - outputs results for 3 options
        % - Raw
        % - Background subtraction tau intensity
        % - Percent area
%%%%%

%% Part 1: Set baseline/backSub for each timecourse batch

clear; close all; clc

%%%%%%%%%%%%%%%%%%%Input
percentile = 99;
%%%%%%%%%%%%%%%%%%%%

anaCh = 1;
nLy = 3;

foldTimecourse = 'Z:\SNM\labMembers\KC\AD_Project\histology_KC1.5mm\TauAge';
cd(foldTimecourse)
load('pixelVars.mat')

% identify all timecourses
f = 'Timecourse*';
d = dir(f);
numTimecourse = size(d,1);
timecourseNames = {};
baselineInt = {};
analysisType = {'raw','backgroundSub','perArea'};

%loop through timecourses individually (A, B, C)
for i = 1:numTimecourse
    cd(foldTimecourse)
    timecourseNames{i} = d(i).name;
    cd(timecourseNames{i})
    load('mouseInfo.mat')
    foldBatchTime = pwd;
    eraseString2 = [foldBatchTime '\'];
    
    %find folders for each mouse
    roiFilename = 'RoiSetCut.zip';
    fUse = findSubF(roiFilename,1,[],0)';
    eraseString = '\RoiSetCut.zip\';

    %sve age, fUse, and geno for later
    age(age>=10) = 9;
    allAge{i} = age;
    allTrueAge{i} = trueAge;
    allGeno{i} = geno;
    allSex{i} = sex;
    allFUse{i} = fUse;

    %set baseline based on young mice
    baselineCount = 1;
    for f = 1:length(fUse)
        % TimecourseIdx: what timecourse batch (A,B,C)
        timecourseIdx{i}(f,1) = i;
        fUse2{i}{f,1} = erase(fUse{f},eraseString);
        miceName{i}{f,1} = erase(fUse2{i}{f,1},eraseString2);

        if baseline(f) == 1
            % Load Image
            cd(erase(fUse{f},roiFilename))
            curIm = tiffreadVolume('MEC.tif');
            
            % normalize image
            curIm = double(curIm(:,:,anaCh));
            curIm = curIm./65535; %highest possible intensity is 1 now
            
            % save normalized image
            imRaw{f} = curIm;
            
            % extract ROI border
            [rows,columns] = size(curIm,1:2);
            regs = extractImageJROI(rows,columns,roiFilename);
            
        
            % Analyze regions
            for ii = 1:size(regs,3)
                
                % calculate region area/intenstity
                rProps = regionprops(regs(:,:,ii),'Area','PixelIdxList');
                rAreaCur = rProps.Area*APP;
                
                % Determine baseline for Background Subtraction
                if ii == 2
                    curImInt = sort(curIm(rProps.PixelIdxList),'descend');
                    baselineInt{i}(baselineCount,1) = prctile(curImInt, percentile);
                    baselineCount = baselineCount + 1;
                end
            end
        end
    end
    backSub(i) = mean(baselineInt{i});
end


%% Part 2: Set up to loop through all mice in all batches
%Make folder for where the graphs will be saved
saveFold = {};
fold = 'Z:\SNM\labMembers\KC\AD_Project\Histology\AD_timecourse_Leica';
cd(fold)
mkdir('FinalVersion')
cd('FinalVersion')
foldVersion = pwd;

% folder where code is saved
baseFolder = 'D:\MATLAB\Code\AnalysisCode\Analysis\HistologyCode\TauInt\TimecourseAll';

% copy self to current directory
copyfile([baseFolder '\quantTauLayersKC_AllTimecourses.m'],'quantTauLayersKC_AllTimecourses.m');

%Make folder for each analysis type
save('backgroundSubtraction.mat','backSub')
for i = 1:length(analysisType)
    cd([foldVersion])
    mkdir(analysisType{i})
    cd([analysisType{i}])
    saveFold{i} = pwd;
    if i == 3
        mkdir('tauPosPixels')
        cd('tauPosPixels')
        tauPixelsFold = pwd;
    end
end

%Combine age, fUse, and geno
combAge = vertcat(allAge{:});
combTrueAge = vertcat(allTrueAge{:});
combGeno = vertcat(allGeno{:});
combFUse2 = vertcat(fUse2{:});
combMiceName = vertcat(miceName{:});
combTimecourseIdx = vertcat(timecourseIdx{:});
combSex = vertcat(allSex{:});
fN = length(combFUse2);
 
% convert age to ageIndex
ageIndex = [];
for i = 1:length(combAge)
    if combAge(i) <= 3
        ageIndex(i) = combAge(i) - 1;
    else
        ageIndex(i) = combAge(i) - 2;
    end
end

ageIndex = ageIndex';

nAge = max(ageIndex);

% Ages
catAge = [2, 3, 5, 6, 7, 8, 9];

% genotypes
% WT=1, AD=2
catGeno = {'WT','AD'};
nGeno = max(combGeno);

% layers
catLayer = {'1','2','3'};

imRaw = cell(fN,1);
roiSeg = cell(fN,1);

% mean tau intensity per mm^2
tauInt = zeros(fN,nLy);

yAxis = {'Tau intensity (A.U.)', '% Area Tau+'};
% yAxis = {'Tau intensity per area (A.U.)', '% Area Tau+'};

foldVersion2 = 'Z:\SNM\labMembers\KC\AD_Project\Histology\AD_timecourse_Leica\quantTayLayersKC_AllTimeCourses\TimecourseABC';
foldMiceInfo = {foldVersion2, foldVersion};
for i = 1:length(foldMiceInfo)
    cd(foldMiceInfo{i})
    save('mouseInfo.mat', 'allAge', 'allTrueAge', 'allGeno', 'allFUse')
    save('cats.mat','catAge','catGeno','catLayer')
    save('mouseInfoComb.mat','combAge','combTrueAge','ageIndex', 'combGeno', 'combFUse2', 'combTimecourseIdx','backSub','combMiceName','combSex')
end


%% Part 3: Loop through tau quantification for each analysis type
% timecourses are combined
%Testing
% f = 1;
rStart = 1;
%%%

% loop through analysisType (raw, backsub, perArea)
for r = rStart:length(analysisType)
    
    % Cycle through folders
    for f = 1:fN
        % Load Image
        
        cd(erase(combFUse2{f},roiFilename))
        % curIm = tiffreadVolume('MAX_MEC.tif');
        curIm = tiffreadVolume('MEC.tif');
        
        % normalize image
        curIm = double(curIm(:,:,anaCh));
        curImOrig = curIm./65535; %highest possible intensity is 1 now

        if r == 1
            curIm = curImOrig;
        else
            %background subtraction
            curIm = curImOrig-backSub(combTimecourseIdx(f));
            curIm(curIm<0) = 0;
            if r == 3
                curIm(curIm>0) = 1;
            end
        end
        
        % save normalized image
        imRaw{f} = curIm;
        
        % extract ROI border
        [rows,columns] = size(curIm,1:2);
        regs = extractImageJROI(rows,columns,roiFilename);
    
        % Analyze regions
        for ii = 1:size(regs,3)
            
            % calculate region area/intenstity
            rProps = regionprops(regs(:,:,ii),'Area','PixelIdxList');
            rAreaCur = rProps.Area*APP;
            
            % mean tau intensity
            curInt = mean(curIm(rProps.PixelIdxList));
            
            % mean tau intensity per mm^2
            if r == 3
                tauInt(f,ii) = curInt*100;
                yAxisType = 2;
            else
                % tauInt(f,ii) = curInt/(rAreaCur/1E6); %If you want tau
                % intensity per area
                tauInt(f,ii) = curInt;
                yAxisType = 1;
            end
        end

        %Visualize tau acculumation in all layers
        % if r == 3
        %     L2 = regs(:,:,2);
        %     BL2 = bwboundaries(L2);
        %     IdxL2 = cat(1,BL2{:});
        %     pixelsPos = find(curIm == 1);
        % 
        %     image1 = curImOrig;
        %     image1(sub2ind(size(image1),pixelsPos)) = 1;
        %     image1(sub2ind(size(image1),IdxL2(:,1),IdxL2(:,2))) = 1;
        % 
        %     image2 = curImOrig;
        %     image2(sub2ind(size(image2),pixelsPos)) = 1;
        %     image2(sub2ind(size(image2),IdxL2(:,1),IdxL2(:,2))) = 0;
        % 
        %     image3 = curImOrig;
        %     image3(sub2ind(size(image3),pixelsPos)) = 0;
        %     image3(sub2ind(size(image3),IdxL2(:,1),IdxL2(:,2))) = 0;
        % 
        %     imagePlain = curImOrig;
        % 
        %     curRGB = cat(3,image1,image2,image3);
        % 
        %     % withRois
        %     subplot(1,2,1);
        %     imshow(curRGB)
        % 
        %     %withoutRois
        %     subplot(1,2,2);
        %     imshow(imagePlain)
        % 
        %     % Link the x-axes and y-axes of both subplots
        %     ax1 = subplot(1,2,1);
        %     ax2 = subplot(1,2,2);
        %     linkaxes([ax1, ax2], 'xy');  % Link both x and y axes
        % 
        %     cd(tauPixelsFold)
        %     savefig([num2str(f) '_TauPositive' num2str(combAge(f)) '_' num2str(combGeno(f)) '.fig']);
        % end

    end

    % Collate results (for MEC layer 2)  
    tauData = cell(nGeno,nAge);
    tauMean = zeros(nGeno,nAge,3);
    tauSEM = zeros(nGeno,nAge,3);
    tauP = zeros(nAge,3);
    
    for ii = 1:size(tauInt,1)
        tauData{combGeno(ii),ageIndex(ii)}(end+1,:) = tauInt(ii,:);
    end
    
    for jj = 1:nAge
        for ii = 1:nGeno
            tauMean(ii,jj,:) = mean(tauData{ii,jj},1);
            tauSEM(ii,jj,:) = std(tauData{ii,jj},0,1)/sqrt(length(tauData{ii,jj}(:,1)));
        end
        
        for lay = 1:3
            [~,tauP(jj,lay)] = ttest2(tauData{1,jj}(:,lay),tauData{2,jj}(:,lay));
        end
    end
            
    cd(saveFold{r})
    save('tauDataGrouped.mat', 'tauInt', 'tauData', 'combAge')
    
    % Plot Results (Line chart: all 3 layers together)
    figure; hold on
    colors = {[.5 .5 1],[0 0 1],[0 0 .5];...
        [1 .5 .5],[1 0 0],[.5 0 0]};
    for ii = 1:nGeno
        for jj = 1:3
            errorbar(catAge,tauMean(ii,:,jj),tauSEM(ii,:,jj),...
                'Color',colors{ii,jj},...
                'DisplayName',[catGeno{ii} ': Layer ' catLayer{jj}])
        end
    end    
    ylabel(yAxis(yAxisType));
    xlabel 'Age (months)';
    legend('Location','Northwest');
    title(analysisType{r});
    savefig('quantTauAllLayers');

    % line chart by layer
    for ll = 1:nLy
        figure; hold on
        colors = {[.5 .5 1],[1 .5 .5]};
        for gg = 1:nGeno
            errorbar(catAge,tauMean(gg,:,ll),tauSEM(gg,:,ll),...
                'Color',colors{1,gg},...
                'DisplayName',[catGeno{gg} ': Layer ' catLayer{ll}])
        end
        titleString = sprintf('Intensity layer %g', ll);
        title([analysisType{r} ': ' titleString]);
        ylabel(yAxis(yAxisType));
        xlabel 'Age (months)';
        legend('Location','Northwest');
        hold off
    
        savefig(titleString);
    end
    
    
    % Seperate Data by layer
    for ll = 1:nLy
        layerData = cell(nGeno,nAge);
        layerDataWithAgeWT = [];
        layerDataWithAgeAD = [];
        countWT = 1;
        countAD = 1;
        for i = 1:nGeno
            for j = 1:nAge
                for k = 1:length(tauData{i,j}(:,ll))
                    layerData{i,j}(k,1) = tauData{i,j}(k,ll);
                    if i == 1
                        layerDataWithAgeWT(countWT,1) = catAge(j);
                        layerDataWithAgeWT(countWT,2) = layerData{i,j}(k,1);
                        countWT = countWT + 1;
                    else
                        layerDataWithAgeAD(countAD,1) = catAge(j);
                        layerDataWithAgeAD(countAD,2) = layerData{i,j}(k,1);
                        countAD = countAD + 1;
                    end
                end
            end
        end
    
        fileName = ['NewTauData' num2str(ll)];
        save([fileName '.mat'], 'layerData', 'layerDataWithAgeWT','layerDataWithAgeAD');
    end
   
    % Scatter for layer 2
    ll = 2;

    load('NewTauData2.mat')
    load('tauDataGrouped.mat')
       
    figure; hold on
    colors = {[1 0 0],[0 1 0],[0 0 1]}; 
    for i = 1:fN
        if combGeno(i) == 2
            scatter(combAge(i),tauInt(i,2),30,colors{1,combTimecourseIdx(i)});
        else
            scatter(combAge(i),tauInt(i,2),30,colors{1,combTimecourseIdx(i)},"filled");
        end
    end
    titleString = sprintf('Layer %g tau intensity', ll);
    title(titleString);
    ylabel(yAxis(yAxisType));
    xlabel 'Age (months)';
    % Annotation.LegendInformation.IconDisplayStyle = 'off';
    % legend('red: A','green: B','blue: C', 'open = AD' , 'filled = WT', 'Location', 'Northwest', 'AutoUpdate', 'off');
    % Manually make legend:
    % First, Batch
    annotation('textbox', [0.15, .9, 0.1, 0.1], 'String', {'red: A', 'green: B', 'blue: C'}, ...
           'FitBoxToText', 'on', 'BackgroundColor', 'white', ...
           'EdgeColor', 'black', 'FontSize', 10);
    % 2nd, geno
    annotation('textbox', [0.8, .9, 0.1, 0.1], 'String', {'open = AD', 'filled = WT'}, ...
       'FitBoxToText', 'on', 'BackgroundColor', 'white', ...
       'EdgeColor', 'black', 'FontSize', 10);
    hold off
    savefig(titleString);
    
    % Box plot for layer 2
    x = [1,2,3,4];
    xGroup = {'WT young', 'WT old', 'AD young', 'AD old'};
   
    %age cutoff: 6.5 month
    fourGroups = cell(1,length(xGroup));
    
    for ii = 1:length(layerDataWithAgeWT)
        fourGroups{1,1} = layerDataWithAgeWT(1:4,2);
        fourGroups{1,2} = layerDataWithAgeWT(5:12,2);
    end
    
    for ii = 1:length(layerDataWithAgeAD)
        fourGroups{1,3} = layerDataWithAgeAD(1:17,2);
        fourGroups{1,4} = layerDataWithAgeAD(18:32,2);
    end
    
    sigPair = [1,4;2,4;3,4];
    % sigPair = combvec(x,x)';
    
    figure;
    barGroup(xGroup,fourGroups,[],sigPair)
    ylabel(yAxis(yAxisType))
    xlabel ('Group')
    title([analysisType{r} ': Layer 2 tau intensity'])
    savefig(['youngVsOldBarGraph'])
    
    save('layer2YvsO.mat', 'fourGroups')

end

close all;