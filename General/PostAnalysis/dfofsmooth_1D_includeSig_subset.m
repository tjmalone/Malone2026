function [activity] = dfofsmooth_1D_includeSig_subset(imageStartNumber,imageEndNumber,...
    abfStartNumber,abfEndNumber,trackStart,trackEnd,binWidth,speedThreshold,dfof,abf)
%this code included dfof_sig just to use the criteria in Dombeck paper that
%at least 20% runs in field have sig trans

%this analysis is down based on dfof generated from the PCAICA code (manually selection of baseline) where the results of
%dfof fluorescence  are stored in the dfof
% keep only the images that were saved during clampex recording
% in this code, Faverage is averaged activity based on segTime. not real
% fluorescence


%% Process inputs

Fuseall = dfof(imageStartNumber:imageEndNumber,:);
nCells = size(Fuseall,2);

% process inputs
abfIms = abfStartNumber:abfEndNumber;
track = trackStart:binWidth:trackEnd;
gaussianWindow = gausswin(3*5/binWidth,1);

% speed at each index of abfIms
speed = diff(abf.y(abf.imageIndex(abfIms)));
fastEnough = [false ;speed>speedThreshold];

% get bin numbers
trackLength = trackEnd-trackStart;
[~,~,nbin] = histcounts(abf.y(abf.imageIndex(abfIms)),track);

% initialize output struct
activity = struct;

% calculate dfof smooth
dfofaverage = zeros(trackLength/binWidth,nCells);
for binNumber = 1:trackLength/binWidth
    dfofaverage(binNumber,:) = mean(Fuseall(nbin==binNumber&fastEnough,:),1,'omitnan');
end

% perform gaussian smoothing
for cellNumber = 1:nCells
    dfofaverage(:,cellNumber) = gaussianSmoothWithNan(dfofaverage(:,cellNumber),gaussianWindow);
end

% store dfof smooth
activity.dfofaverage = dfofaverage;


% calculate dfofM
dfofM_sig = getRunByRunActivityNoSmooth(1,nCells,imageStartNumber,imageEndNumber,...
    abfStartNumber,abfEndNumber,binWidth,trackStart,trackEnd,speedThreshold,dfof,abf);
activity.dfofM = dfofM_sig;

end

