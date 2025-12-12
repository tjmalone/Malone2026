function h = plotCues(cueX,cueY,cueLvl,color,cueMin)
% plot cue templates as rectangles

if nargin<5
    cueMin = 0;
end

locs = contiguous(cueY,1);
barLocs = cueX(locs{1,2});
nCues = size(barLocs,1);

barX = barLocs(:,1);
barY = cueMin*ones(nCues,1);
barW = diff(barLocs,[],2);
barH = (cueLvl-cueMin)*ones(nCues,1);
pos = [barX, barY, barW, barH];

for ii = 1:nCues
    h = rectangle('position',pos(ii,:),'linestyle','-','edgecolor',color);
end





