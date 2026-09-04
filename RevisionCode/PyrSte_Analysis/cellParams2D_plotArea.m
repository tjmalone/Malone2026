function cellParams2D_plotArea(areaum,ksSize)
% plots cell parameters data including area

% plot current ksdensity
[p,x] = ksdensity(areaum,'width',ksSize);

plot(x,p);
xlim([min(x) max(x)]);
title(['Area: width = ' num2str(ksSize)]);
end

