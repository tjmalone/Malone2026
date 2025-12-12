%%

clear; close all; clc

panel = 4;
yLims = [0 0.9];
xLims = [0.5 10.5];
breakLims = [0.2 0.4];

openfig('D:\AD_Project\imagingData\Figures\Figures_current\timecourse\inter_common-cell_timecourse.fig');

ax = nexttile(panel);

figNew = figure;
axNew = copyobj(ax,figNew);
set(axNew, 'Units', 'normalized', 'Position', [0.13 0.11 0.775 0.815]);  % default fill

xlim(xLims)
ylim(yLims)
breakyaxis(breakLims);
