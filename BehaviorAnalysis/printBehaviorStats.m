%%


fprintf(['\nTemporal Correlation:\nWT: p = %.2G, PS19: p = %.2G\nWT: r = %.2G, PS19: r = %.2G\n\n' ...
    'ANOVA:\np = %.2G\n\n'],...
    statsData.pCorr(1),statsData.pCorr(2),statsData.rCorr(1),statsData.rCorr(2),...
    statsData.pAnova(2,1))



