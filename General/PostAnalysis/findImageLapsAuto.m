%% findImageLapsAuto

clear; close all; clc

p1 = 'D:\AD_Project\imagingData';
cd(p1)

d = findSubF('suite2p',4,[],0);

for ii = 1:length(d)
    try
        findImageLaps(d{ii},1);

    catch
        disp(d{ii})
    end
end

cd(p1)
