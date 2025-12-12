%%

clear; clc

nLocs = 1;
useIdx = {[1:3 5:11]};

% d = findSubF('pcaica',2);
d = findSubF('pcaica',2,[],0);

%%

for ii = 1:nLocs
    files = d(useIdx{ii});
    
    mkdir(['loc' num2str(ii)])
    save(['loc' num2str(ii) '\imageFiles.mat'],'files')
end

