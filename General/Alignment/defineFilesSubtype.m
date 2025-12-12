%%

clear; clc

nLocs = 2;

fields = {'famLearn','novLearn','famRecall','novRecall'};

genotype = 'WT';
famLearn = 3;
novLearn = [5:14];
famRecall = [16];
novRecall = [17];

% d = findSubF('suite2p',3);
d = findSubF('suite2p',3,[],0);

%%

for ii = 1:nLocs
    files = struct();
    
    for ff = 1:length(fields)
        if ~isempty(eval(fields{ff}))
            
            files.(fields{ff}) = d(2*(eval(fields{ff})-1)+ii);
        else
            files.(fields{ff}) = {};
        end
    end
    files.learn = [files.famLearn;files.novLearn];
    files.all = [files.famLearn; files.novLearn; files.famRecall;...
        files.novRecall];
    files.genotype = genotype;
    
    mkdir(['loc' num2str(ii)])
    save(['loc' num2str(ii) '\imageFilesNew.mat'],'files')
end

