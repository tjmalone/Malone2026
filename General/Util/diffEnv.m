function envDif = diffEnv(idxs,type)
% determines the difference between two enviroments. User will be prompted
% to select folders.
%

%%

if nargin<2
    type = 1;
end

% set path for analysis code
pth = 'D:\AnalysisCode\PostAnalysis\';

filesL = findSubF('tempL.mat',2,pth,idxs);
filesR = findSubF('tempR.mat',2,pth,idxs);

nF = length(filesL);
envL = cell(1,nF);
envR = cell(1,nF);

for i = 1:nF
    load(filesL{i})
    load(filesR{i})
    
    envL{i} = tempL;
    envR{i} = tempR;
end

if type==1
    envDif = envL{1}+envR{1};
elseif type==2
    envDif = envL{2}+envR{2};
elseif type==3
    envDif = abs(envL{2}-envL{1}) + abs(envR{2}-envR{1});
end

envDif(envDif>0) = 1;

end

