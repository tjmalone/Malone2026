function pth = findSubF(subF,subN,baseF,idxs)
%% findSubF.m
%
% Idenitifies all iterations of a given folder or file type, subN levels
% down from a given folder.
%
% Inputs:
%   (1) subF: folder or file format to find within subfolders
%
%   (2) subN: depth of subfolder layers to search
%
%   (3) baseF: base folder to begin search from. If left empty, will be set
%       to current folder
%

% define base folder
if nargin<3 || isempty(baseF)
    baseF = pwd;
end


%% Select copy folders

p = pwd;
cd(baseF)

% identify all subN-level folders
d = dir([repmat('*/',[1 subN]) '*']);

% valid folder flag
isF = false(size(d));

% identify specified folders or files only
for f = 1:length(d)
    if strcmp(d(f).name,subF)
        isF(f) = 1;
    end
end

% valid folders
foldAna = d(isF);

if nargin<4 || isempty(idxs)
    % display parent directory of all identified folders
    for f = 1:length(foldAna)
        fName = erase(foldAna(f).folder,[pwd '\']);
        
        fprintf('Folder #%d: %s\n', f, fName)
    end
    
    fprintf(newline)
    
    % prompt user to select from folders
    prmt = sprintf('\nSelect files/folders: ');
    sAna = input(prmt);
elseif length(idxs)==1 && idxs==0
    sAna = 1:length(foldAna);
else
    sAna = idxs;
end

nAna = length(sAna);
pth = cell(nAna,1);

% save path for selected folders
for f = 1:nAna
    pth{f} = [foldAna(sAna(f)).folder '\' subF '\'];
end

cd(p)

end