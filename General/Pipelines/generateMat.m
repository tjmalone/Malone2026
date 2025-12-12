function generateMat()
%% generateMat
% generates a .mat file for a set of tifs

%% Get file info

% list all motion corrected tif files
matrixPath = [pwd '\'];
fnMatrix =dir([matrixPath '*tif']);
matrixNumber=length(fnMatrix);
matrixAllPaths=cellfun(@(x)([matrixPath x]),{fnMatrix(1:matrixNumber).name},'uniformoutput',0);

%% Generate mat files

mkdir('matrix')

for n = 1:matrixNumber
    % load tif
    data = loadtiff(matrixAllPaths{n});

    % save mat file
    filenameData=sprintf('%s_%03d.mat','matrix/matrix',n);
    save(filenameData,'data')
end

