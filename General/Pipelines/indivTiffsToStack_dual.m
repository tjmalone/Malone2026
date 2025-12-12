function indivTiffsToStack_dual(folderPath,outFolders,filter,NPerStack,dwnSample,lastIdx,cutPix)
% Converts individual tiff images to tiff stacks.

% Inputs:
%   (1) folder path: the folder containing all your individual images. For
%   example: C:\Users\guy5\Documents\ID20201105\20201201\ 
%   (2) outFolders: the folders to save tiff stacks. Must be cell array of
%   two folders
%   (3) filter: the key word shared by the name of all images. For example:
%   '*Ch2*'; 
%   (4) NPerStack: how many images per stack: for example, we normally use
%   1002 images per stack.
%   (5) dwnSample: the level of downsampling to be performed.
%   (6) cutPix: the first pixel on the x axis saved in tiff stack. Used to
%   crop tiffs with vertical artifacts

tic

% allows overwiting of tiff stacks
options.overwrite = true;

% default down sampling
if nargin<5 || isempty(dwnSample)
    dwnSample = 3;
end

% default crop
if nargin<7
    cutPix = 1;
end


%% Read file names

cd(folderPath);             % go to the folder containing individual tiffs
filenames = ls(filter);     % list all files with the name of the filter

% set max frames
if nargin>=6 && ~isempty(lastIdx) && lastIdx<=size(filenames,1)
    filenames = filenames(1:lastIdx,:);
end

filenamesOut = cell(1,2);   % file names for eaach Z plane
N = zeros(1,2);             % image number for each Z plane (after rolling downsapling)
NS = zeros(1,2);            % stack number for each Z plane

for ii = 1:2;
    filenamesOut{ii} = filenames(ii:2:end,:);
    N(ii) = size(filenamesOut{ii},1)-dwnSample+1;
    NS(ii) = ceil(N(ii)/NPerStack);
end

% return to outer folder (for speed)
cd ../


%% Generate tiff stacks

sz = size(loadtiff([folderPath '/' filenames(1,:)]));

% cycle through Z planes
for ii = 1:2
    
    % generate each tiff stack
    for n = 1:NS
        tStart = tic;
        disp(n)

        % idenitfy tiff names for current stack/plane (extra for downsampling)
        if n<NS(ii)
            filenamesUse=filenamesOut{ii}((n-1)*NPerStack+1:n*NPerStack+dwnSample-1,:);
        else
            filenamesUse=filenamesOut{ii}((n-1)*NPerStack+1:end,:);
        end
        
        % current file number
        useN = size(filenamesUse,1);
        
        % skip final stack if output will only have one frmae
        if useN<=dwnSample
            continue
        end
        
        % load tiffs
        data = zeros(sz(1),sz(2),useN);
        for m=1:useN;
            data(:,:,m) = loadtiff([folderPath '/' filenamesUse(m,:)]);
        end
        
        % crop all tiffs
        dataCrop = data(:,cutPix:end,:);
        
        % Matlab R2016a and later
%         dataMeans = movmean(data(:,cutPix:end,:),dwnSample,3);

        % Earlier than Matlab R2016a
        % calcualte rolling average
        dataMeans = zeros(sz(1),sz(2)+1-cutPix,useN-dwnSample+1);
        for jj = 1:useN-dwnSample+1
            dataMeans(:,:,jj) = mean(dataCrop(:,:,jj:jj+dwnSample-1),3);
        end
        
        % save tiff stack
        filenameStack=sprintf('%s_%03d.tif','tifStack/TS',n);
        saveastiff(int16(dataMeans),[outFolders{ii} '/' filenameStack],options);
        
        display(sprintf('The file was saved successfully. Elapsed time : %.3f s.', toc(tStart)));
    end
end
   
 toc   

end

