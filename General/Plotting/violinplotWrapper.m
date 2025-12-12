function [h,xPositions] = violinplotWrapper(data,colors,indON,varargin)
%% violinplotWrapper
% A wrapper for the violinplot function by Bastian Bechtold to allow for
% additional input formats.
%
% Inputs:
%       X - Category names for x axis. If data is one-dimensional array, X
%           corresponds to the label of each bar. If data is
%           two-dimensional, X corresponds to the name of each group (each
%           row).
%       data - input data. Must be a 2D cell array of 1D numerical arrays.
%       colors - defines the colors for each violin. Must be cell array of
%           numerical triplets (0 to 1). Can be defined for each group or
%           for each bar. If defined for each bar, must be the same size as
%           data.
%       indON - whether to show individual data points. Can be true/false
%           or 0/1. Default: false
%       varargin - arguments fed directly to violinplotBB. Most arguments
%           should be valid, but some such as group order may lead to
%           errors.
%
% Taylor Malone, Yi Gu Lab, 2025
%

% supress legend warnings
warning('off','MATLAB:legend:IgnoringExtraEntries')

% get data sizes
nRows = size(data,1);
nCols = size(data,2);

% preocess indON
if nargin<4 || isempty(indON)
    indON = false;
else
    indON = logical(indON);
end

% set violin labels
vLabs = cellfun(@char,num2cell(string(1:nRows*nCols)),'UniformOutput',false);

% process color inputs
if nargin>2 && ~isempty(colors) && iscell(colors)
    useColor = true;

    % reshape color
    if sum(size(data)>1)==1
        if ~isequal(size(colors),size(data))
            colors = colors';
        end
    end

    % expand color rows
    if size(colors,1)==1
        colors = repmat(colors,nRows,1);
    elseif size(colors,1)~=nRows
        error('Invalid color row number')
    end

    % expand color columns
    if size(colors,2)==1
        colors = repmat(colors,1,nCols);
    elseif size(colors,2)<nCols
        error('Invalid color column number')
    end
else
    useColor = false;
end

% reformat data
violinData = [];
violinCats = {};
violinColors = [];
for rr = 1:nRows
    for cc = 1:nCols
        % generate data
        violinData = cat(1,violinData,data{rr,cc});

        % generate data categories
        szData = length(data{rr,cc});
        violinCats = cat(1,violinCats,repmat(vLabs((rr-1)*nCols+cc),szData,1));

        % generate category colors
        if useColor
            violinColors = cat(1,violinColors,colors{rr,cc});
        end
    end
end

% calculate violin positions and widths (groups are rows, bars are columns)
if nRows>1
    groupWidth = min(0.8,nCols/(nCols+1.5));
    barWidth = groupWidth/nCols;

    xPositions = zeros(nRows,nCols);
    for rr = 1:nRows
        center = rr;
        offsets = ((1:nCols)-(nCols+1)/2)*barWidth;
        xPositions(rr,:) = center + offsets;
    end

    % apply gap and half width
    halfBarWidth = barWidth/2*0.8;
else
    xPositions = 1:nCols;
    halfBarWidth = 0.3;
end

% reshape position and color data
trpPositions = xPositions';
usePositions = trpPositions(:);

% define median marker size
MedianMarkerSize = 36*halfBarWidth/0.3;

% update varargin to process color
if useColor
    varargin = [varargin,'ViolinColor',violinColors];
end

h = violinplotBB(violinData,violinCats,usePositions,'ShowData',indON,...
    'Width',halfBarWidth,'MedianMarkerSize',MedianMarkerSize,...
    'GroupOrder',vLabs,varargin{:});

xlim([min(usePositions)-0.5, max(usePositions)+0.5])


end

