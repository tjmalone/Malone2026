function s = surferror(X,Y,Z,err,varargin)
% surferror draws surf plot with errorbar
% input
%   X : x coordinate (mxn numeric matrix or n vector)
%   Y : y coordinate (mxn numeric matrix or m vector)
%   Z : heights (mxn numeric matrix)
%   err : errorbar heights (mxn numeric matrix)
%
%   Z and err should be same size
%   optional name-value pair arguments of surf are available (see surf in MATLAB doc)
%   2019-12-16 Hyeokjin Cho
%
%   1.0.2 supports vectorized coordinate input
%   2019-12-21 Hyeokjin Cho
    if nargin == 0 % show demo
        surferrorExample
        return
    end
    
    % exceptions
    if isvector(X) && isvector(Y)
        if length(X) == size(Z,2) && length(Y) == size(Z,1)
            [X,Y] = meshgrid(X,Y);
        else
            error('Data dimensions must agree')
        end
    else
        if ~all(size(X)==size(Y) & size(X) == size(Z))
            error('Data dimensions must agree')
        end
    end
    
    s = surf(X,Y,Z,varargin{:});
    originalNextplotsetting = get(gca,'nextplot');
    hold on
    % positive error
    quiver3(X,Y,Z,zeros(size(X)),zeros(size(Y)),err,'showarrowhead','off','autoscale','off','color','k');
    % negative error
    quiver3(X,Y,Z,zeros(size(X)),zeros(size(Y)),-err,'showarrowhead','off','autoscale','off','color','k');
    set(gca,'nextplot',originalNextplotsetting);
    shading interp
end
function surferrorExample
% make xy grid of x:-3~3, y:-3~3 with spacing 1
xrange = -3:3;
yrange = -3:3;
[X,Y] = meshgrid(xrange,yrange);
% bring matlab logo size of 7x7
Z = membrane(1,3);
% set errorbar height - in this example, we just use random number
err = rand(length(yrange),length(xrange));
% plot surferror
surferror(X,Y,Z,err);
end