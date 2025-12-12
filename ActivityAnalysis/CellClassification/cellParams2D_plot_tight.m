function cellParams2D_plot_tight(longAxisum,shortAxisum,areaum,pixelSize)
% plots cell parameters data including long axis, short axis, and area

% combine variables for loop
varL = {longAxisum,shortAxisum,areaum};
ttls = {'Long axis','Short axis','Area'};

L1 = 15.6;
L2 = 17.6;

L3 = 138;
L4 = 178;

figure

% plot ksdensity for each parameter
for ii = 1:length(varL)
    % define current parameter
    A = varL{ii};

    % define ks scale
    if ii==3
        ksSize = 3;
    else
        ksSize = pixelSize;
        % ksSize = 0.3;
    end

    % plot current ksdensity
    subplot(2,3,ii+3); hold on
    [p,x] = ksdensity(A,'width',ksSize);
    % [p,x] = ksdensity(A);

    if ii==1
        xline(L1,'r')
        xline(L2,'r')
    elseif ii==3
        xline(L3,'r')
        xline(L4,'r')
    end

    plot(x,p);
    xlim([min(x) max(x)]);

    % plot current scatter
    subplot(2,3,ii);
    % scatter(A,1:1:length(A),'k','MarkerEdgeAlpha',0.2);
    h = scatter(randsample(A,length(A)),1:1:length(A),'Marker', 'o', 'MarkerEdgeAlpha', 0.2);
    h.SizeData = h.SizeData/10;
    xlim([min(x) max(x)]);
    title(ttls{ii});
end

% % long axis by area
% subplot(2,5,[4:5 9:10])
% h = scatter(longAxisum,areaum,'Marker', 'o', 'MarkerEdgeAlpha', 0.5);
% h.SizeData = h.SizeData/10;
% xlabel('long axis (um)')
% ylabel('area (um^2)')
% axis('square')

end

