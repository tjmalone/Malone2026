function cellParams2D_plot(longAxisum,shortAxisum,areaum,pixelSize)
% plots cell parameters data including long axis, short axis, and area

% combine variables for loop
varL = {longAxisum,shortAxisum,areaum};
ttls = {'Long axis','Short axis','Area'};

figure

% plot ksdensity for each parameter
for ii = 1:length(varL)
    % define current parameter
    A = varL{ii};

    % define ks scale
    if ii==3
       ksSize = 10;
    else
       ksSize = pixelSize;
       ksSize = 0.4;
    end

    % plot current ksdensity
    subplot(2,5,ii+5);
    [p,x] = ksdensity(A,'width',ksSize);
    % [p,x] = ksdensity(A);

    plot(x,p);
    xlim([min(x) max(x)]);

    % plot current scatter
    subplot(2,5,ii);
    % scatter(A,1:1:length(A),'k','MarkerEdgeAlpha',0.2);
    h = scatter(randsample(A,length(A)),1:1:length(A),'Marker', 'o', 'MarkerEdgeAlpha', 0.5);
    h.SizeData = h.SizeData/10;
    xlim([min(x) max(x)]);
    title(ttls{ii});
end

% long axis by area
subplot(2,5,[4:5 9:10])
h = scatter(longAxisum,areaum,'Marker', 'o', 'MarkerEdgeAlpha', 0.5);
h.SizeData = h.SizeData/10;
xlabel('long axis (um)')
ylabel('area (um^2)')
axis('square')

end

