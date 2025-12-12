%%

Alignment = [];


%%

clc


[~,I] = sort(Alignment(:,1));

alignment = Alignment(I,:);
    
disp(size(alignment,1))
disp(length(unique(alignment(:,1))))
disp(length(unique(alignment(:,2))))

[~, I] = unique(alignment(:,1), 'first');
x = 1:length(alignment);
x(I) = [];
disp(alignment(x,1))

[~, I] = unique(alignment(:,2), 'first');
x = 1:length(alignment);
x(I) = [];
disp(alignment(x,2))

%%

save('alignment.mat','alignment')

