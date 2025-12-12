load('alignment.mat');
%this file only include aligned cells in the two FOVs. Put cell indices in
%the first and second columns for FOV1 and 2, respectively.
NCellFOV1=369;%total number of cells in the first fov 
NCellFOV2=343;%total number of cells in the second fov
[~,i]=sort(alignment(:,1));
A=alignment(i,:);
D1=setdiff([1:1:NCellFOV1],alignment(:,1));
D2=setdiff([1:1:NCellFOV2],alignment(:,2));

for n=1:length(D1);
    A(end+1,:)=[D1(n) 0];
end


for n=1:length(D2);
    A(end+1,:)=[0 D2(n)];
end

clear cell_registered_struct;

cell_registered_struct.cell_to_index_map=A;

save('cellRegisteredAlignManual.mat','cell_registered_struct');





