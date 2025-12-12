function params = cellParams2D(allROIs,pixelSize)
%this caculate cell parameters using matlab function regionprops, including
%long axis, short axis, area.

%input:
%allROIs: 3d matrix. 3rd dimention is number of cells, first and second
%dimentions are FOV pixels
%pixelSize: how many um per pixel

areapx=[];%area
longAxispx=[];%long axis
shortAxispx=[];%short axis

for n=1:size(allROIs,3);
    A=allROIs(:,:,n);
    stats=regionprops(A,'Area','MajorAxisLength','MinorAxisLength');
    areapx(n,1)=stats.Area;
    longAxispx(n,1)=stats.MajorAxisLength;
    shortAxispx(n,1)=stats.MinorAxisLength;
end

longAxisum=longAxispx*pixelSize;
areaum=areapx*pixelSize*pixelSize;
shortAxisum=shortAxispx*pixelSize;

%cluster cell
%long axis
A=longAxisum;
S=kmeans(A,2);
%let the small number to be 1, the bigger cluster to be 2;
[~,i]=min(A);
if S(i)==1;
    S;
else
    S(S==1)=3;
    S(S==2)=1;
    S(S==3)=2;
end
S1=S;

%short axis
A=shortAxisum;
S=kmeans(A,2);
%let the small number to be 1, the bigger cluster to be 2;
[~,i]=min(A);
if S(i)==1;
    S;
else
    S(S==1)=3;
    S(S==2)=1;
    S(S==3)=2;
end
S2=S;

%area
A=areaum;
S=kmeans(A,2);
%let the small number to be 1, the bigger cluster to be 2;
[~,i]=min(A);
if S(i)==1;
    S;
else
    S(S==1)=3;
    S(S==2)=1;
    S(S==3)=2;
end
S3=S;

params.allROIs=allROIs;
params.pixelSize=pixelSize;
params.areapx=areapx;
params.longAxispx=longAxispx;
params.shortAxispx=shortAxispx;
params.areaum=areaum;
params.longAxisum=longAxisum;
params.shortAxisum=shortAxisum;
params.area_cluster=S3;
params.longAxis_cluster=S1;
params.shortAxis_cluster=S2;

end

