function badRoi(bROI)

load('allROIs.mat','roi');
i=setdiff(1:size(roi,3),bROI);
roi=roi(:,:,i);
save('allROIsNew.mat','roi');

end