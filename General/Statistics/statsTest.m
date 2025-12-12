%%

load('A1.mat')
load('A2.mat')

[pA1,pM1] = anovaRM2W(A1,A2);
[pA2,pM2] = anovaRM2W_BH(A1,A2);

save('p_bonferroni.mat','pM1')
save('p_bonferroni_holm.mat','pM2')


%%

B1 = trapz(A1');
B2 = trapz(A2');

[~,p2] = ttest2(B1,B2);


%%

C1 = mean(A1,2);
C2 = mean(A2,2);

[~,p3] = ttest2(C1,C2);

[~,p4] = ttest2(A1(:),A2(:));









