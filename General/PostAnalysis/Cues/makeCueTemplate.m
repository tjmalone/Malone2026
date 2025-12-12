%%

clear

envName = '2mEnv3';
trackEnd = 200;
binSize = 5;
nBins = trackEnd/binSize;
tempR = zeros(nBins,1);
tempL = zeros(nBins,1);

R = [110:130];
L = [190:200];

tempR(unique(round(R/binSize))) = 1;
tempL(unique(round(L/binSize))) = 1;

tempRL = max(tempR,tempL);

tempR_2mEnv3 = tempR;
tempL_2mEnv3 = tempL;
tempRL_2mEnv3 = tempRL;

figure; hold on
plot(binSize:binSize:trackEnd,tempR)
plot(binSize:binSize:trackEnd,tempL)

save('tempR.mat','tempR')
save(['tempR_' envName '.mat'],['tempR_' envName])

save('tempL.mat','tempL')
save(['tempL_' envName '.mat'],['tempL_' envName])

save('tempRL.mat','tempRL')
save(['tempRL_' envName '.mat'],['tempRL_' envName])




