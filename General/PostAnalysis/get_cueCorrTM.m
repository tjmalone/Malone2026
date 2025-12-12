function [cueCorr,lag] = get_cueCorrTM(frate,template,lagbins)

%this is similar to the get_cuesscoreYGModified_lessLagsNew but it
%calcualtes the correlation

if all(frate==0)
    cueCorr = 0;
    lag = 0;
else
    [frate_corr,frate_lags] = xcorr(template,frate-mean(frate),lagbins,'coeff');

    [cueCorr,i] = max(frate_corr);
    lag = frate_lags(i);
end

end
