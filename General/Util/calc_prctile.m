function x = calc_prctile(data,value)
    
    rngP = 0:0.01:100;
    perc = prctile(data,0:0.01:100);
    [~, index] = min(abs(perc-value));
    x = rngP(index)/100;
    
end