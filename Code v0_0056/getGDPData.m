function [ Output ] = getGDPData(switchey)

    GDPPC = load('GDPPC.mat');
    switch switchey
        case 1
            GDPPCHPFiltered = hpfilter(log(GDPPC.GDPPC),1600);
            Output = log(GDPPC.GDPPC) - GDPPCHPFiltered;
        case 2
            GDPLogGrowth = log(GDPPC.GDPPC(2:end)) - log(GDPPC.GDPPC(1:end-1));
            Output = GDPLogGrowth - mean(GDPLogGrowth);   
        case 3
            GDPLog = log(GDPPC.GDPPC);
            Output = detrend(GDPLog);
    end

