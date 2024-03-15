function saveStuff(y, states, priorsARMA, proposalsARMA, settings, AR_order, MA_order)
% Save results with proper Naming
Filename = ['Results Bayesian_'];
switch settings.useData
    case 0
        Filename = [Filename 'Synthetic_Order (' num2str(AR_order) ',' num2str(MA_order) ')_Draws ' num2str(settings.draws) '_ ' date ' Save Type 2.mat'];
    case 1
        Filename = [Filename 'Log USGDP GAP HPFilter_Draws ' num2str(settings.draws) '_' date ' Save Type 2.mat'];
    case 2
        Filename = [Filename 'Log Growth USGDP_Draws ' num2str(settings.draws) '_' date ' Save Type 2.mat'];
    case 3
        Filename = [Filename 'Log USGDP minus Linear Trend_Draws ' num2str(settings.draws) '_' date ' Save Type 2.mat'];       
    case 5
        Filename = [Filename 'Uhlig Exampl 1_Draws '  num2str(settings.draws) '_' date ' Save Type 2.mat'];
    case 6
        Filename = [Filename 'US GDP PC HP Detrended Cyclical '  num2str(settings.draws) '_' date ' Save Type 2.mat'];
end;

        arParametersSeries = padcat(states.arParameters);
        maParametersSeries = padcat(states.maParameters);
        sigmaESeries = cat(1,states.sigmaEs);
        pSeries = cat(1,states.ps);
        qSeries = cat(1,states.qs);
        arPacsSeries = padcat(states.arPacs);
        maPacsSeries = padcat(states.maPacs);
        logProposalSeries = cat(1,states.logProposal);

        save(Filename,'arParametersSeries', 'maParametersSeries', 'sigmaESeries', 'pSeries', 'qSeries', 'arPacsSeries',...
            'maPacsSeries', 'logProposalSeries', 'dropped', 'accepted');
        toc


% save(Filename,'-v7.3');
end