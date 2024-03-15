function [y, priorsARMA, proposalsARMA, states, settings, accepted, dropped, modelinfo, parameters] = doBayesian(pilotTuningDraws, acceptanceRateTarget)
    [settings, priorsARMA, proposalsARMA] = getSettings();
    [y, AR_order, MA_order] = getData(settings);
    
    if pilotTuningDraws
        func = @(x)bayesianWrapped(x,y,priorsARMA,proposalsARMA,settings,pilotTuningDraws, acceptanceRateTarget);
        options = optimset('Diagnostics','on', 'Display', 'iter', 'TolFun', 0.01, 'MaxFunEvals', 20);
        [opt,fval,exitflag] = fmincon(func,[0.05 0.05], [], [], [], [], [0 0], [10 10], [], options)
    else
        % Get modelinfo setup here!!!!!!
        [modelinfo, parameters] = exampl1_extract(y);
        
        [states, accepted, dropped, modelinfo, parameters] = bayesianInner(y, priorsARMA, proposalsARMA, settings, modelinfo, parameters);
        try
        Filename = ['Results Bayesian_'];
        Filename = [Filename 'Data RBC Model AR 0.9 -0.23 0.3 MA -0.4 0.6 -0.5 std(z) 0_712 Shock 250 Obs Synthetic Round 1.mat'];
        
        arParametersSeries = padcat(states.arParameters);
        maParametersSeries = padcat(states.maParameters);
        sigmaESeries = cat(1,states.sigmaEs);
        pSeries = cat(1,states.ps);
        qSeries = cat(1,states.qs);
        arPacsSeries = padcat(states.arPacs);
        maPacsSeries = padcat(states.maPacs);
        logProposalSeries = padcat(states.logProposal);

        save(Filename,'arParametersSeries', 'maParametersSeries', 'sigmaESeries', 'pSeries', 'qSeries', 'arPacsSeries',...
            'maPacsSeries', 'logProposalSeries', 'dropped', 'accepted', 'y');
        catch(evenmoreBS)
            disp('Matlab is not so good');
        end;
        try 
        saveStuff(y, states, priorsARMA, proposalsARMA, settings, AR_order, MA_order);
        catch(whateverhappened)
        disp('HDD full?');
        end;
        plotResults(states,settings, accepted, y);
    end;
end