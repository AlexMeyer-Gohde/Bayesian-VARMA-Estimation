 clear all; close all; home; format long g;  rng('shuffle');

numSeries = 1;
seriesCntrStart = 1;
numChains = 1;

consolidatedResults(numSeries, numChains) = struct( ...
    'pRaw', [], ...
    'qRaw', [], ...
    'pPostMax', [], ...
    'qPostMax', []);

[settings, priorsARMA, proposalsARMA] = getSettings();
% [y, AR_order, MA_order] = getData(settings);
% [modelinfo, parameters] = exampl1_extract(y);
datum = date;
if settings.saveProposals
    global PROPOSALS_GLOBAL;
    PROPOSALS_GLOBAL = getEmptyDrawStruct();
    PROPOSALS_GLOBAL(settings.draws) = getEmptyDrawStruct();
end;

for cntrData = seriesCntrStart:(seriesCntrStart + numSeries - 1)
    [y, AR_order, MA_order] = getData(settings);
    [modelinfo, parameters] = exampl1_extract(y);
    for cntrChains = 1:numChains
        clear states;
        
        disp(['Data set Nr.: ' num2str(cntrData) '; Chain Nr.: ' num2str(cntrChains)]);
        tic;
        [states, accepted, dropped, modelinfo, parameters] = bayesianInner(y, priorsARMA, proposalsARMA, settings, modelinfo, parameters);
        disp('Time elapsed for sampling'); toc;

        Filename = ['Results Bayesian_'];
        Filename = [Filename 'Monte Carlo Synthetic Philippe (' num2str(AR_order) ',' num2str(MA_order) ')_Sample Size 100_Draws ' num2str(settings.draws) '_Date ' datum '_Data Series ' num2str(cntrData) '_Chain ' num2str(cntrChains) ' Save Type 2.mat'];
        
        arParametersSeries = padcat(states.arParameters);
        maParametersSeries = padcat(states.maParameters);
        sigmaESeries = cat(1,states.sigmaEs);
        pSeries = cat(1,states.ps);
        qSeries = cat(1,states.qs);
        arPacsSeries = padcat(states.arPacs);
        maPacsSeries = padcat(states.maPacs);
        logPosteriorSeries = cat(1,states.logPosterior);
        
        if settings.saveProposals
            arParametersSeries_PROPOSAL = padcat(PROPOSALS_GLOBAL.arParameters);
            maParametersSeries_PROPOSAL = padcat(PROPOSALS_GLOBAL.maParameters);
            sigmaESeries_PROPOSAL = cat(1,PROPOSALS_GLOBAL.sigmaEs);
            pSeries_PROPOSAL = cat(1,PROPOSALS_GLOBAL.ps);
            qSeries_PROPOSAL = cat(1,PROPOSALS_GLOBAL.qs);
            arPacsSeries_PROPOSAL = padcat(PROPOSALS_GLOBAL.arPacs);
            maPacsSeries_PROPOSAL = padcat(PROPOSALS_GLOBAL.maPacs);
        end;

        if settings.saveProposals
            save(Filename,'arParametersSeries', 'maParametersSeries', 'sigmaESeries', 'pSeries', 'qSeries', 'arPacsSeries','maPacsSeries',...
                'arParametersSeries_PROPOSAL', 'maParametersSeries_PROPOSAL', 'sigmaESeries_PROPOSAL', 'pSeries_PROPOSAL', 'qSeries_PROPOSAL', 'arPacsSeries_PROPOSAL','maPacsSeries_PROPOSAL',...
                'logPosteriorSeries', 'dropped', 'accepted', 'y');
        else            
            save(Filename,'arParametersSeries', 'maParametersSeries', 'sigmaESeries', 'pSeries', 'qSeries', 'arPacsSeries',...
            'maPacsSeries', 'logPosteriorSeries', 'dropped', 'accepted', 'y');
        end;
                
        consolidatedResults(cntrData, cntrChains).pRaw = pSeries;
        consolidatedResults(cntrData, cntrChains).qRaw = qSeries;
        
        %Get Posterior Max in p,q        
        x = 0:1:settings.pMax;
        z = 0:1:settings.qMax;

        bintest = cell(1);
        bintest{1} = x;
        bintest{2} = z;

        [nelements, centers] = hist3([consolidatedResults(cntrData, cntrChains).pRaw(settings.burnIn:end) consolidatedResults(cntrData, cntrChains).qRaw(settings.burnIn:end)],'Edges',bintest);
        [maxPQ, ind] = max(nelements(:));
        [m,n] = ind2sub(size(nelements),ind)
        
        consolidatedResults(cntrData, cntrChains).pPostMax = m-1;
        consolidatedResults(cntrData, cntrChains).qPostMax = n-1;
        if mod(cntrData,10) == 0
            save('consolidated_Temp','consolidatedResults','-v7.3');
        end;        
    end;
end;

save(['Results_Monte Carlo Synthetic Philippe (3,2)_Consolidated_SampleSize 100_Draws 1500000_Date ' datum],'consolidatedResults','-v7.3');