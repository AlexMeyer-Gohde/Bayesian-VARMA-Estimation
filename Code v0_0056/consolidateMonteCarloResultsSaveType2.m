%Consolidate Results from Monte Carlo Experiments Anne Philippe
clear all; close all;

% availSeries = [1:7 11:14 21:26 31:33 41:43 51:53 61:64 71:72 81:83 91:93 101:102];

settings.burnIn = 1000000;
settings.pMax = 10;
settings.qMax = 10;


fileList = dir('Results Bayesian*');

consolidatedResults(length(fileList)) = struct( ...
    'pRaw', [], ...
    'qRaw', [], ...
    'pPostMax', [], ...
    'qPostMax', [],...
    'propPosteriorMass', []);


for cntrRecData = 1:length(fileList)
    for cntrRecChains = 1:numChains        
        disp(['Data set Nr.: ' num2str(cntrRecData) '; Chain Nr.: ' num2str(cntrRecChains)]);
        
        tic;
%         Filename = ['Results Bayesian_Monte Carlo Synthetic Philippe (3,2)_Sample Size 100_Draws 1500000_Date 'num2str(actDate) '-Jan-2014_Data Series ' num2str(cntrRecData) '_Chain ' num2str(cntrRecChains) ' Save Type 2.mat'];
        
        load(fileList(cntrRecData).name);
        
        tmpTime = toc;
        disp(['Finished loading file after ' num2str(tmpTime) ' seconds...']);        
        
        consolidatedResults(cntrRecData).pRaw = pSeries;
        consolidatedResults(cntrRecData).qRaw = qSeries;
        
        %Get Posterior Max in p,q        
        x = 0:1:settings.pMax;
        z = 0:1:settings.qMax;

        bintest = cell(1);
        bintest{1} = x;
        bintest{2} = z;

        [nelements, centers] = hist3([consolidatedResults(cntrRecData).pRaw(settings.burnIn+1:end) consolidatedResults(cntrRecData).qRaw(settings.burnIn+1:end)],'Edges',bintest);
        [maxPQ, ind] = max(nelements(:));
        [m,n] = ind2sub(size(nelements),ind);
        
        consolidatedResults(cntrRecData).pPostMax = m-1;
        consolidatedResults(cntrRecData).qPostMax = n-1;
        
        consolidatedResults(cntrRecData).propPosteriorMass = max(nelements(:)) / length(consolidatedResults(cntrRecData).pRaw(settings.burnIn+1:end));
        
        %Get series for MA and AR Parameters for later calculations
%         consolidatedResults(cntrRecData, cntrRecChains).seriesAR = arParametersSeries;
%         consolidatedResults(cntrRecData, cntrRecChains).seriesMA = maParametersSeries;
%         
%         %Get series of logPosteriors
%         consolidatedResults(cntrRecData, cntrRecChains).seriesLogPosterior = logProposalSeries;
    end;
end;

save(['Results_Monte Carlo Synthetic Philippe (3,2)_Consolidated_SampleSize 100_Draws 2000000_Date 09-Jan-2014.mat'], 'consolidatedResults');