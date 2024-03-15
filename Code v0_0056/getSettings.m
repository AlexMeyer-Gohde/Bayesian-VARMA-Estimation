function [settings, priorsARMA, proposalsARMA] = getSettings()
    %Parameters of Estimation
%Draws: Set the number of draws to be obtained
%Thinning: If the value of this variable is x, only every xth iterations is
%   used as a draw. A value of 1 turns off thinning.
%settings.burnIn: Set the length of the burn in. The burn in is not affected by
%   thinning.
%UseData:
%   Switch to switch between different data
%   0: synthetic
%   1: HP-Filtered US GDP PC Log Deviation+
%   2: US GDP PC Log-Growthrate - avg(log-growthrate)
%   3: US Log-GPD minus linear trend
%pMax:
%   Highest AR-Order to be considered in estimation (for all processes)
%qMax:
%   Highest MA-Order to be considered in estimation (for all processes)
%processCount:
%   Number of ARMA processes
%useSamePriorProposal:
%   If set to 1, the same priors and proposal distributions are used for
%   all processes. Otherwise, it will be possible to define different priors and
%   proposal distributions for each process.

settings.draws = 4000000
settings.thinning = 1;
settings.useData = 11;
settings.pMax = 10;
settings.qMax = 10;
settings.processCount = 1;
settings.burnIn = 1000000;
settings.useSamePriorProposal = 1;
settings.saveProposals = 1;

%Set Lambda to zero if you do not want to use HP Filtering in Likelihood
%Function
settings.HPLambda = 1600;

%Set likelihood function. only relevant if useSolver == 0. 
%   0: Spectral Based for DSGE
%   1: Kalman Filter based (experimental)
%   2: Approximate Likelihood for ARMA following Hamilton (experimental)
%   3: Compare output for tinkering
%   4: Kalman Alternative
%   5: Kalman One-Sided HP
%   6: Kalman No HP
settings.likelihoodFunction = 5;

%Switch use of dsge_solver on and off. Set to 0 if you want to estimate a
%normal ARMA process. If set to 1, spectral likelihood will be used
settings.useSolver = 1;

%Setup Priors. Array of struct
priorsARMA(1) = struct();

priorsARMA(1).priorARParam1 = 0;
priorsARMA(1).priorARParam2 = 0.5;
priorsARMA(1).priorARParam3 = 0/0;
priorsARMA(1).priorARParam4 = 0/0;
priorsARMA(1).priorAR = @(x) truncatedNormalImproperPrior(x,priorsARMA(1).priorARParam1,priorsARMA(1).priorARParam2);
% priorsARMA(1).priorAR = @(x) normpdf(x,priorsARMA(1).priorARParam1,priorsARMA(1).priorARParam2);
% priorsARMA(1).priorAR = @(x) unifpdf(x,priorsARMA(1).priorARParam1,priorsARMA(1).priorARParam2);

priorsARMA(1).priorMAParam1 = 0;
priorsARMA(1).priorMAParam2 = 0.5;
priorsARMA(1).priorMAParam3 = 0/0;
priorsARMA(1).priorMAParam4 = 0/0;
priorsARMA(1).priorMA = @(x) truncatedNormalImproperPrior(x,priorsARMA(1).priorMAParam1,priorsARMA(1).priorMAParam2);
% priorsARMA(1).priorMA = @(x) unifpdf(x,priorsARMA(1).priorMAParam1,priorsARMA(1).priorMAParam2);
% priorsARMA(1).priorMA = @(x) normpdf(x,priorsARMA(1).priorMAParam1,priorsARMA(1).priorMAParam2);

priorsARMA(1).priorSigmaEParam1 = 1;
priorsARMA(1).priorSigmaEParam2 = 1;
priorsARMA(1).priorSigmaEParam3 = 0/0;
priorsARMA(1).priorSigmaEParam4 = 0/0;
priorsARMA(1).priorSigmaE = ...
    @(x) (x>0)*(priorsARMA(1).priorSigmaEParam2^priorsARMA(1).priorSigmaEParam1...
    / gamma(priorsARMA(1).priorSigmaEParam1) * x^(-priorsARMA(1).priorSigmaEParam1 - 1)...
    * exp(-priorsARMA(1).priorSigmaEParam2/x));

priorsARMA(1).priorPParam1 = settings.pMax;
priorsARMA(1).priorPParam2 = 0/0;
priorsARMA(1).priorPParam3 = 0/0;
priorsARMA(1).priorPParam4 = 0/0;
priorsARMA(1).priorP = @(x) unidpdf(x + 1,priorsARMA(1).priorPParam1 + 1);

priorsARMA(1).priorQParam1 = settings.qMax;
priorsARMA(1).priorQParam2 = 0/0;
priorsARMA(1).priorQParam3 = 0/0;
priorsARMA(1).priorQParam4 = 0/0;
priorsARMA(1).priorQ = @(x) unidpdf(x + 1,priorsARMA(1).priorQParam1 + 1);

%Setup proposal distributions for ARMA Coefficients and respective standard
%deviation. 
proposalsARMA(1) = struct();

proposalsARMA(1).proposalARParam1 = 0.05;
proposalsARMA(1).proposalARParam2 = 0/0;
proposalsARMA(1).proposalARParam3 = 0/0;
proposalsARMA(1).proposalARParam4 = 0/0;
proposalsARMA(1).proposalAR = @(mu) vectorizedRTNorm(-1,1,mu,ones(size(mu,1),1)*proposalsARMA(1).proposalARParam1);
proposalsARMA(1).evaluateProposalAR = @(x, mu) evaluateTruncatedNormalPDF(x,-1,1,mu,ones(size(mu,1),1)*proposalsARMA(1).proposalARParam1);
% proposalsARMA(1).proposalAR =  @(mu) normrnd(mu,ones(size(mu,1),1)*proposalsARMA(1).proposalARParam1);
% proposalsARMA(1).evaluateProposalAR = @(x, mu) normpdf(x, mu, ones(size(mu,1),1)*proposalsARMA(1).proposalARParam1);

proposalsARMA(1).proposalMAParam1 = 0.05;
proposalsARMA(1).proposalMAParam2 = 0/0;
proposalsARMA(1).proposalMAParam3 = 0/0;
proposalsARMA(1).proposalMAParam4 = 0/0;
proposalsARMA(1).proposalMA = @(mu) vectorizedRTNorm(-1,1,mu,ones(size(mu,1),1)*proposalsARMA(1).proposalARParam1);
proposalsARMA(1).evaluateProposalMA = @(x, mu) evaluateTruncatedNormalPDF(x,-1,1,mu,ones(size(mu,1),1)*proposalsARMA(1).proposalMAParam1);
% proposalsARMA(1).proposalMA =  @(mu) normrnd(mu,ones(size(mu,1),1)*proposalsARMA(1).proposalMAParam1);
% proposalsARMA(1).evaluateProposalMA = @(x, mu) normpdf(x, mu, ones(size(mu,1),1)*proposalsARMA(1).proposalMAParam1);

proposalsARMA(1).proposalSigmaEParam1 = 0.05;
proposalsARMA(1).proposalSigmaEParam2 = 0/0;
proposalsARMA(1).proposalSigmaEParam3 = 0/0;
proposalsARMA(1).proposalSigmaEParam4 = 0/0;
proposalsARMA(1).proposalSigmaE = @(mu) vectorizedRTNorm(0,1000,mu,ones(size(mu,1),1)*proposalsARMA(1).proposalSigmaEParam1);
proposalsARMA(1).evaluateProposalSigmaE = @(x, mu) evaluateTruncatedNormalPDF(x,0,1000,mu,ones(size(mu,1),1)*proposalsARMA(1).proposalSigmaEParam1);
% proposalsARMA(1).proposalSigmaE =  @(mu) normrnd(mu, proposalsARMA(1).proposalSigmaEParam1);
% proposalsARMA(1).evaluateProposalSigmaE = @(x, mu) normpdf(x, mu, ones(size(mu,1),1)*proposalsARMA(1).proposalSigmaEParam1);

%Proposals AR-Order
proposalsARMA(1).proposalPParam1 = settings.pMax;
proposalsARMA(1).proposalPParam2 = 2;
proposalsARMA(1).proposalPParam3 = 0/0;
proposalsARMA(1).proposalPParam4 = 0/0;
%Uniform Proposal: 
% proposalsARMA(1).proposalP =  @(x) unidrnd(proposalsARMA(1).proposalPParam1 + 1) - 1;

%Discretized Laplace (Troughton Goodsill, Ehler Brooks 2004)
%Initialize Discrete Laplace CDF
proposalsARMA(1).laplaceCDFP = discreteLaplaceCDF(proposalsARMA(1).proposalPParam1, proposalsARMA(1).proposalPParam2);
proposalsARMA(1).laplacePDFP = discreteLaplacePDF(proposalsARMA(1).proposalPParam1, proposalsARMA(1).proposalPParam2);
proposalsARMA(1).proposalP = @(x) sampleDiscreteLaplace(x, proposalsARMA(1).laplaceCDFP);
proposalsARMA(1).evaluateProposalP = @(x) evaluateDiscreteLaplacePDF(x(1),x(2),proposalsARMA(1).laplacePDFP);


%Proposals MA-Order
proposalsARMA(1).proposalQParam1 = settings.qMax;
proposalsARMA(1).proposalQParam2 = 2;
proposalsARMA(1).proposalQParam3 = 0/0;
proposalsARMA(1).proposalQParam4 = 0/0;
%Uniform Proposal: 
% proposalsARMA(1).proposalQ =  @(x) unidrnd(proposalsARMA(1).proposalQParam1 + 1) - 1;

%Discretized Laplace (Troughton Goodsill, Ehler Brooks 2004)
proposalsARMA(1).laplaceCDFQ = discreteLaplaceCDF(proposalsARMA(1).proposalQParam1, proposalsARMA(1).proposalQParam2);
proposalsARMA(1).laplacePDFQ = discreteLaplacePDF(proposalsARMA(1).proposalQParam1, proposalsARMA(1).proposalQParam2);
proposalsARMA(1).proposalQ = @(x) sampleDiscreteLaplace(x, proposalsARMA(1).laplaceCDFQ);
proposalsARMA(1).evaluateProposalQ = @(x) evaluateDiscreteLaplacePDF(x(1),x(2),proposalsARMA(1).laplacePDFQ);


%Make sure process count is consistent
if settings.useSamePriorProposal ~= 1
    if settings.processCount ~= size(priorsARMA,2)
        disp(['settings.processCount does not coincide with number of priors. Resetting to ' num2str(size(priorsARMA,2)) ]);
        settings.processCount = size(priorsARMA,2);
    end;
end;


end