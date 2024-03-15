clear; close all; home; format long g;  rng('shuffle'); 

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
%   all processes. Otherwise, it is possible to define different priors and
%   proposal distributions for each process.

settings.draws = 2000000;
settings.thinning = 1;
settings.useData = 1;
settings.pMax = 10;
settings.qMax = 10;
settings.processCount = 1;
settings.burnIn = 500000;
settings.useSamePriorProposal = 1;

%Setup Priors. Array of struct
priorsARMA(1) = struct();

priorsARMA(1).priorARParam1 = 0;
priorsARMA(1).priorARParam2 = 0.5;
priorsARMA(1).priorARParam3 = 0/0;
priorsARMA(1).priorARParam4 = 0/0;
priorsARMA(1).priorAR = @(x) normpdf(x,priorsARMA(1).priorARParam1,priorsARMA(1).priorARParam2);
% priorsARMA(1).priorAR = @(x) unifpdf(x,priorsARMA(1).priorARParam1,priorsARMA(1).priorARParam2);

priorsARMA(1).priorMAParam1 = 0;
priorsARMA(1).priorMAParam2 = 0.5;
priorsARMA(1).priorMAParam3 = 0/0;
priorsARMA(1).priorMAParam4 = 0/0;
% priorsARMA(1).priorMA = @(x) unifpdf(x,priorsARMA(1).priorMAParam1,priorsARMA(1).priorMAParam2);
priorsARMA(1).priorMA = @(x) normpdf(x,priorsARMA(1).priorMAParam1,priorsARMA(1).priorMAParam2);

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
proposalsARMA(1).proposalAR =  @(mu) normrnd(mu,ones(size(mu,1),1)*proposalsARMA(1).proposalARParam1);
proposalsARMA(1).evaluateProposalAR = @(x, mu) normpdf(x, mu, ones(size(mu,1),1)*proposalsARMA(1).proposalARParam1);

proposalsARMA(1).proposalMAParam1 = 0.05;
proposalsARMA(1).proposalMAParam2 = 0/0;
proposalsARMA(1).proposalMAParam3 = 0/0;
proposalsARMA(1).proposalMAParam4 = 0/0;
proposalsARMA(1).proposalMA =  @(mu) normrnd(mu,ones(size(mu,1),1)*proposalsARMA(1).proposalMAParam1);
proposalsARMA(1).evaluateProposalMA = @(x, mu) normpdf(x, mu, ones(size(mu,1),1)*proposalsARMA(1).proposalMAParam1);

proposalsARMA(1).proposalSigmaEParam1 = 0.05;
proposalsARMA(1).proposalSigmaEParam2 = 0/0;
proposalsARMA(1).proposalSigmaEParam3 = 0/0;
proposalsARMA(1).proposalSigmaEParam4 = 0/0;
proposalsARMA(1).proposalSigmaE =  @(mu) normrnd(mu, proposalsARMA(1).proposalSigmaEParam1);
proposalsARMA(1).evaluateProposalSigmaE = @(x, mu) normpdf(x, mu, ones(size(mu,1),1)*proposalsARMA(1).proposalSigmaEParam1);

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


switch settings.useData
    case 0
        %Generate AR
        sampleSize         =   200;
        sigma_sample        =   1.5;
        AR_coefficients     =   [0; 0; -0.75]
        MA_coefficients     =   [1; -1.5; 0.5625]

%         AR_coefficients     =   [-0.25; 0.23; -0.45]
%         MA_coefficients     =   [1; -0.4; 0.2]
%         MA_coefficients     =   [1; -0.4; 0.2; 0.6; -0.35]
%         MA_coefficients     =   [1]

        disp(['Roots AR Terms: ']);
        abs(roots([1; -AR_coefficients]'))
        disp(['Roots MA Terms: ']);
        abs(roots([MA_coefficients]'))
        AR_order = size(AR_coefficients,1);
        MA_order = size(MA_coefficients,1)-1;
        y = zeros(sampleSize+1000,1);        
        epsilon = normrnd(zeros(sampleSize+1000+ max(MA_order,AR_order)+1,1),sigma_sample);
        
        for AR_gen_step = (max(AR_order,MA_order)+1):sampleSize+1000+max(AR_order,MA_order)+1
            y(AR_gen_step) = transpose(y(AR_gen_step - AR_order:AR_gen_step - 1,1))*AR_coefficients(end:-1:1) +...
                transpose(MA_coefficients)*epsilon(AR_gen_step:-1:AR_gen_step-MA_order);
        end;
        y = y(AR_order+1000+1:end);
    case 1
        %Generate Pure MA
        sampleSize         =   200;
        sigma_sample        =   0.8;
        MA_coefficients     =   [1]
%         MA_pacs = [
%    -0.8049;
%    -0.4430;
%     0.6]

%         MA_coefficients =  [1; -getARParametersFromPACs(MA_pacs, 3)]

        disp(['Roots MA Terms: ']);
        abs(roots([MA_coefficients]'))
        MA_order = size(MA_coefficients,1)-1;
        y = zeros(sampleSize,1);        
        epsilon = normrnd(zeros(sampleSize+MA_order,1),sigma_sample);
        
        for gen_step = 1:sampleSize
            y(gen_step) = transpose(MA_coefficients)*epsilon(gen_step+ MA_order:-1:gen_step);
        end;   
        
    case 2 %Use Log Deviation
        [ y ] = getGDPData(1);
    case 3 %Use Log growth
        [ y ] = getGDPData(2);
    case 4 %Use Log GDP minus linear trend
        [ y ] = getGDPData(3);        
end;    


%Create struct to save state of chain
%Structure:
%.logPosterior:  Value of Log-Posterior at state
%.p: List of p's for each ARMA-process (AR-Order), each process one column
%.q: List of q's for each ARMA-process (MA-Order), each process one column
%.sigmaE: Shock variance for each ARMA-process
%.arParameters: Array with AR-Parameters for each
%ARMA-process in column vectors
%.maParameters: Array with MA-Parameters for each
%ARMA-process in column vectors
% currentState = struct('logPosterior',[], 'ps', [], 'qs', [], 'sigmaEs', [],...
%     'arParameters',     [],...
%     'maParameters',     [],...
%     'arPacs',           [],...
%     'maPacs',           []);

%initialize state structs
clear states;
states(settings.draws) = getEmptyStateStruct();
% draws = getEmptyDrawStruct();
draw = getEmptyDrawStruct();

%set initial state
states(1).arPacs = [0.2; 0.8];
states(1).maPacs = [0.8];
states(1).sigmaEs = 0.8;
states(1).ps = 2;
states(1).qs = 1;
states(1).arParameters = getARParametersFromPACs(states(1).arPacs, states(1).ps);
states(1).maParameters = -getARParametersFromPACs(states(1).maPacs, states(1).qs);
states(1).logProposal = 0.0000000000000000000000005;
states(1).logPosterior = evaluatePosterior(y,states(1),priorsARMA,settings);

accepted=0;
acceptance = zeros(settings.draws,1);

tic;
%Iterate until settings.draws is reached
for i = 2:settings.draws
    [state, draw] = MH_Step_Inner(y, priorsARMA, proposalsARMA, states(i-1), settings);
    accepted = accepted + draw.accepted;
    states(i) = state; 
    if mod(i,100) == 0
%         disp(['Iteration: ' num2str(i)]);
    end;
end
toc;
% Save results with proper Naming
Filename = ['Results Bayesian_'];
switch settings.useData
    case 0
        Filename = [Filename 'Synthetic_Order (' num2str(AR_order) ',' num2str(MA_order) ')_Draws ' num2str(settings.draws) '_SampleSize ' num2str(sampleSize) '_' date '.mat'];
    case 1
        Filename = [Filename 'Log USGDP GAP HPFilter_Draws ' num2str(settings.draws) '_' date '.mat'];
    case 2
        Filename = [Filename 'Log Growth USGDP_Draws ' num2str(settings.draws) '_' date '.mat'];
    case 3
        Filename = [Filename 'Log USGDP minus Linear Trend_Draws ' num2str(settings.draws) '_' date '.mat'];       
end;
    
save(Filename,'-v7.3');

plotResults(states,settings, accepted, y);