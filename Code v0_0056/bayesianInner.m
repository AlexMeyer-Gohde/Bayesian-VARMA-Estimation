function [states, accepted, dropped, modelinfo, parameters] = bayesianInner(y, priorsARMA, proposalsARMA, settings, modelinfo, parameters)

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
global PROPOSALS_GLOBAL
clear states;
states(settings.draws) = getEmptyStateStruct();
% draws = getEmptyDrawStruct();
draw = getEmptyDrawStruct();

%set initial state
% states(1).arPacs = [0.5];
% states(1).maPacs = [0.5];
% states(1).sigmaEs = 1;
% states(1).ps = 1;
% states(1).qs = 1;

states(1).arPacs = [];
states(1).maPacs = [];
states(1).sigmaEs = 1;
states(1).ps = 0;
states(1).qs = 0;

if states(1).ps > 0
    states(1).arParameters = getARParametersFromPACs(states(1).arPacs, states(1).ps);
end;
if states(1).qs > 0
    states(1).maParameters = -getARParametersFromPACs(states(1).maPacs, states(1).qs);
end;
states(1).logProposal = log(0.0000000000000000000000005);
if settings.useSolver == 1
    states(1).logPosterior = evaluatePosterior(y,states(1),priorsARMA,settings, modelinfo, parameters);
else
    states(1).logPosterior = evaluatePosteriorARMA(y,states(1),priorsARMA,settings);
end;

accepted=0;
dropped = 0;
% acceptance = zeros(settings.draws,1);

% progressbar

%Iterate until settings.draws is reached
for i = 2:settings.draws
    [state, draw, modelinfo, parameters] = MH_Step_Inner(y, priorsARMA, proposalsARMA, states(i-1), settings, modelinfo, parameters);
    if settings.saveProposals
        PROPOSALS_GLOBAL(i) = draw;
    end;   
    
    accepted = accepted + draw.accepted;
    dropped = dropped + draw.drawDropped;
    states(i) = state; 
%     progressbar(i/settings.draws)
%     if mod(i,500000) == 0
%         save('temp.mat', '-v7.3');
%     end;
    if mod(i,5000) == 0
        disp(['Iteration ' num2str(i)]);
        disp(['Acceptance rate ' num2str(accepted/i) ]);
    end;
end;
end
