%Pass values in x. Three variables:
%x(1): Sigma Proposal AR and MA
%x(2): Sigma Proposal Sigma_e

function acceptanceRateDistance = bayesianWrapped(x, y, priorsARMA, proposalsARMA, settings, pilotTuningDraws, acceptanceRateTarget)

proposalsARMA(1).proposalARParam1 = x(1);
proposalsARMA(1).proposalMAParam1 = x(1);

proposalsARMA(1).proposalSigmaEParam1 = x(2);

settings.draws = pilotTuningDraws;

%Do things
[states, accepted] = bayesianInner(y, priorsARMA, proposalsARMA, settings);

acceptanceRateDistance = (accepted/pilotTuningDraws - acceptanceRateTarget)^2;
end