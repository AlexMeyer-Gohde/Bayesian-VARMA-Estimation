function [newState, draw, modelinfo, parameters] = MH_Step_Inner(y,  priorsARMA, proposalsARMA, oldState, settings, modelinfo, parameters)
    persistent alphaOld;
    
    newState = oldState;
    
    %Get new draw
    draw = getDrawARMA(oldState, proposalsARMA, settings);
    
    %Evaluate joint log posterior
    if settings.useSolver == 1
        [logPosterior, estimatedVariance, dropDraw, modelinfo, parameters] = evaluatePosterior(y,draw,priorsARMA,settings, modelinfo, parameters);
    else
        [logPosterior, estimatedVariance, dropDraw] = evaluatePosteriorARMA(y,draw,priorsARMA,settings);
    end;
    
    if ~dropDraw
        draw.logPosterior= logPosterior;
        draw.estimatedVariance = estimatedVariance;

        %Compute acceptance probability, ONLY ONE PROCESS!!!!!
        
        %This part is always present: (Note that it is assumed that
        %evaluation of proposal distributions for orders returns
        %log-probabilities
        
        alpha = draw.logPosterior - oldState.logPosterior ...
            - (proposalsARMA(1).evaluateProposalP([oldState.ps(1) draw.ps(1)]) + proposalsARMA(1).evaluateProposalQ([oldState.qs(1) draw.qs(1)])) ...
            + (proposalsARMA(1).evaluateProposalP([draw.ps(1) oldState.ps(1)]) + proposalsARMA(1).evaluateProposalQ([draw.qs(1) oldState.qs(1)]));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %               Add Proposal Ratios                %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        logContribution = evaluateLogProposalRatioARMA(oldState, draw, proposalsARMA);
        
        alpha = alpha + logContribution;
        
        if (alphaOld == alpha)
            if (exp(alpha) ~= 0)
                disp('something fishy is going on here');
            end;
        end;


        if isfinite(alpha)==0
            alpha=(-1/0);
        end
    else
        alpha = (-1/0);
    end;

    alphaOld = alpha;
    drawOld = draw;
    oldStateOld = oldState;
    
    draw.drawDropped = dropDraw;
   
    %Accept draw with probability alpha
    if rand <= exp(min(alpha,0))
       % Accept the candidate
        newState.logPosterior = draw.logPosterior; 
        newState.logProposal = draw.logProposal;
        if draw.ps > 0
            newState.arParameters = draw.arParameters;
            newState.arPacs = draw.arPacs;
        else
            newState.arParameters = 0/0;
            newState.arPacs = 0;
        end;
        if draw.qs > 0
            newState.maParameters = draw.maParameters;
            newState.maPacs = draw.maPacs;
        else
            newState.maParameters = 0/0;
            newState.maPacs = 0;
        end;
        newState.ps = draw.ps;
        newState.qs = draw.qs;
        newState.estimatedVariance = draw.estimatedVariance;
        draw.accepted = 1;                % Note the acceptance
        newState.sigmaEs = draw.sigmaEs;
    else
        %reject the candidate
        draw.accepted  = 0;                % Note the rejection
    end;
end

