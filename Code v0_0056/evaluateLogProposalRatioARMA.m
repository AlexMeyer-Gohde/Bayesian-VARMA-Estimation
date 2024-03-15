function [logContribution] = evaluateLogProposalRatioARMA(oldState, draw, proposalsARMA)
           %AR Proposals        
        switch logical(true)
            case draw.ps > oldState.ps %draw.ps always > 0
                if oldState.ps > 0
                    muEnumerator =  draw.arPacs(1:oldState.ps);
                    muDenominator = [oldState.arPacs; zeros(draw.ps - oldState.ps,1)];
                    
                    logContributionAR = sum(log(proposalsARMA(1).evaluateProposalAR(oldState.arPacs, muEnumerator))) -...
                                        sum(log(proposalsARMA(1).evaluateProposalAR(draw.arPacs,     muDenominator)));
                else
                    muDenominator = zeros(size(draw.arPacs));                    
                    logContributionAR = - sum(log(proposalsARMA(1).evaluateProposalAR(draw.arPacs,   muDenominator)));
                end;
                
            case draw.ps == oldState.ps
                if draw.ps > 0 %then also oldState.ps > 0
                    muEnumerator =  draw.arPacs;
                    muDenominator = oldState.arPacs;
                    
                    logContributionAR = sum(log(proposalsARMA(1).evaluateProposalAR(oldState.arPacs, muEnumerator))) -...
                                        sum(log(proposalsARMA(1).evaluateProposalAR(draw.arPacs,     muDenominator)));
                else %draw.ps == oldState.ps == 0
                    logContributionAR = 0;
                end;
                
            case draw.ps < oldState.ps %thus oldState.ps > 0 in any case
                if draw.ps > 0 
                    muEnumerator =  [draw.arPacs; zeros(oldState.ps - draw.ps,1)];
                    muDenominator = oldState.arPacs(1:draw.ps);   
                    
                    logContributionAR = sum(log(proposalsARMA(1).evaluateProposalAR(oldState.arPacs, muEnumerator))) -...
                                        sum(log(proposalsARMA(1).evaluateProposalAR(draw.arPacs,     muDenominator)));
                else 
                    muEnumerator = zeros(size(oldState.arPacs));                    
                    logContributionAR = sum(log(proposalsARMA(1).evaluateProposalAR(oldState.arPacs, muEnumerator)));
                end;
        end;

        %MA Proposals        
        switch logical(true)
            case draw.qs > oldState.qs %draw.qs always > 0
                if oldState.qs > 0
                    muEnumerator =  draw.maPacs(1:oldState.qs);
                    muDenominator = [oldState.maPacs; zeros(draw.qs - oldState.qs,1)];
                    
                    logContributionMA = sum(log(proposalsARMA(1).evaluateProposalMA(oldState.maPacs, muEnumerator))) -...
                                        sum(log(proposalsARMA(1).evaluateProposalMA(draw.maPacs,     muDenominator)));
                else
                    muDenominator = zeros(size(draw.maPacs));
                    
                    logContributionMA = - sum(log(proposalsARMA(1).evaluateProposalMA(draw.maPacs,   muDenominator)));
                end;
                
            case draw.qs == oldState.qs
                if draw.qs > 0 %then also oldState.qs > 0
                    muEnumerator =  draw.maPacs;
                    muDenominator = oldState.maPacs;
                    
                    logContributionMA = sum(log(proposalsARMA(1).evaluateProposalMA(oldState.maPacs, muEnumerator))) -...
                                        sum(log(proposalsARMA(1).evaluateProposalMA(draw.maPacs,     muDenominator)));
                else %draw.qs == oldState.qs == 0
                    logContributionMA = 0;
                end;
                
            case draw.qs < oldState.qs %thus oldState.qs > 0 in any case
                if draw.qs > 0 
                    muEnumerator =  [draw.maPacs; zeros(oldState.qs - draw.qs,1)];
                    muDenominator = [oldState.maPacs(1:draw.qs)];   
                    
                    logContributionMA = sum(log(proposalsARMA(1).evaluateProposalMA(oldState.maPacs, muEnumerator))) -...
                                        sum(log(proposalsARMA(1).evaluateProposalMA(draw.maPacs,     muDenominator)));
                else 
                    muEnumerator = zeros(size(oldState.maPacs));
                    
                    logContributionMA = sum(log(proposalsARMA(1).evaluateProposalMA(oldState.maPacs, muEnumerator)));
                end;
        end;
        
        logContribution = logContributionMA + logContributionAR;
end