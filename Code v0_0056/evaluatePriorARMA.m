function [priorValue, dropDraw] = evaluatePriorARMA(state, priorsARMA, settings)
    %Use same prior for all processes? Uses Prior on Pacs!!! 
    if settings.useSamePriorProposal == 1
        dropDraw = 0;
        for cntrProcess = 1:settings.processCount
            if ~dropDraw
                priorValue = 0;
                if state.ps > 0
                    priorAR = priorsARMA(1).priorAR(state.arPacs(:,cntrProcess));
                    priorValue = sum(log(priorAR));
                else
                    priorAR = 1/0;
                end;
                if state.qs > 0
                    priorMA = priorsARMA(1).priorMA(state.maPacs(:,cntrProcess));
                    priorValue = priorValue + sum(log(priorMA));
                else
                    priorMA = 1/0;
                end;             
                if (min(priorAR) > 0) && (min(priorMA) > 0) && (priorsARMA(1).priorSigmaE(state.sigmaEs(cntrProcess)) > 0)
                    priorValue = priorValue + ...
                            log(priorsARMA(1).priorP(state.ps(cntrProcess)));
                    priorValue = priorValue + ...
                            log(priorsARMA(1).priorQ(state.qs(cntrProcess)));
                    priorValue = priorValue +  ...
                            log(priorsARMA(1).priorSigmaE(state.sigmaEs(cntrProcess)));
                else
                    dropDraw = 1;
                end;
            end;
        end;
    else %DOES NOT WORK FOR MULTIPLE PROCESSES ATM
        for cntrProcess = 1:settings.processCount
            priorValue = priorValue + ...
                        log(priorsARMA(cntrProcess).priorP(state.ps(cntrProcess)))  +  ...
                        log(priorsARMA(cntrProcess).priorSigmaE(state.sigmaEs(cntrProcess)));
            priorValue = priorValue + ...
                log(priorsARMA(cntrProcess).priorSigmaE(state.sigmaEs(cntrProcess)))+ ...
                log(priorsARMA(cntrProcess).priorP(state.ps(cntrProcess)));
        end;        
    end;
end