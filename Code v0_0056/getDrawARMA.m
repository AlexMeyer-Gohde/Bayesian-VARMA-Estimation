function draw = getDrawARMA(oldState, proposalsARMA, settings)
    %ONLY ONE PROCESS ATM!!!! 
    %Initialize struct
    draw = getEmptyDrawStruct();
    
    %Draw new p's and parameters
    if settings.useSamePriorProposal == 1
        for cntr = 1:settings.processCount
            draw.ps(cntr) = proposalsARMA(1).proposalP(oldState.ps(cntr));            
            draw.qs(cntr) = proposalsARMA(1).proposalQ(oldState.qs(cntr));
            
%             draw.logProposal(cntr) = proposalsARMA(1).evaluateProposalP([oldState.ps(cntr) draw.ps(cntr)]) + ...
%                                      proposalsARMA(1).evaluateProposalQ([oldState.qs(cntr) draw.qs(cntr)]);

            %Draw new PACs and compute Parameters
            if draw.ps(cntr) > 0
                if draw.ps(cntr) > oldState.ps(cntr)
                    if oldState.ps > 0
                        draw.arPacs = proposalsARMA(1).proposalAR([oldState.arPacs(:,cntr); zeros(draw.ps(cntr) - oldState.ps(cntr),1)]);
                    else
                        draw.arPacs = proposalsARMA(1).proposalAR(zeros(draw.ps(cntr) - oldState.ps(cntr),1));
                    end;
                    draw.arParameters(:,cntr) = getARParametersFromPACs(draw.arPacs(:,cntr),draw.ps(cntr));
                else
                    draw.arPacs = [draw.arPacs proposalsARMA(1).proposalAR(oldState.arPacs(1:draw.ps(cntr),cntr))];
                    draw.arParameters(:,cntr) = getARParametersFromPACs(draw.arPacs(:,cntr),draw.ps(cntr));
                end;
            else
                draw.arPacs = [];
                draw.arParameters = [];
            end;
            
            if draw.qs(cntr) > 0
                if draw.qs(cntr) > oldState.qs(cntr)
                    if oldState.qs > 0
                        draw.maPacs = proposalsARMA(1).proposalMA([oldState.maPacs(:,cntr); zeros(draw.qs(cntr) - oldState.qs(cntr),1)]);
                    else
                        draw.maPacs =  proposalsARMA(1).proposalMA(zeros(draw.qs(cntr) - oldState.qs(cntr),1));
                    end;
                    draw.maParameters(:,cntr) = -getARParametersFromPACs(draw.maPacs(:,cntr),draw.qs(cntr));
                else
                    draw.maPacs = proposalsARMA(1).proposalMA(oldState.maPacs(1:draw.qs(cntr),cntr));
                    draw.maParameters(:,cntr) = -getARParametersFromPACs(draw.maPacs(:,cntr),draw.qs(cntr));
                end;
            else
                draw.maPacs = [];
                draw.maParameters = [];
            end;

            %Draw new sigmaE
            draw.sigmaEs = [draw.sigmaEs proposalsARMA.proposalSigmaE(oldState.sigmaEs(cntr))];
        end;
    else
        %NOT YET UPDATED FOR MA!!!!!!!!!!!!!!!!!!!!!
%         for cntr = 1:settings.processCount    
%             draw.ps(cntr) = proposalsARMA(cntr).proposalP(); 
% 
%             %Draw new AR Parameters (to be replaced by PACs)
%             if draw.ps(cntr) > oldState.ps(cntr)
%                 draw.arParameters = [draw.arParameters proposalsARMA(cntr).proposalAR([oldState.arParameters(:,cntr); zeros(draw.ps(cntr) - oldState.ps(cntr),1)])];
%             else
%                 draw.arParameters = [draw.arParameters proposalsARMA(cntr).proposalAR(oldState.arParameters(1:draw.ps(cntr),cntr))];
%             end;
% 
%             %Draw new sigmaE
%             draw.sigmaEs = [draw.sigmaEs proposalsARMA(cntr).proposalSigmaE(oldState.sigmaEs(cntr))];
%         end;
    end;     
end