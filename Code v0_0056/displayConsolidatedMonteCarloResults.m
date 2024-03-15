cntrCorrect = 0;
disp('Without Burn in');
for cntrRecData = 21:100
    for cntrRecChains = 1:1
        disp(['Data set Nr.: ' num2str(cntrRecData) '; Chain Nr.: ' num2str(cntrRecChains)]);
        disp(['Posterior Mode: p=' num2str(consolidatedResults(cntrRecData, cntrRecChains).pPostMax) ' q=' num2str(consolidatedResults(cntrRecData, cntrRecChains).qPostMax) ]);
        if consolidatedResults(cntrRecData, cntrRecChains).pPostMax == 3 && consolidatedResults(cntrRecData, cntrRecChains).qPostMax == 2
            cntrCorrect = cntrCorrect + 1;
        end;
    end;
end;

disp(cntrCorrect);

% x = 0:1:settings.pMax;
% z = 0:1:settings.qMax;
% 
% bintest = cell(1);
% bintest{1} = x;
% bintest{2} = z;
% 
% disp('With Burn in');
% for cntrRecData = 21:100
%     for cntrRecChains = 1:1       
%         pRawBurnin = consolidatedResults(cntrRecData, cntrRecChains).pRaw(settings.burnIn:end);
%         qRawBurnin = consolidatedResults(cntrRecData, cntrRecChains).qRaw(settings.burnIn:end);
%         
%         disp(['Data set Nr.: ' num2str(cntrRecData) '; Chain Nr.: ' num2str(cntrRecChains)]);
%         
%         [nelements, centers] = hist3([pRawBurnin qRawBurnin],'Edges',bintest);
%         [maxPQ, ind] = max(nelements(:));
%         [m,n] = ind2sub(size(nelements),ind);        
%         
%         disp(['Posterior Mode: p=' num2str(m-1) ' q=' num2str(n-1) ]);
%     end;
% end;
