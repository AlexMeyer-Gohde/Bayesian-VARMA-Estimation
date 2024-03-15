%plot Prior vs ConditionalPosterior
burnIn = settings.burnIn;
arOrderMax = 3;
maOrderMax = 0;

% pqtest = [pSeries(burnIn:end) qSeries(burnIn:end)];
% x = 0:1:10;
% z = 0:1:10;
% bintest = cell(1);
% bintest{1} = x;
% bintest{2} = z;
% [nelements, centers] = hist3(pqtest,'Edges',bintest);
% [maxPQ, ind] = max(nelements(:));
% [m,n] = ind2sub(size(nelements),ind);

arPacsSeriesCropped = arPacsSeries(:,burnIn:end);
for cntr = 1:arOrderMax
    pqSieve = (pSeries(burnIn:end) >= cntr);
    
    figure; hold on;
    temp = arPacsSeriesCropped(cntr,pqSieve);
%     temp = temp(pqSieve);
    [f, xi] = ksdensity(temp);
    plot(xi,f,'k','LineWidth',1.5);
    ezplot(@(x) truncatedNormalImproperPrior(x,0,0.5),[-1.1,1.1]);
    title(['AR PAC ' num2str(cntr)]);
    axis('auto');
    legend('Posterior','Prior');
end;

maPacsSeriesCropped = maPacsSeries(:,burnIn:end);
for cntr = 1:maOrderMax
    pqSieve = (qSeries(burnIn:end) >= cntr);
    
    figure; hold on;
     temp = maPacsSeriesCropped(cntr,pqSieve);
%      temp = temp(pqSieve);
    [f, xi] = ksdensity(temp);
    plot(xi,f,'k','LineWidth',1.5);
    ezplot(@(x) truncatedNormalImproperPrior(x,0,0.5),[-1.1,1.1]);
    title(['MA PAC ' num2str(cntr)]);
    legend('Posterior','Prior');
end;
