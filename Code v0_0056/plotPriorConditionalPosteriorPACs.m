%plot Prior vs ConditionalPosterior
burnIn = 1000000;

pqtest = [pSeries(burnIn:end) qSeries(burnIn:end)];
x = 0:1:10;
z = 0:1:10;
bintest = cell(1);
bintest{1} = x;
bintest{2} = z;
[nelements, centers] = hist3(pqtest,'Edges',bintest);
[maxPQ, ind] = max(nelements(:));
[m,n] = ind2sub(size(nelements),ind);
pqSieve = ((pqtest(:,1) == m-1) & (pqtest(:,2) == n-1));

arPacsSeriesCropped = arPacsSeries(:,burnIn:end);
arPacsSeriesCropped = arPacsSeriesCropped(:,pqSieve);
for cntr = 1:m-1
    figure; hold on;
    temp = arPacsSeriesCropped(cntr,:);
%     temp = temp(pqSieve);
    [f, xi] = ksdensity(temp);
    plot(xi,f,'k','LineWidth',1.5);
    ezplot(@(x) truncatedNormalImproperPrior(x,0,0.5),[-1.1,1.1]);
    title(['AR PAC ' num2str(cntr)]);
    axis('auto');
    legend('Conditional Posterior','Prior');
end;

maPacsSeriesCropped = maPacsSeries(:,burnIn:end);
maPacsSeriesCropped = maPacsSeriesCropped(:,pqSieve);
for cntr = 1:n-1
    figure; hold on;
     temp = maPacsSeriesCropped(cntr,:);
%      temp = temp(pqSieve);
    [f, xi] = ksdensity(temp);
    plot(xi,f,'k','LineWidth',1.5);
    ezplot(@(x) truncatedNormalImproperPrior(x,0,0.5),[-1.1,1.1]);
    title(['MA PAC ' num2str(cntr)]);
    legend('Conditional Posterior','Prior');
end;
