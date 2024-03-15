% function displayConsolidatedMonteCarloResultsSaveType2Single(settings, doPlots)
    cumsums = cumsum(ones(1,settings.draws));
    disp('Unconditional Posterior Means AR Coefficients')
    for i = 1:settings.pMax
        disp(['Order: ', num2str(i)]);
        if i <= size(arParametersSeries,1)
            temp=arParametersSeries(i,:);
            if doPlots
                figure;
                temp2 = temp;
                temp2(isnan(temp2))= 0;                
                plot(cumsum(temp2) ./ cumsums);
                title(['Empirical Average AR Parameters ' num2str(i)]);
            end;
            temp = temp(settings.burnIn:end);
            disp(mean(temp(isfinite(temp)==1)));
        else
            disp('NaN');
        end;
    end;

    disp('Unconditional Posterior Means MA Coefficients')
    for i = 1:settings.qMax
        disp(['Order: ', num2str(i)]);
        if i <= size(maParametersSeries,1)
            temp=maParametersSeries(i,:);
            if doPlots
                figure;
                temp2 = temp;
                temp2(isnan(temp2))= 0;                
                plot(cumsum(temp2) ./ cumsums);
                title(['Empirical Average MA Parameters ' num2str(i)]);
            end;
            temp = temp(settings.burnIn:end);
            disp(mean(temp(isfinite(temp)==1)));
        else
            disp('NaN');
        end;
    end;

    disp('Unconditional Posterior Medians AR Coefficients')
    for i = 1:settings.pMax
        disp(['Order: ', num2str(i)]);
        if i <= size(arParametersSeries,1)
            temp=arParametersSeries(i,settings.burnIn:end);
            disp(median(temp(isfinite(temp)==1)));
        else
            disp('NaN');
        end;
    end;
    
    disp('Unconditional Posterior Medians MA Coefficients')
    for i = 1:settings.qMax
        disp(['Order: ', num2str(i)]);
        if i <= size(maParametersSeries,1)
            temp=maParametersSeries(i,settings.burnIn:end);
            disp(median(temp(isfinite(temp)==1)));
        else
            disp('NaN');
        end;
    end;
    
    pqtest = [pSeries(settings.burnIn:end) qSeries(settings.burnIn:end)];
    
    x = 0:1:settings.pMax;
    z = 0:1:settings.qMax;

    bintest = cell(1);
    bintest{1} = x;
    bintest{2} = z;

    [nelements, centers] = hist3(pqtest,'Edges',bintest);
    [maxPQ, ind] = max(nelements(:));
    [m,n] = ind2sub(size(nelements),ind);

    pPostMax = m-1;
    qPostMax = n-1;

    pqSieve = ((pqtest(:,1) == pPostMax) & (pqtest(:,2) == qPostMax));
    
    disp('Conditional Means and Medians AR');
    for i = 1:settings.pMax
        disp(['Order: ', num2str(i)]);
        if i <= size(arParametersSeries,1)
            temp=arParametersSeries(i,settings.burnIn:end);
            if doPlots & (i <= pPostMax)
                figure;
                plot(transpose(cumsum(temp(pqSieve))) ./ cumsum(pqSieve(pqSieve == 1)));
                title(['Conditional Recursive Mean AR Parameter ' num2str(i)]);
            end;
            temp = temp(pqSieve);
            disp(['Mean: ' num2str(mean(temp))]);
            disp(['Median: ' num2str(median(temp))]);
        else
            disp('NaN');
        end;
    end;
    
    disp('Conditional Means and Medians MA');
    for i = 1:settings.qMax
        disp(['Order: ', num2str(i)]);
        if i <= size(maParametersSeries,1)
            temp=maParametersSeries(i,settings.burnIn:end);
            if doPlots & (i <= qPostMax)
                figure;
                plot(transpose(cumsum(temp(pqSieve))) ./ cumsum(pqSieve(pqSieve == 1)));
                title(['Conditional Recursive Mean MA Parameter ' num2str(i)]);
            end;
            temp = temp(pqSieve);
            disp(['Mean: ' num2str(mean(temp))]);
            disp(['Median: ' num2str(median(temp))]);
        else
            disp('NaN');
        end;
    end;
    
    disp('Conditional Mean and Median Sigma');
    temp = sigmaESeries(settings.burnIn:end);
    if doPlots
        figure;
        plot(cumsum(temp(pqSieve)) ./ cumsum(pqSieve(pqSieve == 1)));
        title(['Conditional Recursive Mean Sigma']);
    end;
    temp = temp(pqSieve);
    disp(['Mean: ' num2str(mean(temp))]);
    disp(['Median: ' num2str(median(temp))]);
    
    x = 0:1:settings.pMax;
    z = 0:1:settings.qMax;
    figure;
    hist(pSeries(settings.burnIn:end),x);
    figure;
    hist(qSeries(settings.burnIn:end),z);

    figure;
    bintest = cell(1);
    bintest{1} = x;
    bintest{2} = z;
    hist3(pqtest,'Edges',bintest);
    
%     a=results(dataset,chain).seriesLogPosterior(2:end)==results(dataset,chain).seriesLogPosterior(1:end-1);
%     b=1:length(a);
%     c=cumsum(a')./b;
%     figure; plot(c);

% end