function displayResultsSaveType2(results, settings, doPlots)

    cumsums = cumsum(ones(1,settings.draws));
    disp('Unconditional Posterior Means AR Coefficients')
    for i = 1:settings.pMax
        disp(['Order: ', num2str(i)]);
        if i <= size(results.arParametersSeries,1)
            temp=results.arParametersSeries(i,:);
            if doPlots
                figure;
                temp2 = temp;
                temp2(isnan(temp2))= 0;                
                plot(cumsum(temp2) ./ cumsums);
                title(['Empirical Average AR Parameters ' num2str(i)]);
            end;
            temp = temp(settings.burnIn+1:end);
            disp(mean(temp(isfinite(temp)==1)));
        else
            disp('NaN');
        end;
    end;
    
    disp('Unconditional Posterior Medians AR Coefficients')
    for i = 1:settings.pMax
        disp(['Order: ', num2str(i)]);
        if i <= size(results.arParametersSeries,1)
            temp=results.arParametersSeries(i,settings.burnIn+1:end);
            disp(median(temp(isfinite(temp)==1)));
        else
            disp('NaN');
        end;
    end;
    
    disp('Unconditional Posterior Means MA Coefficients')
    for i = 1:settings.qMax
        disp(['Order: ', num2str(i)]);
        if i <= size(results.maParametersSeries,1)
            temp=results.maParametersSeries(i,:);
            if doPlots
                figure;
                temp2 = temp;
                temp2(isnan(temp2))= 0;                
                plot(cumsum(temp2) ./ cumsums);
                title(['Empirical Average MA Parameters ' num2str(i)]);
            end;
            temp = temp(settings.burnIn+1:end);
            disp(mean(temp(isfinite(temp)==1)));
        else
            disp('NaN');
        end;
    end; 
    
    disp('Unconditional Posterior Medians MA Coefficients')
    for i = 1:settings.qMax
        disp(['Order: ', num2str(i)]);
        if i <= size(results.maParametersSeries,1)
            temp=results.maParametersSeries(i,settings.burnIn+1:end);
            disp(median(temp(isfinite(temp)==1)));
        else
            disp('NaN');
        end;
    end;
    
    
    disp('Unconditional Posterior Mean Sigma_e');
    temp = results.sigmaESeries(settings.burnIn+1:end);
    disp(mean(temp));    
    
    disp('Unconditional Posterior Median Sigma_e');
    disp(median(temp));  
    
    if doPlots
        figure;
        plot(transpose(cumsum(results.sigmaESeries)) ./ cumsums);
        title(['Empirical Average Sigma_e ']);
    end;


    pqtest = [results.pSeries(settings.burnIn+1:end) results.qSeries(settings.burnIn+1:end)];    
    x = 0:1:settings.pMax;
    z = 0:1:settings.qMax;
    bintest = cell(1);
    bintest{1} = x;
    bintest{2} = z;
    [nelements, centers] = hist3(pqtest,'Edges',bintest);
    [maxPQ, ind] = max(nelements(:));
    [m,n] = ind2sub(size(nelements),ind);

    results.pPostMax = m-1;
    results.qPostMax = n-1;
    
    pqSieve = ((pqtest(:,1) == results.pPostMax) & (pqtest(:,2) == results.qPostMax));
    
    disp('Conditional Means and Medians AR');
    for i = 1:settings.pMax
        disp(['Order: ', num2str(i)]);
        if i <= size(results.arParametersSeries,1)
             temp=results.arParametersSeries(i,:);
            if doPlots & (i <= results.pPostMax)
                figure;
                plot(transpose(cumsum(temp(pqSieve))) ./ cumsum(pqSieve(pqSieve == 1)));
                title(['Empirical Average AR Parameters (Conditional)' num2str(i)]);
            end;
            temp=results.arParametersSeries(i,settings.burnIn+1:end);
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
        if i <= size(results.maParametersSeries,1)
            temp=results.maParametersSeries(i,:);
            if doPlots && (i <= results.qPostMax)
                figure;
                plot(transpose(cumsum(temp(pqSieve))) ./ cumsum(pqSieve(pqSieve == 1)));
                title(['Empirical Average MA Parameters (Conditional)' num2str(i)]);
            end;
	    temp=results.maParametersSeries(i,settings.burnIn+1:end);
            temp = temp(pqSieve);
            disp(['Mean: ' num2str(mean(temp))]);
            disp(['Median: ' num2str(median(temp))]);
        else
            disp('NaN');
        end;
    end;
    
    disp('Conditional Posterior Mean Sigma_e');
    temp = results.sigmaESeries(settings.burnIn+1:end);
    temp = temp(pqSieve);
    disp(mean(temp));    
    
    disp('Conditional Posterior Median Sigma_e');
    disp(median(temp));  
    
    if doPlots
        figure;
        plot(cumsum(temp) ./ cumsum(pqSieve(pqSieve == 1)));
        title(['Empirical Average Sigma_e (Conditional) ']);
    end;
    
    figure;
    hist(results.pSeries(settings.burnIn+1:end),x);
    figure;
    hist(results.qSeries(settings.burnIn+1:end),z);

    figure;
    hist3(pqtest,'Edges',bintest);
    
    a=results.sigmaESeries(2:end)==results.sigmaESeries(1:end-1);
    b=1:length(a);
    c=cumsum(a')./b;
    figure; plot(c);

end