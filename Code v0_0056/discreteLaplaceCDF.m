function y = discreteLaplaceCDF(k_max,lambda)
    %Returns a matrix of values from the CDF of the discretized Laplace
    %distribution following ehlers brooks (2004) and Troughton Goodsill
     
    y = zeros(k_max+1,k_max+1);
    vecKPrime = transpose(0:1:k_max);
    summing = tril(ones(k_max+1,k_max+1));
    
    for cntrK = 1:k_max+1
        vecK = ones(k_max+1,1)*(cntrK-1);
        pdfVal = abs(vecKPrime-vecK); %abs(k'-k)
        pdfVal = -lambda * pdfVal; %-lambda*(abs(k'-k))
        pdfVal = exp(pdfVal);%exp(-lambda*(abs(k'-k)))
        constant = sum(pdfVal);
        
        y(:,cntrK) = (summing*pdfVal)/constant;
    end;    
end