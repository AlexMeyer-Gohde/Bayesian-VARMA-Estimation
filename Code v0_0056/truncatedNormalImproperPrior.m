function val = truncatedNormalImproperPrior(x, mu, sigma)
    val = zeros(size(x));
    for cntr = 1:max(size(x))
        if (x(cntr) > -1) && (x(cntr) < 1)
            val(cntr) = normpdf(x(cntr),mu,sigma);   
        else
            val(cntr) = 0;
        end;
    end;
end