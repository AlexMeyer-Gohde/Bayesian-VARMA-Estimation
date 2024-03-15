function plotConditionalPQ(p,q,process,states,settings)
    figure;
    disp('Unconditional Posterior Means AR Coefficients')
    temp1 = states(settings.burnIn:end);
    disp(['Conditional on p=', num2str(p) , ', q=', num2str(q)]);
    temp=padcat(temp1.arParameters((temp1.ps(process) == p) & (temp1.qs(process) == q)));
    if p <= size(temp,1)
        temp=temp(p,settings.burnIn:end);
        disp(mean(temp(isfinite(temp)==1)));
    else
        disp('NaN');
    end;
end