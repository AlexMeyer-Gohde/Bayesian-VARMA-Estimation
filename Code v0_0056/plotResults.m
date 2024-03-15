function plotResults(states, settings, accepted, y)

%Plot results
disp('Unconditional Posterior Means AR Coefficients')
for i = 1:settings.pMax
    disp(['Order: ', num2str(i)]);
    temp=padcat(states(settings.burnIn:end).arParameters);
    if i <= size(temp,1)
        temp=temp(i,:);
        disp(mean(temp(isfinite(temp)==1)));
    else
        disp('NaN');
    end;
end;

disp('Unconditional Posterior Means MA Coefficients')
for i = 1:settings.qMax
    disp(['Order: ', num2str(i)]);
    temp=padcat(states(settings.burnIn:end).maParameters);
    if i <= size(temp,1)
        temp=temp(i,:);
        disp(mean(temp(isfinite(temp)==1)));
    else
        disp('NaN');
    end;
end;

disp('Unconditional Posterior Medians AR Coefficients')
for i = 1:settings.pMax
    disp(['Order: ', num2str(i)]);
    temp=padcat(states(settings.burnIn:end).arParameters);
    if i <= size(temp,1)
        temp=temp(i,:);
        disp(median(temp(isfinite(temp)==1)));
    else
        disp('NaN');
    end;
end;

temp = cat(1,states(settings.burnIn:end).sigmaEs);	
disp(['Unconditional Posterior Mean \sigma_e: ' num2str(mean(temp)) ]);

disp(['Unconditional Posterior Median \sigma_e: ' num2str(median(temp)) ]);


temp = cat(1,states(settings.burnIn:end).estimatedVariance);	
disp(['Standard deviation data vector: ' num2str(sqrt(var(y)))]);
disp(['Unconditional Posterior Mean estimated Standard Deviation: ' num2str(mean(sqrt(temp(settings.burnIn:end))) )]);

disp(['Unconditional Posterior Median estimated Standard Deviation: ' num2str(sqrt(median(temp(settings.burnIn:end))) )]);


p_count_raw = cat(1,states(settings.burnIn:end).ps);
q_count_raw = cat(1,states(settings.burnIn:end).qs);
AcceptanceRate = accepted/settings.draws
p_count = p_count_raw .* accepted;
x = 0:1:settings.pMax;
z = 0:1:settings.qMax;
figure;
hist(p_count_raw,x);
figure;
hist(q_count_raw,z);

figure;
pqtest = [p_count_raw q_count_raw];
bintest = cell(1);
bintest{1} = x;
bintest{2} = z;
hist3(pqtest,'Edges',bintest);

temp=padcat(states.logPosterior);
a=temp(2:end)==temp(1:end-1);
b=1:length(a);
c=cumsum(a')./b;
figure; plot(c);

end