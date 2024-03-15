% This script plots stuff. p_mode needs to be set in advance!

%DN -> DN: Check this!!!!! The treatment of the BurnIn might not be
%correct!!!!!!!!!!!!!!!!!


AR_Params_TMP = AR_Params_State(p_count_raw==p_mode,:);
for i=1:p_mode
    [f,xi] = ksdensity(AR_Params_TMP(BurnIn:end,i));
    h=figure('Name', ['AR Coefficient ' num2str(i) ' (conditional on p=' num2str(p_mode) ')']);
%     title(['AR Coefficient ' num2str(i) ' prior and posterior (conditional on p=' num2str(p_mode) ')']);
    hold on;
    fplot(prior_ar_params, [-3 3 ]);
    plot(xi,f,'r');
    hold off;
end;



disp(['Posterior Means AR Coefficients conditional on p=' num2str(p_mode) ]);
for i = 1:p_mode
    disp(['Order: ', num2str(i)]);
    temp=AR_Params_TMP(BurnIn:end,i);
    disp(mean(temp(isfinite(temp)==1)));
end;

disp(['Posterior Medians AR Coefficients conditional on p=' num2str(p_mode) ]);
for i = 1:p_mode
    disp(['Order: ', num2str(i)]);
    temp=AR_Params_TMP(BurnIn:end,i);
    disp(median(temp(isfinite(temp)==1)));
end;


disp(['90 percent credible sets']);
for i = 1:p_mode
    disp(['Order: ', num2str(i)]);
    temp=AR_Params_TMP(BurnIn:end,i);
    disp(prctile((temp(isfinite(temp)==1)),[5 95]));
end;


sigmaTMP = cell2mat(draws(:,3));
sigmaTMP = sigmaTMP(p_count_raw==p_mode);
disp(['\sigma_e']);
disp(prctile(sigmaTMP,[5 95]));

disp(['Conditional Posterior Mean \sigma_e: ' num2str(mean(sigmaTMP(BurnIn:end))) ]);

disp(['Conditional Posterior Median \sigma_e: ' num2str(median(sigmaTMP(BurnIn:end))) ]);


figure('Name', 'Prior and Posterior \sigma_e');
hold on;
% title('Prior and Posterior \sigma_e');
[f,xi] = ksdensity(cell2mat(draws(:,3)));
plot(xi,f,'r');
fplot(prior_sigma_e, [0 max(cell2mat(draws(:,3)))]);
hold off;