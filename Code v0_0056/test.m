disp('Marginal Posteriors AR Coefficients');
figure;
for i = 1:p_max
    subplot(round(sqrt(p_max)),round(sqrt(p_max))+round(1-rem(sqrt(p_max),1)),i);
    hist(AR_Params_State(:,i),[-4:0.01:4]);
end;