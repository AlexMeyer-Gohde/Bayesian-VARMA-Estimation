credible_set_percentile=0.1;
skip = 300; %Must divide number_of_draws less burnin without remainder
burnin = 1000000;
number_of_draws = 4000000;
number_of_draws = number_of_draws - burnin;
    temp=IRF; %sort(IRF,3);
    IRFs.mean_IRF=mean(temp,3);
    IRFs.median_IRF=median(temp,3);
    IRFs.lower_credible_IRF=prctile(temp,100*credible_set_percentile,3)%temp(:,:,round(credible_set_percentile*number_of_draws/skip+1));
    IRFs.upper_credible_IRF=prctile(temp,100*(1-credible_set_percentile),3)%temp(:,:,round((1-credible_set_percentile)*number_of_draws/skip));
    IRFs.IRF=IRF;





VARNAMES = ['capital    ',
            'consumption',
            'output     ',
            'labor      ',
            'interest   ',
            'investment ',
            'technology '];
[modelinfo, parameters] = exampl1_extractGenerateARMAShockData();
parameters.Z.p=3; %ar order here
parameters.Z.q=0; %ma order here
parameters.Z.P=[1.1681 -0.0725 -0.1215];
parameters.Z.Q=[];
parameters.Z.Sigma=0.5733^2;
[modelinfo,parameters]=solve_recursive_dsge(modelinfo,parameters); IRF_posterior_mean = DSGE_IRF(modelinfo, parameters,40);
parameters.Z.p=1; %ar order here
parameters.Z.q=0; %ma order here
parameters.Z.P=[0.95]; %ar parameters here
parameters.Z.Q=[]; %ma parameters here
parameters.Z.Sigma=0.712^2;
[modelinfo,parameters]=solve_recursive_dsge(modelinfo,parameters);IRF_Hansen = DSGE_IRF(modelinfo, parameters,40);
        
        [variables horizon]=size(IRF_Hansen)
        for i=1:variables
            figure; hold on;
            plot(0:horizon-1,IRF_Hansen(i,:),'r',0:horizon-1,IRF_posterior_mean(i,:),'b',0:horizon-1,IRFs.median_IRF(i,:),'g',0:horizon-1,IRFs.lower_credible_IRF(i,:),':g',0:horizon-1,IRFs.upper_credible_IRF(i,:),':g')
            title(VARNAMES(i,:));
            legend('Hansen','Posterior Mode Model','Posterior Mode IRF','Posterior IRF 10% Bound','Posterior IRF 90% Bound','Location','Best');
            plot(0:horizon-1,zeros(1,horizon),'k')
        end