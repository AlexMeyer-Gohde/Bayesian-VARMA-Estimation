% function [IRFs]=posterior_IRF(IRF_horizon,number_of_draws,credible_set_percentile,burnin,data_set_name,skip)
noninvertible=0;
do.IRF=1;
do.COR=0;
do.HPLAMBDA = 0; %set to zero to switch off HP filtering
horizon.IRF=40;
horizon.COR=10;
skip = 30; %Must divide number_of_draws less burnin without remainder
burnin = 1000000;
number_of_draws = 4000000;
credible_set_percentile = 0.1;
data_set_name = 'X:\US GDP Shenanigans\One-Sided HP Filtered Likelihood Widened Proposals\One Sided Filtered 5 5\Results Bayesian_Monte Carlo Synthetic Philippe (-Inf,-Inf)_Sample Size 100_Draws 4000000_Date 28-Apr-2014_Data Series 1_Chain 1 Save Type 2.mat';

%correct number of draws for burnin
number_of_draws = number_of_draws - burnin;

[modelinfo, parameters] = exampl1_extractGenerateARMAShockData();
IRF=zeros((modelinfo.num_exog+modelinfo.num_endo)*modelinfo.num_exog,horizon.IRF,number_of_draws/skip);
COV=zeros((modelinfo.num_endo)*modelinfo.num_exog,1+(1+2*horizon.COR)*modelinfo.num_endo,number_of_draws/skip);
%extract draws: Daniel?

load(data_set_name);

for draw=1:skip:number_of_draws
    state.ps=pSeries(draw+burnin);
    state.qs=qSeries(draw+burnin);
    state.arParameters=arParametersSeries(:,draw+burnin);
    state.maParameters=maParametersSeries(:,draw+burnin);
    state.sigmaEs=sigmaESeries(draw+burnin);

    if do.COR && do.IRF
        [IRF(:,:,(draw-1)/skip +1),COV(:,:,(draw-1)/skip +1)] = DSGE_IRF_COR_Extract(modelinfo,state,horizon,parameters,do,noninvertible);
    elseif do.IRF
        [IRF(:,:,(draw-1)/skip +1),~] = DSGE_IRF_COR_Extract(modelinfo,state,horizon,parameters,do,noninvertible);
    else
        [~,COV(:,:,(draw-1)/skip +1)] = DSGE_IRF_COR_Extract(modelinfo,state,horizon,parameters,do,noninvertible);
    end
end
if do.COR && do.IRF
    if noninvertible
        save(['IRF_COV_noninv_' date '.mat'],'IRF','COV'); 
    else
        save(['IRF_COV_' date '.mat'],'IRF','COV');
    end
elseif do.IRF
    if noninvertible
        save(['IRF_noninv_' date '.mat'],'IRF','COV');
    else
        save(['IRF_' date '.mat'],'IRF','COV');
    end
else
    save(['COV_' date '.mat'],'COV');
end



