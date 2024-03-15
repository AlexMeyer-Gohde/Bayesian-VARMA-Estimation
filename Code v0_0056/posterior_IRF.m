% function [IRFs]=posterior_IRF(IRF_horizon,number_of_draws,credible_set_percentile,burnin,data_set_name,skip)
noninvertible=1
[modelinfo, parameters] = exampl1_extractGenerateARMAShockData();
IRF=zeros((modelinfo.num_exog+modelinfo.num_endo)*modelinfo.num_exog,IRF_horizon,number_of_draws/skip);
%extract draws: Daniel?

load(data_set_name);

for draw=1:skip:number_of_draws
    state.ps=pSeries(draw+burnin);
    state.qs=qSeries(draw+burnin);
    state.arParameters=arParametersSeries(:,draw+burnin);
    state.maParameters=maParametersSeries(:,draw+burnin);
    state.sigmaEs=sigmaESeries(draw+burnin);
    % if min(isfinite(state.arParameters)) == 0 || min(isfinite(state.maParameters)) == 0 || isfinite(state.sigmaEs) == 0
    %     disp('something fishy just happened');
    % end;
    IRF(:,:,(draw-1)/skip +1) = DSGE_IRF_Extract(modelinfo,state,IRF_horizon,parameters,noninvertible);
    % if min(IRF(isfinite)) == 0
    %     disp('something fishy just happened');
    % end;
end


temp=sort(IRF,3);
IRFs.mean_IRF=mean(temp,3);
IRFs.median_IRF=median(temp,3);
IRFs.lower_credible_IRF=temp(:,:,round(credible_set_percentile*number_of_draws/skip+1));
IRFs.upper_credible_IRF=temp(:,:,round((1-credible_set_percentile)*number_of_draws/skip));
IRFs.IRF=IRF;

