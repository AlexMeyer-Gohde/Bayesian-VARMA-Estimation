credible_set_percentile=0.1;
skip = 30; %Must divide number_of_draws less burnin without remainder
burnin = 1000000;
number_of_draws = 4000000;
number_of_draws = number_of_draws - burnin;

    temp=COV;%sort(COV,3);
    COVs.mean_COV=mean(temp,3);
    COVs.median_COV=median(temp,3);
    COVs.lower_credible_COV=prctile(temp,100*credible_set_percentile,3);%temp(:,:,round(credible_set_percentile*number_of_draws/skip+1));
    COVs.upper_credible_COV=prctile(temp,100*(1-credible_set_percentile),3);%temp(:,:,round((1-credible_set_percentile)*number_of_draws/skip));



VARNAMES = ['capital    ',
            'consumption',
            'output     ',
            'labor      ',
            'interest   ',
            'investment ',
            'technology '];
       %[modelinfo, parameters] = exampl1_extractGenerateARMAShockData();
        CORRELATION_HORIZON=10;
[obs variables]=size(modelinfo.data);
obs_location=[3];
ENDOGENOUS_VARIABLE_NUMBER=6;
parameters.Z.p=3; %ar order here
parameters.Z.q=0; %ma order here
parameters.Z.P=[1.1689 -0.0732 -0.1224];
parameters.Z.Q=[];
parameters.Z.Sigma=0.5873^2;
[modelinfo,parameters]=solve_recursive_dsge(modelinfo,parameters); COR_posterior_mean = DSGE_COR(modelinfo, parameters,CORRELATION_HORIZON,1600);
parameters.Z.p=1; %ar order here
parameters.Z.q=0; %ma order here
parameters.Z.P=[0.95]; %ar parameters here
parameters.Z.Q=[]; %ma parameters here
parameters.Z.Sigma=0.712^2;
[modelinfo,parameters]=solve_recursive_dsge(modelinfo,parameters);COR_Hansen = DSGE_COR(modelinfo, parameters,CORRELATION_HORIZON,1600);
        


for j=0:CORRELATION_HORIZON
    Gamma_1_data{j+1}=modelinfo.data(1:obs-j,:)'*modelinfo.data(1+j:obs,:);
end
for j=1:2*length(Gamma_1_data)-1
    for i=1:variables
        if j-length(Gamma_1_data)<0
            Xcorr_data{i}(:,j)=Gamma_1_data{abs(j-(length(Gamma_1_data)+1))}(:,i)./((Gamma_1_data{1}(i,i).^(1/2))*(diag(Gamma_1_data{1}).^(1/2)));
            %Xcorr_data{i}(:,j)=(obs/(obs-abs(j-length(Gamma_1_data))))*Gamma_1_data{abs(j-(length(Gamma_1_data)+1))}(:,i)./((Gamma_1_data{1}(i,i).^(1/2))*(diag(Gamma_1_data{1}).^(1/2)));
        else
            Xcorr_data{i}(:,j)=Gamma_1_data{abs(j-(length(Gamma_1_data))+1)}(i,:)'./((Gamma_1_data{1}(i,i).^(1/2))*(diag(Gamma_1_data{1}).^(1/2)));
            %Xcorr_data{i}(:,j)=(obs/(obs-abs(j-length(Gamma_1_data))))*Gamma_1_data{abs(j-(length(Gamma_1_data))+1)}(i,:)'./((Gamma_1_data{1}(i,i).^(1/2))*(diag(Gamma_1_data{1}).^(1/2)));
        end
    end
end
for jjj=1:length(obs_location)
    for iii=1:length(obs_location)
        data_correlations{iii,jjj}=Xcorr_data{jjj}(iii,1:2*CORRELATION_HORIZON+1);
    end
end

for iii=1:ENDOGENOUS_VARIABLE_NUMBER
    for jjj=1:ENDOGENOUS_VARIABLE_NUMBER
        figure;
        if iii==jjj
        c_select=1+[CORRELATION_HORIZON:2*CORRELATION_HORIZON]*ENDOGENOUS_VARIABLE_NUMBER+jjj;
        xaxis=0:CORRELATION_HORIZON;
        if ismember(iii,obs_location)*ismember(jjj,obs_location);
            plot(xaxis, data_correlations{obs_location==iii, obs_location==jjj}(CORRELATION_HORIZON+1:end),'g',xaxis, COR_Hansen(iii,c_select),'r',xaxis,COR_posterior_mean(iii,c_select),'k',xaxis,COVs.median_COV(iii,c_select),'b',xaxis,COVs.lower_credible_COV(iii,c_select),':b',xaxis,COVs.upper_credible_COV(iii,c_select),':b')
            legend('Data','Hansen','Posterior Mode Model','Posterior Mode','Posterior 10% Bound','Posterior 90% Bound','Location','Best');
        else
            plot(xaxis, COR_Hansen(iii,c_select),'r',xaxis,COR_posterior_mean(iii,c_select),'k',xaxis,COVs.median_COV(iii,c_select),'b',xaxis,COVs.lower_credible_COV(iii,c_select),':b',xaxis,COVs.upper_credible_COV(iii,c_select),':b')
            legend('Hansen','Posterior Mode Model','Posterior Mode','Posterior 10% Bound','Posterior 90% Bound','Location','Best');
        end
        hold on; plot(xaxis,zeros(1,length(xaxis)),':k'); hold off
        eval(sprintf('title(''Autocorrelations of %s'')',VARNAMES(jjj,:)))
        ylabel('Correlation Coefficient');
        xlabel('j');
        axis([0,CORRELATION_HORIZON,-1,1])
        else
            c_select=1+[0:2*CORRELATION_HORIZON]*ENDOGENOUS_VARIABLE_NUMBER+jjj;
        xaxis=-CORRELATION_HORIZON:CORRELATION_HORIZON;
        if ismember(iii,obs_location)*ismember(jjj,obs_location);
            plot(xaxis, data_correlations{ obs_location==iii, obs_location==jjj},'g',xaxis, COR_Hansen(iii,c_select),'r',xaxis,COR_posterior_mean(iii,c_select),'k',xaxis,COVs.median_COV(iii,c_select),'b',xaxis,COVs.lower_credible_COV(iii,c_select),':b',xaxis,COVs.upper_credible_COV(iii,c_select),':b')
            legend('Data','Hansen','Posterior Mode Model','Posterior Mode','Posterior 10% Bound','Posterior 0% Bound','Location','Best');
        else
            plot(xaxis, COR_Hansen(iii,c_select),'r',xaxis,COR_posterior_mean(iii,c_select),'k',xaxis,COVs.median_COV(iii,c_select),'b',xaxis,COVs.lower_credible_COV(iii,c_select),':b',xaxis,COVs.upper_credible_COV(iii,c_select),':b')
            legend('Hansen','Posterior Mode Model','Posterior Mode','Posterior 10% Bound','Posterior 90% Bound','Location','Best');
        end
         hold on; plot(xaxis,zeros(1,length(xaxis)),':k'); hold off
        eval(sprintf('title(''Cross-Correlations of %s at t+j with %s at t'')',VARNAMES(iii,:),VARNAMES(jjj,:)))
        ylabel('Correlation Coefficient');
        xlabel('j');
        axis([-CORRELATION_HORIZON,CORRELATION_HORIZON,-1,1])
        end

    end
end