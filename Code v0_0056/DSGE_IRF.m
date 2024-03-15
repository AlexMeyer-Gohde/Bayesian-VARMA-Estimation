function IRF = DSGE_IRF(modelinfo, parameters,IRF_horizon)
max_p_q=max(parameters.Z.p,parameters.Z.q);
z = zeros(modelinfo.num_exog,IRF_horizon+max_p_q);
y = zeros(modelinfo.num_endo,IRF_horizon+1);
IRF=zeros((modelinfo.num_exog+modelinfo.num_endo)*modelinfo.num_exog,IRF_horizon);
if modelinfo.num_exog>1
	Std_Dev=chol(parameters.Z.Sigma);
else
	Std_Dev=(parameters.Z.Sigma)^(1/2);
end
for exog_select=1:modelinfo.num_exog
    epsilon = zeros(modelinfo.num_exog,IRF_horizon+max_p_q);
    epsilon(:,end-IRF_horizon+1)=Std_Dev(:,exog_select);%epsilon(exog_select,end-IRF_horizon+1)=1;
    for time_step =1:IRF_horizon
        z(:,time_step+max_p_q+1) = parameters.Z.P*vec(z(:,time_step+max_p_q:-1:time_step+max_p_q-parameters.Z.p+1)) +...
            [eye(modelinfo.num_exog), parameters.Z.Q]*vec(epsilon(:,time_step+max_p_q:-1:time_step+max_p_q-parameters.Z.q));
        if parameters.Z.q>0
            y(:,time_step+1)=parameters.X.Lambda*y(:,time_step)+parameters.X.Phi*vec(z(:,time_step+max_p_q+1:-1:time_step+max_p_q-parameters.Z.p+2))...
                +parameters.X.Theta*vec(epsilon(:,time_step+max_p_q:-1:time_step+max_p_q-parameters.Z.q+1));
        else
            y(:,time_step+1)=parameters.X.Lambda*y(:,time_step)+parameters.X.Phi*vec(z(:,time_step+max_p_q+1:-1:time_step+max_p_q-parameters.Z.p+2));
        end
    end;
    IRF((modelinfo.num_exog+modelinfo.num_endo)*(exog_select-1)+1:(modelinfo.num_exog+modelinfo.num_endo)*exog_select,:)=[y(:,end-IRF_horizon+1:end);
        z(:,end-IRF_horizon+1:end)];
end;
end
function x=vec(x)
x=x(:);
end