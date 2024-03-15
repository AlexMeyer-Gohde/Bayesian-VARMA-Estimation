function [loglikelihood, estimatedVariance] = evaluateLikelihoodARMAOnlyOLD( model_info, parameters )
%EVALUATE_LIKELIHOOD only univariate ARMA Model without DSGE Solver before
%for old format
%Input:         
%               model_info.data: vectorized data set
%               model_info.num_observations: number of observations
%               model_info.num_observables: number of observables =model_info.num_exog
%               parameters.W.Phi: AR component of W (num_exog x num_exog*p)
%               parameters.W.Theta: MA component of W (num_exog x num_exog*q)
%               
%
%Output:        loglikelihood: the value of the log-likelihood function with the
%                   given parameters and data set 

max_order=max(parameters.W.p+1,parameters.W.q+1);% max_order=max(max_order,3);

% grid_size=2*model_info.num_observations;
grid_size=max(2^10,2^(nextpow2(2*model_info.num_observations)+1));

grid_points=0:2*pi/grid_size:2*pi-2*pi/grid_size;
s_Y=zeros(model_info.num_exog^2,grid_size);

dd_1=repmat(grid_points, [ max_order 1]);
dd_2=repmat((0:max_order-1)',[1 length(grid_points)]);
dd=dd_1.*dd_2;

if model_info.num_exog>1
    EE_MINUS=zeros([max_order grid_size model_info.num_exog model_info.num_exog]); EE_MINUS(:,:,1:model_info.num_exog+1:model_info.num_exog^2)=repmat(exp(-1i*dd), [1 1 model_info.num_exog]);EE_MINUS=permute(EE_MINUS, [3 1 4 2]); EE_MINUS=reshape(EE_MINUS, 2*max_order, 2*grid_size);
    for n = 1 : grid_size
        try 
            if parameters.W.p==0
                W_e_minus=parameters.W.Theta*EE_MINUS(1:(parameters.W.q+1)*model_info.num_exog,(n-1)* model_info.num_exog+1:n* model_info.num_exog);
            else
                W_e_minus=(eye(model_info.num_exog)-parameters.W.Phi*EE_MINUS(model_info.num_exog+1:(parameters.W.p+1)*model_info.num_exog,(n-1)* model_info.num_exog+1:n* model_info.num_exog))\...
                    parameters.W.Theta*EE_MINUS(1:(parameters.W.q+1)*model_info.num_exog,(n-1)* model_info.num_exog+1:n* model_info.num_exog);
            end
        catch whateverhappened
            disp('boom');
        end;

        Temp=W_e_minus*parameters.Sigma*W_e_minus';
        s_Y(:,n)=Temp(:);
    end

    s_Y=ifft(transpose(s_Y)/(2*pi),'symmetric')*2*pi;
    s_Y=reshape(transpose(s_Y),[model_info.num_observables model_info.num_observables*grid_size]);

else
    EE_MINUS=exp(-1i*dd);
        try 
            if parameters.W.p==0
                W_e_minus=parameters.W.Theta*EE_MINUS(1:(parameters.W.q+1)*model_info.num_exog,:);
            else
                W_e_minus=(ones(1,grid_size)-parameters.W.Phi*EE_MINUS(model_info.num_exog+1:(parameters.W.p+1)*model_info.num_exog,:)).\(parameters.W.Theta*EE_MINUS(1:(parameters.W.q+1)*model_info.num_exog,:));
            end
        catch whateverhappened
            disp('boom');
        end;

        s_Y=parameters.Sigma*(W_e_minus.*conj(W_e_minus));

    s_Y=ifft(s_Y/(2*pi),'symmetric')*2*pi;
    s_Y=reshape(transpose(s_Y),[model_info.num_observables model_info.num_observables*grid_size]);
end
estimatedVariance = s_Y(1);
loglikelihood=block_levinson_loglik(model_info.data, s_Y(:,1:model_info.num_observations));

