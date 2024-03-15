function COR = DSGE_COR(model_info, parameters,COR_horizon,HP_flag)
model_info.num_observables=model_info.num_endo;
parameters.X.Upsilon=eye(model_info.num_endo);
max_order=max(parameters.Z.p+1,parameters.Z.q+1);% max_order=max(max_order,3);

grid_size=max(2^10,2^(nextpow2(2*COR_horizon)+1));
grid_points=0:2*pi/grid_size:2*pi-2*pi/grid_size;
s_Y=zeros(model_info.num_observables^2,grid_size);
if HP_flag
    HP_LAMBDA=HP_flag;
    hp = 4*HP_LAMBDA*(1 - cos(grid_points)).^2 ./ (1 + 4*HP_LAMBDA*(1 - cos(grid_points)).^2);
end

dd_1=repmat(grid_points, [ max_order 1]);
dd_2=repmat((0:max_order-1)',[1 length(grid_points)]);
dd=dd_1.*dd_2;

if model_info.num_exog>1 || model_info.num_endo>1
    EE_MINUS=zeros([max_order grid_size model_info.num_exog model_info.num_exog]); EE_MINUS(:,:,1:model_info.num_exog+1:model_info.num_exog^2)=repmat(exp(-1i*dd), [1 1 model_info.num_exog]);EE_MINUS=permute(EE_MINUS, [3 1 4 2]);EE_MINUS=reshape(EE_MINUS, [model_info.num_exog*max_order, model_info.num_exog*grid_size]);
    for n = 1 : grid_size
        try
            if parameters.Z.p==0
                if parameters.Z.q > 0 %p == 0, q > 0
                    W_e_minus=parameters.X.Upsilon*((eye(model_info.num_endo)-parameters.X.Lambda*exp(-1i*grid_points(n)))\...
                        ((parameters.X.Phi*([eye(model_info.num_exog) parameters.Z.Q]*EE_MINUS(1:(parameters.Z.q+1)*model_info.num_exog,(n-1)* model_info.num_exog+1:n* model_info.num_exog)))...
                        +parameters.X.Theta*EE_MINUS(1:(parameters.Z.q)*model_info.num_exog,(n-1)* model_info.num_exog+1:n* model_info.num_exog)));
                else %p == 0, q == 0
                    W_e_minus=parameters.X.Upsilon*((eye(model_info.num_endo)-parameters.X.Lambda*exp(-1i*grid_points(n)))\...
                        (parameters.X.Phi*([eye(model_info.num_exog) parameters.Z.Q]*EE_MINUS(1:(parameters.Z.q+1)*model_info.num_exog,(n-1)* model_info.num_exog+1:n* model_info.num_exog))));
                end;
            else %p > 0
                if parameters.Z.q > 0 %p > 0, q > 0
                    W_e_minus=parameters.X.Upsilon*((eye(model_info.num_endo)-parameters.X.Lambda*exp(-1i*grid_points(n)))\...
                        (parameters.X.Phi*EE_MINUS(1:(parameters.Z.p)*model_info.num_exog,(n-1)* model_info.num_exog+1:n* model_info.num_exog)/...
                        (eye(model_info.num_exog)-parameters.Z.P*EE_MINUS(model_info.num_exog+1:(parameters.Z.p+1)*model_info.num_exog,(n-1)* model_info.num_exog+1:n* model_info.num_exog))*...
                        ([eye(model_info.num_exog) parameters.Z.Q])*EE_MINUS(1:(parameters.Z.q+1)*model_info.num_exog,(n-1)* model_info.num_exog+1:n* model_info.num_exog)...
                        +parameters.X.Theta*EE_MINUS(1:(parameters.Z.q)*model_info.num_exog,(n-1)* model_info.num_exog+1:n* model_info.num_exog)));
                else  %p > 0, q == 0
                    W_e_minus=parameters.X.Upsilon*((eye(model_info.num_endo)-parameters.X.Lambda*exp(-1i*grid_points(n)))\...
                        ((parameters.X.Phi*EE_MINUS(1:(parameters.Z.p)*model_info.num_exog,(n-1)* model_info.num_exog+1:n* model_info.num_exog))/...
                        (eye(model_info.num_exog)-parameters.Z.P*EE_MINUS(model_info.num_exog+1:(parameters.Z.p+1)*model_info.num_exog,(n-1)* model_info.num_exog+1:n* model_info.num_exog))*...
                        (([eye(model_info.num_exog) parameters.Z.Q])*EE_MINUS(1:(parameters.Z.q+1)*model_info.num_exog,(n-1)* model_info.num_exog+1:n* model_info.num_exog))));
                end;
            end
        catch whateverhappened
            disp('boom');
        end;
        Temp=W_e_minus*parameters.Z.Sigma*W_e_minus';
        if HP_flag
            s_Y(:,n)=hp(n)^2*Temp(:);
        else
            s_Y(:,n)=Temp(:);
        end
    end
    s_Y=ifft(transpose(s_Y)/(2*pi),'symmetric')*2*pi;
    
    s_Y=reshape(transpose(s_Y), [model_info.num_observables model_info.num_observables*grid_size]);
    
else
    EE_MINUS=exp(-1i*dd);
    try
        if parameters.Z.p==0
            W_e_minus=parameters.X.Upsilon*((ones(1,grid_size)-parameters.X.Lambda*exp(-1i*grid_points)).\...
                (parameters.Z.Q*EE_MINUS(1:(parameters.Z.q+1)*model_info.num_exog,:)...
                +parameters.X.Theta*EE_MINUS(1:(parameters.Z.q)*model_info.num_exog,:)));
        else
            W_e_minus=parameters.X.Upsilon*((ones(1,grid_size)-parameters.X.Lambda*exp(-1i*grid_points)).\...
                (parameters.X.Phi*EE_MINUS(1:(parameters.Z.p)*model_info.num_exog,:)./...
                (ones(1,grid_size)-parameters.Z.P*EE_MINUS(model_info.num_exog+1:(parameters.Z.p+1)*model_info.num_exog,:)).*...
                parameters.Z.Q*EE_MINUS(1:(parameters.Z.q+1)*model_info.num_exog,:)...
                +parameters.X.Theta*EE_MINUS(1:(parameters.Z.q)*model_info.num_exog,:)));
        end
    catch whateverhappened
        disp('boom');
    end;
    if HP_flag
        s_Y=parameters.Z.Sigma*hp.^2.*(W_e_minus.*conj(W_e_minus));
    else
        s_Y=parameters.Z.Sigma*(W_e_minus.*conj(W_e_minus));
    end
   
    
    s_Y=ifft(s_Y/(2*pi),'symmetric')*2*pi;
    s_Y=reshape(transpose(s_Y), [model_info.num_observables model_info.num_observables*grid_size]);
end
COR=[s_Y(:,end-COR_horizon*model_info.num_endo+1:end) s_Y(:,1:(COR_horizon+1)*model_info.num_endo)];
COR_Factor=repmat(1./(diag(s_Y(:,1:model_info.num_endo)).^(1/2)*diag(s_Y(:,1:model_info.num_endo)).^(1/2)'), [1 2*COR_horizon+1]);
COR=[diag(s_Y(:,1:model_info.num_endo)) COR_Factor.*COR];