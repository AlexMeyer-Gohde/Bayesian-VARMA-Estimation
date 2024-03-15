function [model_info,parameters]=solve_recursive_dsge(model_info,parameters)

% %%%%%%%%%%%%%%%%%%%%%%%%%
% 
% solve_recursive_dsge
% This function returns the unique, stable, recursive solution (if it exists)
% x_t=\Lambda x_{t-1}+\sum_{i=0}^{p-1}\Phi_iz_{t-i}+\sum_{i=0}^q\Theta_i\epsilon_{t-i}
% to the problem
% 0=E_t\left[Ax_{t+1}+Bx_t+Cx_{t-1}+Dz_t\right]
% with
% z_t=\sum_{i=1}^{p}P_iz_{t-i}+\sum_{i=0}^qQ_i\epsilon_{t-i}
% and
% \epsilon\sim\left(0,\Omega)
% 
% inputs
% 
% 
% outputs
% 
% 
% created 10/21/2013
% 
% %%%%%%%%%%%%%%%%%%%%%%%

try
    [parameters.X.Lambda,~] = solab([eye(model_info.num_endo) zeros(model_info.num_endo); zeros(model_info.num_endo) -parameters.X.A ],[ zeros(model_info.num_endo)  eye(model_info.num_endo) ;parameters.X.C parameters.X.B],model_info.num_endo);
    A_tilde=(parameters.X.A*parameters.X.Lambda+parameters.X.B)\parameters.X.A;
    if parameters.Z.p>1
        [Phi_0]=solve_p_generalized_sylvester(model_info,parameters,A_tilde);
        A_tilde_Phi_0=A_tilde*Phi_0;
        parameters.X.Phi=zeros(model_info.num_endo,model_info.num_exog*parameters.Z.p);
        for i=1:parameters.Z.p-1
            parameters.X.Phi(:,model_info.num_exog*(parameters.Z.p-1-i)+1:model_info.num_exog*(parameters.Z.p-i))=-A_tilde_Phi_0*parameters.Z.P(:,(parameters.Z.p-i)*model_info.num_exog+1:(parameters.Z.p-i+1)*model_info.num_exog)-A_tilde*parameters.X.Phi(:,(parameters.Z.p-i)*model_info.num_exog+1:(parameters.Z.p-i+1)*model_info.num_exog);
        end
        parameters.X.Phi=parameters.X.Phi(:,1:end-model_info.num_exog);parameters.X.Phi=[Phi_0 parameters.X.Phi];
    elseif parameters.Z.p==1
        [Phi_0]=solve_p_generalized_sylvester(model_info,parameters,A_tilde);
        A_tilde_Phi_0=A_tilde*Phi_0;
        parameters.X.Phi=Phi_0;
    else %if parameters.Z.p==0
        Phi_0=-(parameters.X.A*parameters.X.Lambda+parameters.X.B)\parameters.X.D;
        A_tilde_Phi_0=A_tilde*Phi_0;
        parameters.X.Phi=Phi_0;
    end
    if parameters.Z.q>0
        parameters.X.Theta=zeros(model_info.num_endo,model_info.num_exog*(parameters.Z.q+1));
        for i=1:parameters.Z.q
            parameters.X.Theta(:,model_info.num_exog*(parameters.Z.q-i)+1:model_info.num_exog*(parameters.Z.q+1-i))=-A_tilde_Phi_0*parameters.Z.Q(:,(parameters.Z.q-i)*model_info.num_exog+1:(parameters.Z.q-i+1)*model_info.num_exog)-A_tilde*parameters.X.Theta(:,(parameters.Z.q-i+1)*model_info.num_exog+1:(parameters.Z.q-i+2)*model_info.num_exog);
        end
        parameters.X.Theta=parameters.X.Theta(:,1:end-model_info.num_exog);
    end
catch err
    disp('boom in solver');
end

