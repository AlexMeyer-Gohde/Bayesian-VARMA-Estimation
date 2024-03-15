function [ L ] = kalman_log_likelihood_alt_mv_first_difference(model_info, parameters)
% function [L ] = log_likelihood(parameters, model_info)
% This function computes the log-likelihood of an
% ARMA(p,q) model:
%
%  y[t] = phi(1)y[t-1] + ... + phi(p)y[t-p] + ...
%         + e[t] + theta(1)e[t-1] + ... + theta(q)e[t-q]
%
% REFERENCES:
%
%
% Wrtitte 4/17/2014 by Alexander Meyer-Gohde
% try

    
    p = parameters.Z.p;
    q = parameters.Z.q;
    T = length(model_info.data);
    ne=model_info.num_exog;
    nx=model_info.num_endo;
    ny=model_info.num_observables;
    % Build matrices of state space representation:
    %
    %      x[t+1] = A x[t] + R eps[t]
    %      y[t]   = Zx[t] + H eps[t]
    %
    Z=zeros(ny,ny+nx+ne*(q+max(p,1)));
    Z(:,end-ny:end)=-eye(ny);
    Z(:,1:nx)=parameters.X.Upsilon;
    A=zeros(ny+nx+ne*(q+max(p,1)),ny+nx+ne*(q+max(p,1)));
    A(1:nx,1:nx)=parameters.X.Lambda; %%% 5.5.15 AMG there appears to have been a typo here :-(
    R=zeros(ny+nx+ne*(q+max(p,1)),ne);
    R(1:nx,1:ne)=parameters.X.Phi(:,1:ne); %%% 8.5.15 AMG there appears to have been a typo here :-(
    R(nx+1:nx+ne,1:ne)=eye(ne);   %%% 8.5.15 AMG there appears to have been a typo here :-(
    %%%%%%%%%%%%%%%%%%Continue here%%%%%%%%%%%%%%%%%
    if q>1 && p>1
        temp_pA=[parameters.X.Phi(:,ne+1:end)]; temp_qA=[parameters.X.Theta(:,ne+1:end)];temp_p=eye(ne*(p-1));temp_q=eye(ne*(q-1));
    elseif q>1
        temp_pA=[];temp_qA=[parameters.X.Theta(:,ne+1:end)];temp_p=[];temp_q=eye(ne*(q-1));
    elseif p>1
        temp_pA=[parameters.X.Phi(:,ne+1:end)];temp_qA=[];temp_p=eye(ne*(p-1));temp_q=[];
    else
        temp_pA=[];temp_qA=[];temp_p=[];temp_q=[];
    end
    
    if q>0 && p>0
        if isempty(temp_pA)==0
            A(1:nx,nx+1:nx+ne*(p-1))=temp_pA;
        end
        if isempty(temp_qA)==0
            A(1:nx,nx+ne*(p)+1:nx+ne*(p)+ne*(q-1))=temp_qA;
        end
        A(1:nx,nx+1:nx+ne*(q+p))=A(1:nx,nx+1:nx+ne*(q+p))+parameters.X.Phi(:,1:ne)*[parameters.Z.P parameters.Z.Q];
        A(nx+1:nx+ne,nx+1:nx+ne*(q+p))=[parameters.Z.P parameters.Z.Q];
        A(nx+ne+1:nx+ne*p,nx+1:nx+ne*(p-1))=temp_p;
        A(nx+ne*(p+1)+1:nx+ne*(p+q),nx+ne*p+1:nx+ne*(p+q-1))=temp_q;
        R(1:nx,1:ne)=R(1:nx,1:ne)+parameters.X.Theta(:,1:ne);
        R(nx+ne*p+1:nx+ne*p+ne,1:ne)=eye(ne);        
    elseif p>0
        if isempty(temp_pA)==0
            A(1:nx,nx+1:nx+ne*(p-1))=temp_pA;
        end
        A(1:nx,nx+1:nx+ne*(p))=A(1:nx,nx+1:nx+ne*(p))+parameters.X.Phi(:,1:ne)*[parameters.Z.P ];
        A(nx+1:nx+ne,nx+1:nx+ne*(p))=[parameters.Z.P ];
        A(nx+ne+1:nx+ne*p,nx+1:nx+ne*(p-1))=temp_p;
    elseif q>0
        if isempty(temp_qA)==0
            A(1:nx,nx+ne+1:nx+ne*(q))=temp_qA;
        end
        A(1:nx,nx+ne+1:nx+ne+ne*(q))=A(1:nx,nx+ne+1:nx+ne+ne*(q))+parameters.X.Phi(:,1:ne)*[parameters.Z.Q];
        A(nx+1:nx+ne,nx+ne+1:nx+ne+ne*(q))=[parameters.Z.Q];
        A(nx+ne*(2)+1:nx+ne*(1+q),nx+ne*1+1:nx+ne*(q))=temp_q;
        R(nx+ne+1:nx+ne+ne,1:ne)=eye(ne);
        R(1:nx,1:ne)=R(1:nx,1:ne)+parameters.X.Theta(:,1:ne);
    end
    A(end-ny+1:end,1:nx)=parameters.X.Upsilon;
    
    Q=R*parameters.Z.Sigma*R';
    % Initialize state and covariance matrix
    xhat_tt  = zeros(ny+nx+ne*(q+max(p,1)),1);
    Sigma_tt = reshape( (eye((ny+nx+ne*(q+max(p,1)))^2)-kron(A,A) ) \ Q(:), ny+nx+ne*(q+max(p,1)),ny+nx+ne*(q+max(p,1)) );
    logsum = 0;
    sumsq  = 0;
    
    for i=0:T-1     % Start Kalman Filter Recursion
        xhat_t1t  = A*xhat_tt;
        Sigma_t1t = A*Sigma_tt*A' + Q;
        yhat_t1t  = Z*xhat_t1t;
        omega     = Z*Sigma_t1t*Z';
        delta_t1  = (Sigma_t1t*Z')/omega;
        innov     = model_info.data(:,i+1) - yhat_t1t;
        xhat_t1t1 = xhat_t1t + delta_t1*innov;
        Sigma_t1t1= Sigma_t1t - delta_t1*Z*Sigma_t1t;
        
        % Add likelihood terms:
        logsum = logsum + log(det(omega));
        sumsq  = sumsq  + (innov'/omega)*innov;

        % Update estimates;
        xhat_tt  = xhat_t1t1;
        Sigma_tt = Sigma_t1t1;
    end
    
    L = logsum + sumsq;
    L=-(1/2)*(ny*1.83787706640935+L);

% catch error
%     disp('Error in Kalman Filter');
end % function log_likelihood