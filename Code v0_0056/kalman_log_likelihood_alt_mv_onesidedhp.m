function [ L ] = kalman_log_likelihood_alt_mv_onesidedhp(model_info, parameters, HP_LAMBDA )
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
% Wrtitte 4/9/2014 by Alexander Meyer-Gohde
% try
    if HP_LAMBDA <= 0
        disp('Lambda nonsensical');
        error('Lambda nonsensical for Kalman Filter Likelihood');
    end;
    
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
    Z=zeros(ny,2*ny+nx+ne*(q+max(p,1)));
    Z(:,1:ny)=-eye(ny);
    Z(:,2*ny+1:2*ny+nx)=parameters.X.Upsilon;%This matrix selects the observables as a linear combo of enodgenous variables
    A=zeros(2*ny+nx+ne*(q+max(p,1)),2*ny+nx+ne*(q+max(p,1)));
    A(1:ny,1:ny)=2*eye(ny); A(1:ny,ny+1:2*ny)=-eye(ny);A(ny+1:2*ny,1:ny)=eye(ny);%%% 8.5.15 AMG there appears to have been a typo here :-(
    A(2*ny+1:2*ny+nx,2*ny+1:2*ny+nx)=parameters.X.Lambda;
    R=zeros(2*ny+nx+ne*(q+max(p,1)),ny+ne);
    R(1:ny,1:ny)=(1/HP_LAMBDA^(1/2))*eye(ny);
    R(2*ny+1:2*ny+nx,ny+1:ny+ne)=parameters.X.Phi(:,1:ne);
    R(2*ny+nx+1:2*ny+nx+ne,ny+1:ny+ne)=eye(ne);
    
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
            A(2*ny+1:2*ny+nx,2*ny+nx+1:2*ny+nx+ne*(p-1))=temp_pA;
        end
        if isempty(temp_qA)==0
            A(2*ny+1:2*ny+nx,2*ny+nx+ne*(p)+1:2*ny+nx+ne*(p)+ne*(q-1))=temp_qA;
        end
        A(2*ny+1:2*ny+nx,2*ny+nx+1:2*ny+nx+ne*(q+p))=A(2*ny+1:2*ny+nx,2*ny+nx+1:2*ny+nx+ne*(q+p))+parameters.X.Phi(:,1:ne)*[parameters.Z.P parameters.Z.Q];
        A(2*ny+nx+1:2*ny+nx+ne,2*ny+nx+1:2*ny+nx+ne*(q+p))=[parameters.Z.P parameters.Z.Q];
        A(2*ny+nx+ne+1:2*ny+nx+ne*p,2*ny+nx+1:2*ny+nx+ne*(p-1))=temp_p;
        A(2*ny+nx+ne*(p+1)+1:2*ny+nx+ne*(p+q),2*ny+nx+ne*p+1:2*ny+nx+ne*(p+q-1))=temp_q;
        R(2*ny+1:2*ny+nx,ny+1:ny+ne)=R(2*ny+1:2*ny+nx,ny+1:ny+ne)+parameters.X.Theta(:,1:ne);
        R(2*ny+nx+ne*p+1:2*ny+nx+ne*p+ne,ny+1:ny+ne)=eye(ne);        
    elseif p>0
        if isempty(temp_pA)==0
            A(2*ny+1:2*ny+nx,2*ny+nx+1:2*ny+nx+ne*(p-1))=temp_pA;
        end
        A(2*ny+1:2*ny+nx,2*ny+nx+1:2*ny+nx+ne*(p))=A(2*ny+1:2*ny+nx,2*ny+nx+1:2*ny+nx+ne*(p))+parameters.X.Phi(:,1:ne)*[parameters.Z.P ];
        A(2*ny+nx+1:2*ny+nx+ne,2*ny+nx+1:2*ny+nx+ne*(p))=[parameters.Z.P ];
        A(2*ny+nx+ne+1:2*ny+nx+ne*p,2*ny+nx+1:2*ny+nx+ne*(p-1))=temp_p;
    elseif q>0
        if isempty(temp_qA)==0
            A(2*ny+1:2*ny+nx,2*ny+nx+ne+1:2*ny+nx+ne*(q))=temp_qA;
        end
        A(2*ny+1:2*ny+nx,2*ny+nx+1:2*ny+nx+ne*(q))=A(2*ny+1:2*ny+nx,2*ny+nx+1:2*ny+nx+ne*(q))+parameters.X.Phi(:,1:ne)*[parameters.Z.Q];
        A(2*ny+1:2*ny+nx,2*ny+nx+ne+1:2*ny+nx+ne+ne*(q))=A(2*ny+1:2*ny+nx,2*ny+nx+ne+1:2*ny+nx+ne+ne*(q))+parameters.X.Phi(:,1:ne)*[parameters.Z.Q];
        A(2*ny+nx+1:2*ny+nx+ne,2*ny+nx+ne+1:2*ny+nx+ne+ne*(q))=[parameters.Z.Q];
        A(2*ny+nx+ne*(2)+1:2*ny+nx+ne*(1+q),2*ny+nx+ne*1+1:2*ny+nx+ne*(q))=temp_q;
        R(2*ny+nx+ne+1:2*ny+nx+ne+ne,ny+1:ny+ne)=eye(ne);
        R(2*ny+1:2*ny+nx,ny+1:ny+ne)=R(2*ny+1:2*ny+nx,ny+1:ny+ne)+parameters.X.Theta(:,1:ne);
    end
    
    Q=R*[eye(ny),zeros(ny,ne);zeros(ne,ny),parameters.Z.Sigma]*R';
    Q_2=Q(2*ny+1:end,2*ny+1:end);
    % Initialize state and covariance matrix
    xhat_tt_2  = zeros(nx+ne*(q+max(p,1)),1);
    Sigma_tt_2 = reshape( (eye((nx+ne*(q+max(p,1)))^2)-kron(A(2*ny+1:end,2*ny+1:end),A(2*ny+1:end,2*ny+1:end)) ) \ Q_2(:), nx+ne*(q+max(p,1)),nx+ne*(q+max(p,1)) );
    Sigma_tt_1=1e5*eye(2*ny);
    xhat_tt_1=[2*model_info.data(:,1)-model_info.data(:,2);3*model_info.data(:,1)-2*model_info.data(:,2)];
    xhat_tt=[xhat_tt_1;xhat_tt_2];
    Sigma_tt=[Sigma_tt_1, zeros(2*ny,nx+ne*(q+max(p,1))); zeros(nx+ne*(q+max(p,1)),2*ny),Sigma_tt_2];
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