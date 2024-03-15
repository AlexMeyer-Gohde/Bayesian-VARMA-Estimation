function [ L ] = kalman_log_likelihood(parameters, model_info )
        % function [L ] = log_likelihood(parameters, model_info)
        % This function computes the log-likelihood of an
        % ARMA(p,q) model:
        %
        %  y[t] = phi(1)y[t-1] + ... + phi(p)y[t-p] + ... 
        %         + e[t] + theta(1)e[t-1] + ... + theta(q)e[t-q] 
        %
	% REFERENCES: 
	% 1) E.J. Hannan and L. Kavalieris. "A Method for Autoregressive-Moving
	%       Average Estimation" Biometrika, Vol 71, No2, Aug 1984.
	% 2) E.J. Hannan and A.J. McDougall. "Regression Procedures for ARMA
	%       Estimation" Journal of the American Statistical Association, Vol
	%       83, No 409, June 1988.
	% 3) R.H. Jones : Maximum Likelihood Fitting of ARMA Models to
	%       Time Series With Missing Observations" Technometrics, Vol 22, No3, 
	%       August 1980.
	% 
	% Written by Constantino Hevia. August 2008
	%
	% Modified 1/14/2014 by Alexander Meyer-Gohde

	r = max(parameters.W.p,parameters.W.q+1);
        
        % Build matrices of state space representation:
        %
        %      x[t+1] = A x[t] + R eps[t+1]
        %      y[t]   = Z'x[t]
        %
        
        A = zeros(r,r); R = zeros(r,1); Z = zeros(r,1);
        A(1:parameters.W.p,1)=parameters.W.Phi(:); A(1:r-1,2:r)=eye(r-1);
         R(1:parameters.W.q+1) = parameters.W.Theta(:);
        Z(1) = 1; Q = parameters.Sigma*R*R';
        
        
        
        % Initialize state and covariance matrix
        xhat_tt  = zeros(r,1);
%         Sigma_tt = reshape( (eye(r^2)-kron(A,A) ) \ Q(:), r, r );
        
        
        [Sigma_t1t,L,G] = dare(A',Z,Q,0);
      omega     = Z'*Sigma_t1t*Z;
      delta_t1  = Sigma_t1t*Z/omega;
%       Sigma_t1t1= Sigma_t1t - delta_t1*Z'*Sigma_t1t;
        logsum = 0;
        sumsq  = 0;
        T = length(model_info.data);
        for i=0:T-1     % Start Kalman Filter Recursion
            xhat_t1t  = A*xhat_tt;
           % Sigma_t1t = A*Sigma_tt*A' + Q; 
            yhat_t1t  = Z'*xhat_t1t;
           % omega     = Z'*Sigma_t1t*Z;
            %delta_t1  = Sigma_t1t*Z/omega;
            innov     = model_info.data(i+1) - yhat_t1t;
            xhat_t1t1 = xhat_t1t + delta_t1*innov;
            %Sigma_t1t1= Sigma_t1t - delta_t1*Z'*Sigma_t1t;
            
            % Add likelihood terms:
            logsum = logsum + log(omega);
            sumsq  = sumsq  + innov^2/omega;
            
            % Update estimates;
            xhat_tt  = xhat_t1t1; 
        end
        
        L = logsum + sumsq;
        L=-(1/2)*(T*1.83787706640935+L);
    end % function log_likelihood