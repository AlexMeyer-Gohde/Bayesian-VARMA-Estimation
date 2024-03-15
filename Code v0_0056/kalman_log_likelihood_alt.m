function [ L ] = kalman_log_likelihood_alt(parameters, model_info )
        % function [L ] = log_likelihood(parameters, model_info)
        % This function computes the log-likelihood of an
        % ARMA(p,q) model:
        %
        %  y[t] = phi(1)y[t-1] + ... + phi(p)y[t-p] + ... 
        %         + e[t] + theta(1)e[t-1] + ... + theta(q)e[t-q] 
        %
	% REFERENCES: 
	% 1) Pearlman (1980) An algorithm for the exact likelihood of a
	% high-order autoregressive-moving average process
	%
	% Wrtitte 3/11/2014 by Alexander Meyer-Gohde

	r = max(parameters.W.p,parameters.W.q);
        
    p = parameters.W.p;
    q = parameters.W.q;
    T = length(model_info.data);
        % Build matrices of state space representation:
        %
        %      x[t+1] = A x[t] + R eps[t]
        %      y[t]   = Zx[t] + H eps[t]
        %
        

        if (p==0) && (q==0)
            logsum = 0;
            sumsq  = 0;
            omega = parameters.Sigma;
            for i=0:T-1     % Start Kalman Filter Recursion
                innov     = model_info.data(i+1);

                % Add likelihood terms:
                logsum = logsum + log(omega);
                sumsq  = sumsq  + innov^2/omega;
            end
        else
            A = zeros(r,r); R_1 = zeros(r,1); R_2=R_1;Z = zeros(1,r); H=1;
            A(1:p,1)=parameters.W.Phi(:); A(1:r-1,2:r)=eye(r-1);
            R_1(1:q) = parameters.W.Theta(2:end);R_2(1:p)=parameters.W.Phi(:);R=R_1+R_2;
            Z(1) = 1; Q = parameters.Sigma*(R*R');
            
            C=R*parameters.Sigma*H';%Covariance between measurement and observation noise
            H_Sigma_H_trans_plus_ZC_plus_ZC_trans=H*parameters.Sigma*H'+Z*C+C'*Z';
        
            C_trans=C';
            % Initialize state and covariance matrix
            xhat_tt  = zeros(r,1);
            Sigma_tt = reshape( (eye(r^2)-kron(A,A) ) \ Q(:), r, r );

            logsum = 0;
            sumsq  = 0;

            for i=0:T-1     % Start Kalman Filter Recursion
                xhat_t1t  = A*xhat_tt;
                Sigma_t1t = A*Sigma_tt*A' + Q; 
                yhat_t1t  = Z*xhat_t1t;
                omega     = Z*Sigma_t1t*Z'+H_Sigma_H_trans_plus_ZC_plus_ZC_trans;
                delta_t1  = (Sigma_t1t*Z'+C)/omega;
                innov     = model_info.data(i+1) - yhat_t1t;
                xhat_t1t1 = xhat_t1t + delta_t1*innov;
                Sigma_t1t1= Sigma_t1t - delta_t1*Z*Sigma_t1t-delta_t1*C_trans;

                % Add likelihood terms:
                logsum = logsum + log(omega);
                sumsq  = sumsq  + innov^2/omega;

                % Update estimates;
                xhat_tt  = xhat_t1t1; 
                Sigma_tt = Sigma_t1t1;
            end
        end;
        L = logsum + sumsq;
        L=-(1/2)*(1.83787706640935+L);
    end % function log_likelihood