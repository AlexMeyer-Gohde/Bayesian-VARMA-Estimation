    function [ L sigmahat] = kalmanLogLikelihood( phi, theta, p, q, sigmaE, y)
        % function [L sigmahat] = log_likelihood(coefs)
        % This function computes the negative of the log-likelihood of an
        % ARMA(p,q) model:
        %
        %  y[t] = phi(1)y[t-1] + ... + phi(p)y[t-p] + ... 
        %         + e[t] + theta(1)e[t-1] + ... + theta(q)e[t-q] 
        %
        % where phi    = coefs(1:p)
        %       theta  = coefs(p+1:p+q)
        
        %Taken from ARMA_MLE.M (WORLDBANK)        
        %
        %  y[t] = phi(1)y[t-1] + ... + phi(p)y[t-p] + e[t] + theta(1)e[t-1] + 
        %         ... + theta(q)e[t-q] 
        %
        % The algorithm computes the exact likelihood function using the Kalman
        % filter. 
        %
        %
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
        T = max(size(y));
        r = max(p,q+1);
        
%         phi   = coefs(1:p); 
%         theta = coefs(p+1:p+q);
        
        % Build matrices of state space representation:
        %
        %      x[t+1] = A x[t] + R eps[t+1]
        %      y[t]   = Z'x[t]
        %
        
        A = zeros(r,r); R = zeros(r,1); Z = zeros(r,1);
        A(1:p,1)=phi; A(1:r-1,2:r)=eye(r-1);
        R(1) = 1; R(2:q+1) = theta;
        Z(1) = 1; Q = R*R';
        
        % Initialize state and covariance matrix
        xhat_tt  = zeros(r,1);
        Sigma_tt = reshape( (eye(r^2)-kron(A,A) ) \ Q(:), r, r );
        
        logsum = 0;
        sumsq  = 0;
        
        for i=0:T-1     % Start Kalman Filter Recursion
        if ~isfinite(xhat_tt)
            disp('BS !!!');
        end;
            xhat_t1t  = A*xhat_tt;
            Sigma_t1t = A*Sigma_tt*A' + Q; 
            yhat_t1t  = Z'*xhat_t1t;
            %Commented out by DN 06.01.13 because tinkering
            omega     = Z'*Sigma_t1t*Z;
            delta_t1  = Sigma_t1t*Z/omega;
            innov     = y(i+1) - yhat_t1t;
            xhat_t1t1 = xhat_t1t + delta_t1*innov;
            Sigma_t1t1= Sigma_t1t - delta_t1*Z'*Sigma_t1t;
            
            % Add likelihood terms:
            logsum = logsum + log(omega);
            sumsq  = sumsq  + innov^2/omega;
            
            % Update estimates;
            xhat_tt  = xhat_t1t1; 
            Sigma_tt = Sigma_t1t1;
        end
        
        %Original Likelihood        
%         L = - logsum - T*sumsq;

        %Modified Likelihood (DN 06.01.13)        
        L = -logsum - T*log(sumsq) - T*log(sigmaE) -T\(sigmaE^2);
        
        if ~isfinite(L) | ~isreal(L)
            disp('BS !!!');
        end;
        sigmahat = sqrt(sumsq/T);
        
    end % function log_likelihood
