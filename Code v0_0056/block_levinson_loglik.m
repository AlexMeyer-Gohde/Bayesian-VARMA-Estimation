function [x] = block_levinson_loglik(y, L, varargin)
%Block Levinson recursion for efficiently determining the log-likelihood.
%   BLOCK_LEVINSON(Y, L) yields
%   x=-0.5*(d*N*ln(2*pi)+ln(det(G))+y'*inv(G)*y)
%   where G is a symmetric Block-Toeplitz matrix given by
%
%G=[L(1) L(2)' ...  ... L(N)';
%   L(2) L(1) L(2)' ... L(N-1)';
%   ...
%   L(N) L(N-1) ...     L(1)]
%   G is not inputed in full because the row vector [L(1) L(2)... L(N)] contains all the necessary information.
%   The row vector [L(1) L(2)... L(N)] contains the first N entries of the
%   autocovariance generating function of the underlying stochastic process
%
%   INPUT:
%           y -- 1xN vector of observations
%           L -- dxdN or dNxd block vector of autocovariances
%
%   Created by Alexander Meyer-Gohde
%   Created on Feb. 24th, 2009
%
%
%   Akaike, Hirotugu (1973), "Block Toeplitz Matrix Inversion", SIAM J.
%   Appl. Math., 24 (2), pp. 234-241


sz = size(L);
ln2pi = 1.83787706640935;
if nargin ==2
    if sz(1)>=sz(2) %block column input
        
        d=sz(2);
        N = sz(1) / d;
        Pp=L(1:d,1:d);%Get first block
        
        r=L(d+1:N*d,:);%Get first block column
        
        a_hat=L(d+1:N*d,:)';
        a_hat = reshape( a_hat, [ d d (N-1) ] );
        a_hat = a_hat(:,:,(N-1):-1:1);
        a_hat = permute( a_hat, [ 1 3 2 ] );
        a_hat = reshape( a_hat, [ d*(N-1) d ] );
        
    else %block row input
        
        d = sz(1);
        N = sz(2) / d;
        Pp=L(1:d,1:d); %Get first block
        
        r = reshape( L(:,d+1:N*d), [ d d (N-1) ] );
        r = permute( r, [ 1 3 2 ] );
        r = reshape( r, [ d*(N-1) d ] ); %Get first block column
        
        a_hat=L(:,d+1:N*d)';
        a_hat = reshape( a_hat, [ d (N-1) d ] );
        a_hat = a_hat(:,(N-1):-1:1,:);
        a_hat = reshape( a_hat, [ d*(N-1) d ] ); %Get last block column
        
    end
else %if third input argument is present, use inverse FFT population moments format [L(1) L(2) ... L(N+1) L(N)' L(N-1)' ... L(2)']
    d = sz(1);
    N = sz(2) / (2*d);
    
    Pp=L(1:d,1:d);%Get first block
    
    r = reshape( L(:,d+1:N*d), [ d d (N-1) ] );
    r = permute( r, [ 1 3 2 ] );
    r = reshape( r, [ d*(N-1) d ] );%Get first block column
    
    a_hat = reshape( L(:,d*(N+1)+1:d*2*N), [ d d (N-1) ] );
    a_hat = permute( a_hat, [ 1 3 2 ] );
    a_hat = reshape( a_hat, [ d*(N-1) d ] );%Get last block column
    
    
end


s_inv = Pp;                     % initial value for s_inv
q_inv = s_inv;                  % initial value for q_inv
[chol_s_inv,chol_error]=chol(s_inv);
if chol_error~=0         % Catch a failure to invert s_inv and/or q_inv
    x=-10^50;
    return
end
half_lndet=sum(log(diag(chol_s_inv)));% first component of 0.5*ln(det(G))
s=chol_s_inv\(chol_s_inv'\eye(d,d));% initial value for s
q=s;                            % initial value for q
x = y(1:d)'*s*y(1:d);           % first component of y'*inv(G)*y

si=-r(1:d,:)*q;
qf=-a_hat((N-2)*d+1:(N-1)*d,:)*s;
for n = 2:N-1
    s_inv = (eye(d,d)-si(:,1:d)*qf(:,(n-2)*d+1:(n-1)*d))*s_inv;
    [chol_s_inv,chol_error]=chol(s_inv);
    if chol_error~=0         % Catch a failure to invert s_inv and/or q_inv
        x=-10^50;
        return
    end
    half_lndet=half_lndet+sum(log(diag(chol_s_inv)));
    q_inv = (eye(d,d)-qf(:,(n-2)*d+1:(n-1)*d)*si(:,1:d))*q_inv;
    [chol_q_inv,chol_error]=chol(q_inv);
    if chol_error~=0         % Catch a failure to invert s_inv and/or q_inv
        x=-10^50;
        return
    end
    s=chol_s_inv\(chol_s_inv'\eye(d,d));
    q=chol_q_inv\(chol_q_inv'\eye(d,d));
    sy=s*y((n-1)*d+1:n*d);
    siY=si*y(1:(n-1)*d);
    x=x+y((n-1)*d+1:n*d)'*sy+siY'*sy+sy'*siY+siY'*s*siY;
    
    si_old=si;
    si=[zeros(d) si]-(r((n-1)*d+1:(n)*d,:)+si*r(1:(n-1)*d,:))*q*[eye(d) qf];
    qf=[qf zeros(d)]-(a_hat((N-n-1)*d+1:(N-n)*d,:)+qf*a_hat((N-n)*d+1:(N-1)*d,:))*s*[si_old eye(d)];
    
    
end

s_inv = (eye(d,d)-si(:,1:d)*qf(:,(N-2)*d+1:(N-1)*d))*s_inv;
[chol_s_inv,chol_error]=chol(s_inv);
if chol_error~=0         % Catch a failure to invert s_inv and/or q_inv
    x=-10^50;
    return
end
half_lndet=half_lndet+sum(log(diag(chol_s_inv)));
s=chol_s_inv\(chol_s_inv'\eye(d,d));
sy=s*y((N-1)*d+1:N*d);
siY=si*y(1:(N-1)*d);
x=x+y((N-1)*d+1:N*d)'*sy+siY'*sy+sy'*siY+siY'*s*siY;
x=-0.5*(N*d)*ln2pi-0.5*x-half_lndet;