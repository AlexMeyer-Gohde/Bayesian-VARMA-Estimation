function [model_info, parameters] = exampl1_extract()

N_bar     = 1.0/3;  % Steady state employment is a third of total time endowment
Z_bar     = 1; % Normalization
rho       = .36; % Capital share
delta     = .025; % Depreciation rate for capital
R_bar     = 1.01; % One percent real interest per quarter
eta       =  1.0; % constant of relative risk aversion = 1/(coeff. of intertemporal substitution)
psi       = .95; % autocorrelation of technology shock
sigma_eps = 0.626066984205494; % Standard deviation of technology shock.  Units: Percent.

% Calculating the steady state:

betta   = 1.0/R_bar;  % Discount factor beta
YK_bar  = (R_bar + delta - 1)/rho;  % = Y_bar / K_bar
K_bar   = (YK_bar / Z_bar)^(1.0/(rho-1)) * N_bar;
I_bar   = delta * K_bar;
Y_bar   = YK_bar * K_bar;
C_bar   = Y_bar - delta*K_bar;
A       =  C_bar^(-eta) * (1 - rho) * Y_bar/N_bar; % Parameter in utility function

% Declaring the matrices. 


VARNAMES = ['capital    ',
            'consumption',
            'output     ',
            'labor      ',
            'interest   ',
            'investment ',
            'technology '];

% Translating into coefficient matrices.  
% The equations are, conveniently ordered:
% 1) 0 = - I i(t) - C c(t) + Y y(t)
% 2) 0 = I i(t) - K k(t) + (1-delta) K k(t-1)
% 3) 0 = rho k(t-1) - y(t) + (1-rho) n(t) + z(t)
% 4) 0 = -eta c(t) + y(t) - n(t)
% 5) 0 = - rho Y/K k(t-1) + rho Y/K y(t) - R r(t)
% 6) 0 = E_t [ - eta c(t+1) + r(t+1) + eta c(t) ]
% 7) z(t+1) = psi z(t) + epsilon(t+1)
% CHECK: 7 equations, 7 variables.
%
% Endogenous state variables "x(t)": k(t)
% Endogenous other variables "y(t)": c(t), y(t), n(t), r(t), i(t)
% Exogenous state variables  "z(t)": z(t).
% Switch to that notation.  Find matrices for format
% 0 = AA x(t) + BB x(t-1) + CC y(t) + DD z(t)
% 0 = E_t [ FF x(t+1) + GG x(t) + HH x(t-1) + JJ y(t+1) + KK y(t) + LL z(t+1) + MM z(t)]
% z(t+1) = NN z(t) + epsilon(t+1) with E_t [ epsilon(t+1) ] = 0,

% for k(t):
AA = [ 0
       - K_bar
       0
       0
       0 ];

% for k(t-1):
BB = [ 0
       (1-delta)*K_bar
       rho
       0
       - rho * YK_bar ];

%Order:   consumption  output      labor     interest  investment
CC = [    -C_bar,      Y_bar,      0,        0,        -I_bar % Equ. 1)
          0,           0,          0,        0,        I_bar  % Equ. 2)
          0,           -1,         1-rho,    0,        0      % Equ. 3)      
          -eta,        1,          -1,       0,        0      % Equ. 4)
          0,           rho*YK_bar, 0,        - R_bar,  0 ];   % Equ. 5)

DD = [ 0
       0
       1
       0
       0 ];

FF = [ 0 ];

GG = [ 0 ];

HH = [ 0 ];

JJ = [ -eta,  0,  0,  1,  0];

KK = [ eta,   0,  0,  0,  0];

LL = [ 0 ];

MM = [ 0 ];


ny=5;
nx=1;
model_info.num_exog=1;
model_info.num_endo=ny+nx;
% model_info.num_observations = length(data);
model_info.num_observables = 1;
% [a b] = size(data);
% if b>a
%     model_info.data = transpose(data);
% else
%     model_info.data = data;
% end;

parameters.X.A=[zeros(ny,nx+ny);FF JJ];
parameters.X.B=[AA CC; GG KK];
parameters.X.C=[BB zeros(ny, ny); HH zeros(nx, ny)];
parameters.X.D=[DD; MM];
parameters.X.Upsilon=[0 0 1 0 0 0];


parameters.Z.p=3 %ar order here
parameters.Z.q=3; %ma order here
parameters.Z.P=[0.9 -0.23 0.3]; %ar parameters here
parameters.Z.Q=[ -0.4 0.6 -0.5]; %ma parameters here
parameters.Z.Sigma = sigma_eps^2;

end
