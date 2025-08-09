%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set_initial_values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code used in Piecuch 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define initial process and parameter values as described in methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by CGP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%v

% mean parameters
mu=normrnd(HP.eta_tilde_mu,sqrt(HP.zeta_tilde_mu_2));
nu=normrnd(HP.eta_tilde_nu,sqrt(HP.zeta_tilde_nu_2));
alpha=[];

% variance parameters
pi_2=min([1 1/randraw('gamma', [0,1/HP.chi_tilde_pi_2,HP.xi_tilde_pi_2], [1,1])]); % use min to prevent needlessly large values
delta_2=min([1 1/randraw('gamma', [0,1/HP.chi_tilde_delta_2,HP.xi_tilde_delta_2], [1,1])]); % use min to prevent needlessly large values
sigma_2=min([1 1/randraw('gamma', [0,1/HP.chi_tilde_sigma_2,HP.xi_tilde_sigma_2], [1,1])]); % use min to prevent needlessly large values
tau_2=min([1 1/randraw('gamma', [0,1/HP.chi_tilde_tau_2,HP.xi_tilde_tau_2], [1,1])]); % use min to prevent needlessly large values
gamma_2=min([1 1/randraw('gamma', [0,1/HP.chi_tilde_gamma_2,HP.xi_tilde_gamma_2], [1,1])]); % use min to prevent needlessly large values

% inverse length scale parameters
phi=exp(normrnd(HP.eta_tilde_phi,sqrt(HP.zeta_tilde_phi_2)));
lambda=exp(normrnd(HP.eta_tilde_lambda,sqrt(HP.zeta_tilde_lambda_2)));

% spatial fields
a=(mvnrnd(zeros(M,1),gamma_2*eye(M)))';        
b=(mvnrnd(mu*ones(N,1),pi_2*exp(-lambda*D)))';        
l=(mvnrnd(nu*ones(M,1),tau_2*eye(M)))';        

% AR(1) parameter
r=HP.u_tilde_r+(HP.v_tilde_r-HP.u_tilde_r )*rand(1);

% process
y_0=zeros(N,1);
y=zeros(N,K);
