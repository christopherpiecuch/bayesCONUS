%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function HP = set_hyperparameters(N,K,M,LON,LAT,TGR_DATA,priorChoice)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code used in Piecuch 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates a structure filled with the hyperparameters (parameters of the
% prior distributions) following the rationale outlined in
% "Piecuch_model_description.pdf"
% INPUT: 
%   N           Number of target locations (here 211)
%   K           Number of time steps (here 118, 1900-2017)
%   M           Number of tide gauges (here 53)
%   LON         Longitudes of target locations
%   LAT         Latitudes of target locations
%   TGR_DATA    Array of tide gauge RSL time series
%   priorChoice Flag to choose wide/normal/narrow priors on inverse ranges
%               (value can be 1,2,3 ... see below usage)
% OUTPUT:
%   HP          Structure of hyperparameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by CGP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function HP = set_hyperparameters(N,K,M,LON,LAT,TGR_DATA,priorChoice)

t=1:K;
m=nan(M,1);
s=nan(M,1);
r=nan(M,1);
e=nan(M,1);
y0=nan(M,1);
l=nan(M,1);
for n=1:M
    y=[]; 
    y=squeeze(TGR_DATA(n,:)); 
    i=[]; 
    i=find(~isnan(y));
        p=[]; q=[]; [p,q]=polyfit(t(i),y(i),1);
        w=(inv(q.R)*inv(q.R)')*q.normr^2/q.df; 
        m(n)=p(1);
        s(n)=w(1,1);
        [a b]=aryule(y(i)-p(1)*t(i)-p(2),1);
        d(n)=std(y(i)-p(1)*t(i)-p(2));
        r(n)=-a(2);
        e(n)=sqrt(b);
        l(n)=p(1)*mean(t)+p(2);
        y0(n)=p(2)-l(n);
    clear y i p w a b
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set variance inflation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var_infl=5^2;
var_infl2=10^2;
var0_infl=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the hyperparameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y_0
HP.eta_tilde_y_0 = y0;
HP.zeta_tilde_y_0_2 = diag((d+0.05).^2);


% r
HP.u_tilde_r                = 0.1;% Lower bound of r prior
HP.v_tilde_r                = 0.9; % Upper bound of r prior

% mu
HP.eta_tilde_mu             = nanmean(m);%+nanmean(u);% Mean of mu prior
HP.zeta_tilde_mu_2          = var_infl2*(nanvar(m));%+nanvar(u)); % Variance of mu prior

% nu
HP.eta_tilde_nu             = nanmean(l); % Mean of nu prior
HP.zeta_tilde_nu_2          = var_infl*nanvar(l); % Variance of nu prior

% pi_2
HP.xi_tilde_pi_2            = 1/2; % Shape of pi_2 prior
HP.chi_tilde_pi_2           = 1/2*nanvar(m); % Inverse scale of pi_2 prior

% delta_2
HP.xi_tilde_delta_2         = 1/2; % Shape of delta_2 prior
HP.chi_tilde_delta_2        = 1/2*1e-4; % Guess (1 cm)^2 error variance

% sigma_2
HP.xi_tilde_sigma_2         = 1/2; % Shape of sigma_2 prior
HP.chi_tilde_sigma_2        = 1/2*nanmean(e.^2); % Inverse scale of sigma_2 prior

% tau_2
HP.xi_tilde_tau_2           = 1/2; % Shape of tau_2 prior
HP.chi_tilde_tau_2          = 1/2*nanvar(l); % Inverse scale of tau_2 prior

% gamma_2
HP.xi_tilde_gamma_2         = 1/2; % Shape of tau_2 prior
HP.chi_tilde_gamma_2        = 1/2*(1e-3)^2; % Guess (1 mm/yr)^2 error variance

priorRangeMean              = -6.9;
% priorChoice = 1 --> (0.35)^2 -- middle of the road
% priorChoice = 2 --> (0.70)^2 -- wide priors
% priorChoice = 3 --> (0.10)^2 -- narrow priors

if priorChoice==1
    priorRangeVar = (0.35)^2;
elseif priorChoice==2
    priorRangeVar = (0.70)^2;
elseif priorChoice==3
    priorRangeVar = (0.10)^2;
else
    error('Invalid choice for priorChoice input.  Must be 1,2,3')
end

% phi
HP.eta_tilde_phi            = priorRangeMean; % "Mean" of phi prior
HP.zeta_tilde_phi_2         = priorRangeVar; % "Variance" of phi prior

% lambda (this one's strongly constrained; 95% within 500,2000 km)
HP.eta_tilde_lambda         = priorRangeMean; % "Mean" of phi prior
HP.zeta_tilde_lambda_2      = priorRangeVar; % "Variance" of phi prior


