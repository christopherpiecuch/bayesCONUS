%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init_vals_pickup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code used in Piecuch 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and initialize values in case of pickup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by CGP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nn_nonnan = find(~isnan(R),1,'last');

mu=MU(nn_nonnan);
nu=NU(nn_nonnan);
pi_2=PI_2(nn_nonnan);
delta_2=DELTA_2(nn_nonnan);
sigma_2=SIGMA_2(nn_nonnan);
gamma_2=GAMMA_2(nn_nonnan);
tau_2=TAU_2(nn_nonnan);
phi=PHI(nn_nonnan);
lambda=LAMBDA(nn_nonnan);

a=squeeze(A(nn_nonnan,:))';
b=squeeze(B(nn_nonnan,:))';
l=squeeze(ELL(nn_nonnan,:))';
r=squeeze(R(nn_nonnan))';
y_0=squeeze(Y_0(nn_nonnan,:))';

y=squeeze(Y(nn_nonnan,:,:))';

