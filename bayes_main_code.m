%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function bayes_main_code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by CGP 09 August 2025
% This code represents the RSL submodel of the larger hierarchical 
% Bayesian model described in Piecuch et al. (2018)
% * It only considers tide-gauge data; no GPS data, paleo proxies, or GIA models are included
% * Code below corresponds to their Equations 1-4, 6, and 8 with u=0 (not
%    separating SSH from VLM) and wg=0 (not distinguishing GIA from
%    non-GIA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bayes_main_code(experimentName,runnum)

% make sure you have a save directory
if exist('bayes_model_solutions/')==0
    mkdir bayes_model_solutions/
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Say hello
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause(0.1)
disp('Hello.  Things have started.')
pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some preliminary input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
% name of experiment and file to be saved
save_name=[experimentName,'_runNum_',num2str(runnum)];  % file name

%%%
% number of draws to perform
% for analysis
NN_burn=200000; % 200,000 warm-up draws
NN_post=200000; % 200,000 post-warm-up draws
thin_period=200; % thin chains keeping 1 of 200
% for testing
%NN_burn=100; 
%NN_post=100; 
%thin_period=1;
startFromPickup=0; % start from existing pickup or not; default=0

%%%
% geographic bounds
% whole earth
inp.la1=-90; % Southern latitudinal bounds of study region
inp.la2=90; % Northern latitudinal bounds "
inp.lo1=-180; % Western longitudinal bounds "
inp.lo2=180; % Eastern longitudinal bounds "

%%%
% minimum length of tide gauge record to be involved in the estimation
inp.minnum=30;  % units are years

%%%
% PSMSL coastline codes to be considered (see http://www.psmsl.org/data/obtaining/)
inp.coastcode=[940 960 823]; % PSMSL ID for CONUS

%%%
% start year of inference
inp.t0 = 1900;

%%%
% prior on inverse range parameters
inp.priorChoice=1;			% prior range of ~500,2000 km (default)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define iteration parameters based on input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NN_burn_thin=NN_burn/thin_period; % Total number of burn-in to keep
NN_post_thin=NN_post/thin_period; % Total number of post-burn-in to keep
NN=NN_burn+NN_post;  % Total number of draws to take 
NN_thin=NN_burn_thin+NN_post_thin;% Total number of draws to keep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare and load PSMSL annual tide gauge RSL records
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
% load tide gauge data
[TGR_DATA,TGR_LON,TGR_LAT,TGR_NAME,TGR_COAST] = prepare_tgr_data(inp.la1,...
    inp.la2,inp.lo1,inp.lo2,inp.minnum,inp.coastcode,inp.t0);

%%%
% include 0.5-degree regular grid cells within study domain to estimate
% the solution at
[G_LON,G_LAT] = grab_coastal_cells(inp.la1,inp.la2,inp.lo1,inp.lo2,...
    TGR_LON,TGR_LAT,[],[],0.5);  
g_lon=G_LON; g_lat=G_LAT;
G_LON=[]; G_LAT=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define space and time parameters related to data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up dimensions of process
[N_tgr,K_tgr]=size(TGR_DATA);
[N_gia,K_gia]=size(G_LON);

% bring together all target locations
N = N_tgr+N_gia;
K = K_tgr;
LON = [TGR_LON G_LON'];
LAT = [TGR_LAT G_LAT' ];

% define some useful space and time values
D=EarthDistances([LON' LAT']); % distances between locations
T=1:K; T=T-mean(T); % time steps (centered on zero)
T0=T(1)-1; % initial time step (relative to centered series)
M=sum(sum(~isnan(TGR_DATA'))~=0); % number of locations with at least one datum
% define names of all target locations
for n=1:N
    if n<=N_tgr
        NAME(n).name = TGR_NAME(n).name;
    elseif n>N_tgr
        NAME(n).name = [num2str(abs(G_LON(n-N_tgr))),'W_',...
            num2str(abs(G_LAT(n-N_tgr))),'N'];
    end
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cluster locations east and west of 100W
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = determine_clusters(LON);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the seeds of the random number generator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rng(runnum*sum(clock))
rng(runnum)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate space for the sample arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha=[];
initialize_output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the hyperparameter values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HP = set_hyperparameters(N,K,M,LON,LAT,TGR_DATA,inp.priorChoice);

%%%%%%%%%%%%%%%%%%%%
% Set initial values
%%%%%%%%%%%%%%%%%%%%
set_initial_values

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up selection matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%
H_master=double(~isnan(TGR_DATA));
M_k=sum(H_master);
for k=1:K
    gauges_with_data(k).indices=find(H_master(:,k)~=0);
    selMat(k).H=zeros(M_k(k),N);
    selMat(k).F=zeros(M_k(k),M);
    for m_k=1:M_k(k)
       selMat(k).H(m_k,gauges_with_data(k).indices(m_k))=1;
       selMat(k).F(m_k,gauges_with_data(k).indices(m_k))=1;
    end
    Z(k).z=squeeze(TGR_DATA(gauges_with_data(k).indices,k));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up identity matrices and vectors of zeros or ones
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_N=eye(N);
I_M=eye(M);
ONE_N=ones(N,1);
ONE_M=ones(M,1);
ZERO_N=zeros(N,1);
ZERO_M=zeros(M,1);

for k=1:K
   I_MK(k).I=eye(M_k(k));
   ONE_MK(k).ONE=ones(M_k(k),1);
   ZERO_MK(k).ZERO=zeros(M_k(k),1);
end
III=inv(HP.zeta_tilde_y_0_2);

nn_start=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Bayes code takes a while to run.  So, the code periodically saves out
% the solution "so far".  This code immediately below is to, in the case
% that the model crashes, start from a pickup file, so as to prevent having
% to start from scratch.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if startFromPickup
    load(['bayes_model_solutions/experiment_',save_name,'.mat'],...
            'MU','NU','PI_2','DELTA_2','SIGMA_2','TAU_2','GAMMA_2',...
            'PHI','LAMBDA','A','B','ELL','R','Y_0','Y', 'HP','*DATA','*LON','*LAT',...
            '*NAME','N','K','M','D','nn','inp','g_lon','g_lat','G','OMEGA_2','RHO','ALPHA')
    % set start time
    nn_start=nn+1; clear nn
    % reset initial values
    init_vals_pickup
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through the Gibbs sampler with Metropolis steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
disp('starting.')
for nn=nn_start:NN, 
    if mod(nn,10)==0, 
        toc, 
        disp([num2str(nn),' of ',num2str(NN),' iterations done.']), 
        tic, 
    end
    nn_thin=[]; nn_thin=ceil(nn/thin_period);
    
    if (mod(nn,(NN_post/10))==0)
        % periodically output, just in case of stall or crash
        disp('hello') 
        save(['bayes_model_solutions/experiment_',save_name,'.mat'],...
            'MU','NU','PI_2','DELTA_2','SIGMA_2','TAU_2','GAMMA_2',...
            'PHI','LAMBDA','A','B','ELL','R','Y_0','Y', 'HP','*DATA','*LON','*LAT',...
            '*NAME','N','K','M','D','nn','inp','G','ALPHA','RHO','OMEGA_2','G','ALPHA','RHO','OMEGA_2')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define matrices to save time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SigmaMat=sigma_2*(C.*exp(-phi*D)); invSigmaMat=inv(SigmaMat);
    SMat=C.*exp(-phi*D); invSMat=inv(SMat);
    PiMat=pi_2*exp(-lambda*D); invPiMat=inv(PiMat); 
    LMat=exp(-lambda*D); invLMat=inv(LMat);
    OmegaMat=omega_2*exp(-rho*D); invOmegaMat=inv(OmegaMat); 
    RMat=exp(-rho*D); invRMat=inv(RMat);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(y_K|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V_Y_K=[]; PSI_Y_K=[];
    V_Y_K=delta_2^(-1)*(selMat(K).H'*(Z(K).z-selMat(K).F*(l+a*T(K))))+...
    	invSigmaMat*(r*y(:,K-1)+(T(K)-r*T(K-1))*b+(T(K)^2-r*T(K-1)^2)*g);
    PSI_Y_K=inv(1/delta_2*selMat(K).H'*selMat(K).H+invSigmaMat);
    y(:,K)=mvnrnd(PSI_Y_K*V_Y_K,PSI_Y_K)';
    clear V_Y_K PSI_Y_K   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(y_k|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for kk=(K-1):-1:1
    	V_Y_k=[]; PSI_Y_k=[];
      	if kk==1
        	V_Y_k=1/delta_2*(selMat(1).H'*(Z(1).z-selMat(1).F*(l+a*T(1))))+...
                invSigmaMat*(r*(y_0+y(:,2))+(1+r^2)*T(1)*b-r*(T0+T(2))*b+(1+r^2)*T(1)^2*g-r*(T0^2+T(2)^2)*g);
        else
         	V_Y_k=1/delta_2*(selMat(kk).H'*(Z(kk).z-selMat(kk).F*(l+a*T(kk))))+...
            	invSigmaMat*(r*(y(:,kk-1)+y(:,kk+1))+(1+r^2)*T(kk)*b-r*(T(kk-1)+T(kk+1))*b+(1+r^2)*T(kk)^2*g-r*(T(kk-1)^2+T(kk+1)^2)*g);
        end
       	PSI_Y_k=inv(1/delta_2*selMat(kk).H'*selMat(kk).H+(1+r^2)*invSigmaMat);
       	y(:,kk)=mvnrnd(PSI_Y_k*V_Y_k,PSI_Y_k)';
      	clear V_Y_k PSI_Y_k 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(y_0|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V_Y_0=[]; PSI_Y_0=[];
    V_Y_0=III*HP.eta_tilde_y_0+invSigmaMat*(r*y(:,1)-r*(T(1)-r*T0)*b-r*(T(1)^2-r*T0^2)*g);
    PSI_Y_0=inv(III+r^2*invSigmaMat);
    y_0=mvnrnd(PSI_Y_0*V_Y_0,PSI_Y_0)';
    clear V_Y_0 PSI_Y_0
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(a|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  	V_A=[]; PSI_A=[]; SUM_K1=ZERO_M; SUM_K2=zeros(M);
    for kk=1:K
    	SUM_K1=SUM_K1+T(kk)*selMat(kk).F'*(Z(kk).z-...
         	selMat(kk).H*y(:,kk)-selMat(kk).F*l);
        SUM_K2=SUM_K2+T(kk)^2*(selMat(kk).F'*selMat(kk).F);
   	end
    V_A=(1/delta_2)*SUM_K1;
    PSI_A=inv(1/gamma_2*eye(M)+1/delta_2*SUM_K2);
    a=mvnrnd(PSI_A*V_A,PSI_A)';
    clear V_A PSI_A SUM_K*

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(delta_2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUM_K=0;
    for kk=1:K
      	xxx=[]; xxx=(Z(kk).z-selMat(kk).H*y(:,kk)-selMat(kk).F*(l+a*T(kk)));
      	SUM_K=SUM_K+xxx'*xxx;
   	end
    delta_2=1/randraw('gamma', [0,1/(HP.chi_tilde_delta_2+1/2*SUM_K),...
     	(HP.xi_tilde_delta_2+1/2*sum(M_k))], [1,1]);    
    clear SUM_K
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(r|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V_R=0; PSI_R=0;
    for kk=1:K
     	if kk==1
         	V_R=V_R+((y_0-b*T0-g*T0^2)')*invSigmaMat*(y(:,1)-b*T(1)-g*T(1)^2);
         	PSI_R=PSI_R+((y_0-b*T0-g*T0^2)')*invSigmaMat*(y_0-b*T0-g*T0^2);
        else
         	V_R=V_R+((y(:,kk-1)-b*T(kk-1)-g*T(kk-1)^2)')*invSigmaMat*(y(:,kk)-b*T(kk)-g*T(kk)^2);
          	PSI_R=PSI_R+((y(:,kk-1)-b*T(kk-1)-g*T(kk-1)^2)')*invSigmaMat*(y(:,kk-1)-b*T(kk-1)-g*T(kk-1)^2);
        end        
   	end
    PSI_R=inv(PSI_R);
    sample=normrnd(PSI_R*V_R,sqrt(PSI_R));
    if sample>HP.v_tilde_r;
        sample=HP.v_tilde_r;
    end
    if sample<HP.u_tilde_r
        sample=HP.u_tilde_r;
    end   
    r=sample;
    clear V_R PSI_R sample

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(sigma_2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUM_K=0;
    for kk=1:K
     	if kk==1
         	DYKK=[];
          	DYKK=y(:,1)-r*y_0-(T(1)-r*T0)*b-(T(1)^2-r*T0^2)*g;
         	SUM_K=SUM_K+(DYKK')*invSMat*DYKK;           
        else
         	DYKK=[];
           	DYKK=y(:,kk)-r*y(:,kk-1)-(T(kk)-r*T(kk-1))*b-(T(kk)^2-r*T(kk-1)^2)*g;
         	SUM_K=SUM_K+(DYKK')*invSMat*DYKK;           
        end
    end
   	sigma_2=1/randraw('gamma', [0,1/(HP.chi_tilde_sigma_2+1/2*SUM_K),...
     	(HP.xi_tilde_sigma_2+N*K/2)], [1,1]);
   	clear SUM_K DYKK
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(phi|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Phi_now=log(phi);
    Phi_std=0.05;
    Phi_prp=normrnd(Phi_now,Phi_std);
    R_now=C.*exp(-exp(Phi_now)*D);
    R_prp=C.*exp(-exp(Phi_prp)*D);
    invR_now=inv(R_now);
    invR_prp=inv(R_prp);
    sumk_now=0;
    sumk_prp=0;
   	for kk=1:K
      	if kk==1
         	DYYK=y(:,1)-r*y_0-(T(1)-r*T0)*b-(T(1)^2-r*T0^2)*g;
          	sumk_now=sumk_now+(DYYK')*invR_now*DYYK;
          	sumk_prp=sumk_prp+(DYYK')*invR_prp*DYYK;
        else
         	DYYK=y(:,kk)-r*y(:,kk-1)-(T(kk)-r*T(kk-1))*b-(T(kk)^2-r*T(kk-1)^2)*g;
         	sumk_now=sumk_now+(DYYK')*invR_now*DYYK;
         	sumk_prp=sumk_prp+(DYYK')*invR_prp*DYYK;
        end
    end
        
 	ins_now=-1/(2*HP.zeta_tilde_phi_2)*(Phi_now-HP.eta_tilde_phi)^2-1/(2*sigma_2)*sumk_now;
   	ins_prp=-1/(2*HP.zeta_tilde_phi_2)*(Phi_prp-HP.eta_tilde_phi)^2-1/(2*sigma_2)*sumk_prp;
  	MetFrac=det(R_prp*invR_now)^(-K/2)*exp(ins_prp-ins_now);
   	success_rate=min(1,MetFrac);
   	if rand(1)<=success_rate
     	Phi_now=Phi_prp; 
    end
  	phi=exp(Phi_now);
  	clear Phi_now Phi_std Phi_prp mat_now mat_prp ins_* sumk MetFrac success_rate R_*
    % redefine matrices since you just updated sigma_2 and phi
    SigmaMat=sigma_2*(C.*exp(-phi*D)); invSigmaMat=inv(SigmaMat);
    SMat=C.*exp(-phi*D); invSMat=inv(SMat);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(l|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V_L=[]; PSI_L=[]; SUM_K1=ZERO_M; SUM_K2=zeros(M,M);
    for kk=1:K
     	SUM_K1=SUM_K1+(selMat(kk).F')*(Z(kk).z-selMat(kk).H*y(:,kk)-...
            selMat(kk).F*a*T(kk));
        SUM_K2=SUM_K2+(selMat(kk).F')*selMat(kk).F;
    end
    V_L=nu/tau_2*ONE_M+1/delta_2*SUM_K1;
    PSI_L=inv(1/tau_2*I_M+1/delta_2*SUM_K2);
    l=mvnrnd(PSI_L*V_L,PSI_L)';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(nu|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V_NU=[]; PSI_NU=[];
   	V_NU=HP.eta_tilde_nu/HP.zeta_tilde_nu_2+1/tau_2*(ONE_M'*l);
    PSI_NU=inv(1/HP.zeta_tilde_nu_2+M/tau_2);
    nu=normrnd(PSI_NU*V_NU,sqrt(PSI_NU));
    clear V_NU PSI_NU
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(tau_2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tau_2=1/randraw('gamma', [0,1/(HP.chi_tilde_tau_2+...
     	1/2*(((l-nu*ONE_M)')*(l-nu*ONE_M))),(HP.xi_tilde_tau_2+M/2)], [1,1]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(gamma_2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gamma_2=1/randraw('gamma', [0,1/(HP.chi_tilde_gamma_2+...
     	1/2*a'*a),(HP.xi_tilde_gamma_2+M/2)], [1,1]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(mu|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V_MU=[]; PSI_MU=[];
    V_MU=HP.eta_tilde_mu/HP.zeta_tilde_mu_2+ONE_N'*invPiMat*(b);
   	PSI_MU=inv(1/HP.zeta_tilde_mu_2+ONE_N'*invPiMat*ONE_N);
    mu=normrnd(PSI_MU*V_MU,sqrt(PSI_MU));
    clear V_MU PSI_MU  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(alpha|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V_ALPHA=[]; PSI_ALPHA=[];
    V_ALPHA=HP.eta_tilde_alpha/HP.zeta_tilde_alpha_2+ONE_N'*invOmegaMat*(g);
    PSI_ALPHA=inv(1/HP.zeta_tilde_alpha_2+ONE_N'*invOmegaMat*ONE_N);
    alpha=normrnd(PSI_ALPHA*V_ALPHA,sqrt(PSI_ALPHA));
    clear V_ALPHA PSI_ALPHA  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(pi_2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    inside1=[]; inside2=[];
    inside1=N/2;
    inside2=1/2*((b-mu*ONE_N)'*invLMat)*(b-mu*ONE_N);
    pi_2=1/randraw('gamma', [0,1/(HP.chi_tilde_pi_2+inside2),...
      	(HP.xi_tilde_pi_2+inside1)], [1,1]);
   	clear inside*
    PiMat=pi_2*exp(-lambda*D); invPiMat=inv(PiMat); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(omega_2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    inside1=[]; inside2=[];
    inside1=N/2;
    inside2=1/2*((g-alpha*ONE_N)'*invRMat)*(g-alpha*ONE_N);
    omega_2=1/randraw('gamma', [0,1/(HP.chi_tilde_omega_2+inside2),...
      	(HP.xi_tilde_omega_2+inside1)], [1,1]);
   	clear inside*
    OmegaMat=omega_2*exp(-rho*D); invOmegaMat=inv(OmegaMat); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(lambda|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Lambda_now=log(lambda);
    Lambda_std=0.05;
    Lambda_prp=normrnd(Lambda_now,Lambda_std);
    R_now=exp(-exp(Lambda_now)*D);
    R_prp=exp(-exp(Lambda_prp)*D);
    invR_now=inv(R_now);
    invR_prp=inv(R_prp);
    sumk_now=0;
    sumk_prp=0;
    sumk_now=(b-mu*ONE_N)'*invR_now*(b-mu*ONE_N);
    sumk_prp=(b-mu*ONE_N)'*invR_prp*(b-mu*ONE_N);
 	ins_now=-1/(2*HP.zeta_tilde_lambda_2)*(Lambda_now-HP.eta_tilde_lambda)^2-1/(2*pi_2)*sumk_now;
   	ins_prp=-1/(2*HP.zeta_tilde_lambda_2)*(Lambda_prp-HP.eta_tilde_lambda)^2-1/(2*pi_2)*sumk_prp;
  	MetFrac=det(R_prp*invR_now)^(-1/2)*exp(ins_prp-ins_now);
   	success_rate=min(1,MetFrac);
   	if rand(1)<=success_rate
     	Lambda_now=Lambda_prp; 
    end
  	lambda=exp(Lambda_now);
  	clear Lambda_now Lambda_std Lambda_prp mat_now mat_prp ins_* sumk MetFrac success_rate R_*
    % update matrices since you just updated pi_2 and lambda
    PiMat=pi_2*exp(-lambda*D); invPiMat=inv(PiMat); 
    LMat=exp(-lambda*D); invLMat=inv(LMat);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(rho|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Rho_now=log(rho);
    Rho_std=0.05;
    Rho_prp=normrnd(Rho_now,Rho_std);
    R_now=exp(-exp(Rho_now)*D);
    R_prp=exp(-exp(Rho_prp)*D);
    invR_now=inv(R_now);
    invR_prp=inv(R_prp);
    sumk_now=0;
    sumk_prp=0;
    sumk_now=(g-alpha*ONE_N)'*invR_now*(g-alpha*ONE_N);
    sumk_prp=(g-alpha*ONE_N)'*invR_prp*(g-alpha*ONE_N);
 	ins_now=-1/(2*HP.zeta_tilde_rho_2)*(Rho_now-HP.eta_tilde_rho)^2-1/(2*omega_2)*sumk_now;
   	ins_prp=-1/(2*HP.zeta_tilde_rho_2)*(Rho_prp-HP.eta_tilde_rho)^2-1/(2*omega_2)*sumk_prp;
  	MetFrac=det(R_prp*invR_now)^(-1/2)*exp(ins_prp-ins_now);
   	success_rate=min(1,MetFrac);
   	if rand(1)<=success_rate
     	Rho_now=Rho_prp; 
    end
  	rho=exp(Rho_now);
  	clear Rho_now Rho_std Rho_prp mat_now mat_prp ins_* sumk MetFrac success_rate R_*
    % update matrices since you just updated rho
    OmegaMat=omega_2*exp(-rho*D); invOmegaMat=inv(OmegaMat); 
    RMat=exp(-rho*D); invRMat=inv(RMat);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(b|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   	V_B=[]; PSI_B=[]; SUM_K=ZERO_N;
    for kk=1:K
     	if kk==1
         	SUM_K=SUM_K+(T(1)-r*T0)*(y(:,1)-r*y_0-g*(T(1)^2-r*T0^2));
        else
          	SUM_K=SUM_K+(T(kk)-r*T(kk-1))*(y(:,kk)-r*y(:,kk-1)-g*(T(kk)^2-r*T(kk-1)^2));
        end
   	end
    V_B=invPiMat*(mu*ONE_N)+invSigmaMat*SUM_K;
    PSI_B=inv(invPiMat+invSigmaMat*sum((T-r*[T0 T(1:K-1)]).^2));
    b=mvnrnd(PSI_B*V_B,PSI_B)';
    clear V_B PSI_B SUM_K

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(g|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   	V_G=[]; PSI_G=[]; SUM_K=ZERO_N;
    for kk=1:K
     	if kk==1
         	SUM_K=SUM_K+(T(1)^2-r*T0^2)*(y(:,1)-r*y_0-b*(T(1)-r*T0));
        else
          	SUM_K=SUM_K+(T(kk)^2-r*T(kk-1)^2)*(y(:,kk)-r*y(:,kk-1)-b*(T(kk)-r*T(kk-1)));
        end
   	end
    V_G=invOmegaMat*(alpha*ONE_N)+invSigmaMat*SUM_K;
    PSI_G=inv(invOmegaMat+invSigmaMat*sum((T.^2-r*[T0 T(1:K-1)].^2).^2));
    g=mvnrnd(PSI_G*V_G,PSI_G)';
    clear V_G PSI_G SUM_K

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now update arrays
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    update_all_arrays

end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delete the burn-in period values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delete_burn_in

%%%%%%%%%%%%%
% save output
%%%%%%%%%%%%%
save(['bayes_model_solutions/experiment_',save_name,'.mat'],...
    'MU','NU','PI_2','DELTA_2','SIGMA_2','TAU_2','GAMMA_2',...
    'PHI','LAMBDA','A','B','ELL','R','Y_0','Y', 'HP','*DATA','*LON','*LAT',...
    '*NAME','N','K','M','D','nn','inp','g_lon','g_lat','G','ALPHA','RHO','OMEGA_2')
