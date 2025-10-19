%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function bayes_posterior_prediction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by CGP 09 August 2025
% This code performs posterior prediction from tide gauges
% onto regular grid; see Chapter 2.2 in Rasmussen and Williams
% (2006) Gaussian Processes for Machine Learning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bayes_posterior_prediction(experimentName,runnum)

% loop through files
for qqq=1:numel(runnum)
    expNam=(['bayes_model_solutions/experiment_',experimentName,'_runNum_',num2str(runnum(qqq))]);
    
    % load solution
    load([expNam,'.mat'])

    % get coastal cells
    [g_lon,g_lat]=grab_coastal_cells(-90,90,-180,180,TGR_LON,TGR_LAT,[],[],0.5);  

    % do some posterior prediction
    % trends
    Dall=EarthDistances([[LON'; g_lon] [LAT'; g_lat]]);  
    Ng=numel(g_lon);
    ONE=ones(N,1);
    ONEG=ones(Ng,1);
    EYE=eye(N);
    EYEG=eye(Ng);
    BG=zeros(numel(MU),Ng);
    for kk=1:numel(MU), disp(num2str(kk))
        MATMAT=PI_2(kk)*exp(-LAMBDA(kk)*Dall);
        y=B(kk,:)';
        muy=MU(kk)*ONE;
        mux=MU(kk)*ONEG;
        BMAT=MATMAT(1:N,1:N);
        AMAT=MATMAT(((N+1):end),((N+1):end));
        CMAT=MATMAT(((N+1):end),(1:N));
        thing=AMAT-CMAT*inv(BMAT)*CMAT';
        thing=0.5*(thing+thing');
        x=mvnrnd(mux+CMAT*inv(BMAT)*(y-muy),thing);
        BG(kk,:)=x;
    end

    GG=zeros(numel(ALPHA),Ng);
    for kk=1:numel(ALPHA), disp(num2str(kk))
        MATMAT=OMEGA_2(kk)*exp(-RHO(kk)*Dall);
        y=G(kk,:)';
        muy=ALPHA(kk)*ONE;
        mux=ALPHA(kk)*ONEG;
        BMAT=MATMAT(1:N,1:N);
        AMAT=MATMAT(((N+1):end),((N+1):end));
        CMAT=MATMAT(((N+1):end),(1:N));
        thing=AMAT-CMAT*inv(BMAT)*CMAT';
        thing=0.5*(thing+thing');
        x=mvnrnd(mux+CMAT*inv(BMAT)*(y-muy),thing);
        GG(kk,:)=x;
    end

    % innovations
    E=nan*Y;
    T=1900:2024;
    T=T-nanmean(T);
    T0=min(T)-1;

    % initial conditions
    Y_0G=zeros(numel(MU),Ng);
    Y_0G=ones(numel(MU),1)*median(BG)*T0+ones(numel(MU),1)*median(GG)*T0^2+(ones(numel(MU),1)*( abs(std(BG)*T0)+sqrt(median(SIGMA_2)) )).*randn(numel(MU),Ng);
    for kk=1:numel(MU)
        for tt=1:numel(T)
            if tt==1
                E(kk,tt,:)=squeeze(Y(kk,tt,:))'-B(kk,:)*T(tt)-G(kk,:)*T(tt)^2-R(kk)*(squeeze(Y_0(kk,:))-B(kk,:)*T0-G(kk,:)*T0^2);
            else
                E(kk,tt,:)=squeeze(Y(kk,tt,:))'-B(kk,:)*T(tt)-G(kk,:)*T(tt)^2-R(kk)*(squeeze(Y(kk,tt-1,:))'-B(kk,:)*T(tt-1)-G(kk,:)*T(tt-1)^2);
            end
        end
    end
    
    % innovations
    EG=nan*ones(numel(MU),numel(T),numel(g_lon));
    for kk=1:numel(MU), disp(num2str(kk))
        for tt=1:numel(T)
            MATMAT=SIGMA_2(kk)*exp(-PHI(kk)*Dall);
            y=squeeze(E(kk,tt,:));
            muy=0*ONE;
            mux=0*ONEG;
            BMAT=MATMAT(1:N,1:N);
            AMAT=MATMAT(((N+1):end),((N+1):end));
            CMAT=MATMAT(((N+1):end),(1:N));
            thing=AMAT-CMAT*inv(BMAT)*CMAT';
            thing=0.5*(thing+thing');
            x=mvnrnd(CMAT*inv(BMAT)*y,thing);
            EG(kk,tt,:)=x;
        end
    end

    % sea level
    YG=nan*ones(numel(MU),numel(T),numel(g_lon));
    for kk=1:numel(MU), disp(num2str(kk))
        for tt=1:numel(T)
            if tt==1
                YG(kk,tt,:)=(BG(kk,:)*T0)'+(GG(kk,:)*T0^2)'+R(kk)*(Y_0G(kk,:)'-(BG(kk,:)*T0)'-(GG(kk,:)*T0^2)')+squeeze(EG(kk,tt,:));
            else
                YG(kk,tt,:)=(BG(kk,:)*T(tt))'+(GG(kk,:)*T(tt)^2)'+R(kk)*(squeeze(YG(kk,tt-1,:))-(BG(kk,:)*T(tt-1))'-(GG(kk,:)*T(tt-1)^2)')+squeeze(EG(kk,tt,:));
            end
        end
    end

    % format output
    Y=Y;
    YG=YG;
    A=A;
    B=B;
    G=G;
    BG=BG;
    GG=GG;
    LON=LON;
    LAT=LAT;
    GLON=g_lon;
    GLAT=g_lat;
    DATA=TGR_DATA;
    NAME=TGR_NAME;
    YEAR=1900:2024;
    L=ELL;

    save([expNam,'_gridded.mat'],'Y','YG','B','A','BG','LON','LAT','GLON','GLAT','DATA','NAME','YEAR','L','G','GG')
    clearvars -except qqq expNam experimentName runnum
end