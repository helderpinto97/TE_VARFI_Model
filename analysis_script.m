%% EXEMPLARY SCRIPT FOR COMPUTING MULTISCALE ENTROPY (CE_x) BASED ON MULTIVARIATE ARFI MODELS
clear; close all; clc;

%% 1) open data and set analysis parameters
load('example_multivariate_series.mat')
Yo=[series_HP; series_SAP; series_RESP]; % all series in a matrix (series in row)

tau_vect = 1:12; % range of time scales to be explored
q_trunc = 50;    % lag at which truncate (Filtro d)
pmax=20; % maximum  order for model order selection
ncoeff=48; % number of coeffs of FIR lowpass filter for linear MSE
d_min=-0.5; d_max=1; %range for d_estimation()

f_taus = 1./(2*tau_vect); %cutoff frequencies of rescaling filters
nscales=length(f_taus);  %number of time scales spanned

jj=1; jj_label='HP'; % index of target
ii=2; ii_label='SAP'; % index of first driver
kk=3; kk_label='RESP'; % index of second driver


%% 2) estimate differencing parameter and filter each series to remove long-range correlations

% Number of processes
M=size(Yo,1);
for m=1:M
    %d estimation (Whittle semiparametric estimator)
    est_d=d_estimation(Yo(m,:),d_min,d_max);
    ed(m)=est_d;
    %d filter
    [Ytmp,~]=remove_d(Yo(m,:)',est_d); %filtered data
    Y(m,:)=Ytmp';
end

% model identification
[~,p_ARd,~,~] = mos_idMVAR(Y,pmax,0); % model order selection (Bayesian Information criterion)
[eAm_ARd,eSu_ARd,~,~]=idMVAR(Y,p_ARd,0); % identification (ordinary least squares)


%% 3) Approximation of the VARFI process with a finite-order VAR process
eCpol=[eye(M) -eAm_ARd]; 
[eAm_ARFI,eApol] = arfi2ar2(eCpol,ed,q_trunc); 

%% 4) Computation of TE Decompositions for VAR process
for s=1:nscales
    disp(['scale ' int2str(s)]);
    tau=tau_vect(s);

    % MA parameters resulting from the change of scale
    if tau==1
        q=0; b=1;
    else
        q=ncoeff; % number of filter coeffs
        ft=1/(2*tau); %cutoff frequency
        Wn=2*ft; %normalized cutoff frequency (fNyquist=1)
        b=fir1(q,Wn,'noscale'); %Hamming window, linear phase (symmetry of b coeffs)
    end
    Bk=zeros(M,M,q+1);
    for l=1:q+1
        Bk(:,:,l)=b(l)*eye(M);
    end
    Bm=[];
    for kp=1:q+1
        Bm=[Bm Bk(:,:,kp)];
    end
    B0=Bm(1:M,1:M);
    Bm=Bm(1:M,M+1:end);

    
    % ISS parameters
    [A,C,K,V,Vy] = varma2iss(eAm_ARFI,Bm,eSu_ARd,B0); % max(abs(eig(A-K*C)))
    
    %%% parameters after downsampling
    [Ad,Kd,Vd] = iss_ds(A,C,K,V,tau);
    Cd=C; %Rtau=R;
    
    [VR, lambda0] = iss_PCOV(Ad,Cd,Kd,Vd,jj);
    Sj=lambda0(jj,jj);
    Sj_j=VR;
    
    tmp = iss_PCOV(Ad,Cd,Kd,Vd,[jj ii]);
    Sj_ji=tmp(1,1);
    
    tmp = iss_PCOV(Ad,Cd,Kd,Vd,[jj kk]);
    Sj_jl=tmp(1,1);
    
    tmp = iss_PCOV(Ad,Cd,Kd,Vd,[jj ii kk]);
    Sj_ijl=tmp(1,1);
       
    % Interaction Information Decomposition
    Til_j(s)=0.5*log(Sj_j/Sj_ijl); %Joint transfer
    Ti_j(s)=0.5*log(Sj_j/Sj_ji); % Individual transfer
    Tl_j(s)=0.5*log(Sj_j/Sj_jl);
    Ij_il(s)=-Ti_j(s)-Tl_j(s)+Til_j(s); %Interaction transfer (NET SYNERGY)
    
    % Partial Information Decomposition
    Ril_j(s)=min(Ti_j(s),Tl_j(s)); % Redundant transfer
    Ui_j(s)=Ti_j(s)-Ril_j(s); % Unique transfer
    Ul_j(s)=Tl_j(s)-Ril_j(s); % Unique transfer
    Sil_j(s)=Til_j(s)-Ui_j(s)-Ul_j(s)-Ril_j(s); %Synergistic transfer
    
end

%% Information Measures Profiles 

% Plot the time series
serie = {'HP','SBP','RESP'};  %labels
for m = 1:size(Yo,1)
    figure(1); subplot(3,1,m); 
    plot(Yo(m,:)); title([serie{m}]); xlim([1 size(Yo,2)]);
    figure(2); subplot(3,1,m); 
    plot(Y(m,:),'r'); title([serie{m} ', after removing d and normalization']); xlim([1 size(Yo,2)]);
end

%% IID Measures Plots

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1)
plot(Til_j,'-ok'); hold on;
plot(Ti_j,'-ob');
plot(Tl_j,'-or'); 
plot(Ij_il,'-og'); hold off;
xlabel('\tau');
ylabel('TEs');
legend('T_{i,l \rightarrow j}','T_{i \rightarrow j}','T_{l \rightarrow j}','I_{i,j \rightarrow j}')
title('IID Measures');


subplot(1,2,2)
plot(Ui_j,'-ok'); hold on;
plot(Ul_j,'-ob');
plot(Ril_j,'-or');
plot(Sil_j,'-og'); hold off;
legend('U_{i \rightarrow j}','U_{l \rightarrow j}','R_{i,l \rightarrow j}','S_{i,l \rightarrow j}')
title('PID Measures');
xlabel('\tau');
ylabel('TEs');
sgtitle(['MSE of series j=' jj_label '(i=' ii_label ',k=' kk_label ')']);
