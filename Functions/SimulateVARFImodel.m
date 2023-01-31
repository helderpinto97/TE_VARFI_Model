% Simulation of VARFI model for various configurations of d (Using a VAR model of superior order) and
% calculation of MSE and TE (Cardiovascular System)
% 
% Input:
% 
% tau_vector - vector of the time scales
% 
% M - Dimension of the process to simulate
% 
% d_matrix - Configurations of long memory vector (in rows)
% 
% par - Structure with the model parameters (poles, coupling interactions etc)
% 
% par.poles{1}=([]);
% par.poles{2}=([]); ...
% par.Su=ones(1,M);
% par.coup=[] matrix of the coupling interaction 
% 
% index_proc - vector with the index of the targeta and the source
% processes (dimension 3 and target in first position)
% 
% filename - name for results (.mat file)

function [Results]=SimulateVARFImodel(tau_vector,M,d_matrix,par,index_proc,filename)

% Error

if length(par.poles)~=M
    error(['The length of poles cell must be equal to ' int2str(M)]);

elseif length(par.Su)~=M
    error(['The length of Su_diag must be ' int2str(M)]);
end

if sum(d_matrix<0)>0
    error('The d values must be positive');
end

if length(index_proc)>=4 || length(index_proc)<2
    error('The index_proc must be equal 2, for bivariate models (target and source) or 3 in multivariate models (one target and 2 two sources)');
end

if sum(M==[2 3])==0
    error('Function only works for bivariate and trivariate models');
end

q_trunc = 50;    % lag at which truncate (Filtro d)
ncoeff_FIR=48; % number of coeffs of FIR lowpass filter for linear MSE

nscales=length(tau_vector); % number of times scales to use
nsim=size(d_matrix,1);% Number simulations (n rows of d_matrix)

% Theorical VAR Parameters
[Am,Su,~,~]=theoreticalVAR(M,par);

% C polynomial
Cpol=[eye(M) -Am'];

% Pre allocate the cell for the results
Complexity=cell(nscales,nsim);

for j=1:nsim
    clc;disp(['d configuration' ' ' int2str(j)]);
    [Am_ARFI,~] = arfi2ar2(Cpol,d_matrix(j,:),q_trunc);
    
    for is=1:nscales
    disp(['scale ' int2str(is) ' of ' int2str(nscales) '...'])
    tau=round(tau_vector(is));

    % LINEAR MSE on high-order VAR (computes both information storage IS and complexity IC)
    % core function: FIR filter for rescaling, VAR to VARMA to SS models,
    % downsampled SS process, estimation of partial variances as in [Faes et al Entropy 2017, ref. 19 of main paper]
    % and finally computation of global complexity (Eq. 3) and individual complexity (Eqs. 4a,b,c) for all target series
    [~,IC]=ar2mse2(Am_ARFI,Su,tau,ncoeff_FIR);
    Complexity{is,j}=IC;
    end
end

if M==3
   Results=Multivariate_aux_fun(Complexity,nsim,nscales,index_proc);
else
   Results=Bivariate_aux_fun(Complexity,nsim,nscales,index_proc);
end

% Save parameters of model in Structure
Results.par=par;

% Save d configurations on the Structure
Results.d_matrix=d_matrix;

% Save the Results
if isempty(filename)~=1
    save(filename,'-struct','Results');
end

end
