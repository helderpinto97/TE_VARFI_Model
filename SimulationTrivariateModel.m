
%% Simulation 1 Trivariate Model -- Complexity and Transfer Entropy

% The initial VAR model:
% 
% $$V_{n}=2 \rho_{v} \cdot \cos 2 \pi f_{v} \cdot V_{n-1}-\rho_{v}^{2}
% \cdot V_{n-2}+U_{v, n}$$
% 
% $$W_{n}=2 \rho_{w} \cdot \cos 2 \pi f_{w} \cdot W_{n-1}-\rho_{w}^{2}
% \cdot W_{n-2}+a \cdot Y_{n-2}+d \cdot V_{n-1}+U_{w, n}$$
% 
% $$Y_{n}=2 \rho_{y} \cdot \cos 2 \pi f_{y} \cdot Y_{n-1}-\rho_{w}^{2}
% \cdot Y_{n-2}+b \cdot W_{n-1}+c \cdot V_{n-1}+U_{y, n}$$
% 
% where $\boldsymbol{U}_{n}=\left[\mathcal{U}_{v, n}, U_{w, n}, U_{y, n}\right]$ is a vector of zero mean white Gaussian noises of unit variance and lated with each other $(\Lambda=I)$.
% 


% Various Configurations of d (long memory parameter)
clear all; close all; clc;

% Parameters
tau_vector = 1:15; % range of time scales to be explored
M=3; % Process Dimension

% Long memory parameters to simulate (Matrix)
%d_matrix=[0 0 0;0.1 0.2 0.7; 0.1 0.3 0.3; 0.1 0.65 0.45];
d_matrix=[0 0 0; 0.1 0.2 0.3];
%d_matrix=[0 0 0];

% Parameters of Coupling (Example only- For the simulation vary the value of c)
a=0.1;
b=0.4;
c=1;
e=c; % e corresponds to d on the model

% Poles of the VAR Model
par.poles{1}=([0.9 0.25]); % Oscillations RESP
par.poles{2}=([0.8 0.1]); % Self-oscillation ABP
par.poles{3}=([0.8 0.1]); % Self-oscillations RR

% Coupling interaction
par.coup=[1 2 1 e; 3 2 2 a; 1 3 1 c; 2 3 1 b];
par.Su=ones(1,M);

% Index of the Target process in the first position.
index_proc=[3 1 2];

% If you wanna save the results -> Uncomment the filename.
filename='Variousdconfigurations.mat';
Results=SimulateVARFImodel(tau_vector,M,d_matrix,par,index_proc,filename);


%% Plots

[f1,f2,f3,l]=plot_results(Results,tau_vector);

%% Simulation 2 Trivariate Model -- Complexity and Transfer Entropy
% Fix long memory parameters dv and dw
% close all;

dv=repelem(0.1,19);
dw=repelem(0.55,19);
dy=linspace(0,0.7,20);
d_matrix=[[0 dv]' [0 dw]' dy'];

% Simulations
filename='Fixdvanddw1.mat';
Results=SimulateVARFImodel(tau_vector,M,d_matrix,par,index_proc,filename);

% Plots
[f5,f6,f7,~]=plot_results(Results,tau_vector);


%% Simulation 3 Trivariate Model -- Complexity and Transfer Entropy
% Fix long memory parameters dv and dy
% close all;

dv=repelem(0.1,19);
dw=linspace(0,0.7,20);
dy=repelem(0.45,19);
d_matrix=[[0 dv]' dw' [0 dy]'];

% Simulations
filename='Fixdvanddy.mat';
Results=SimulateVARFImodel(tau_vector,M,d_matrix,par,index_proc,filename);

% Plots
[f8,f9,f10,~]=plot_results(Results,tau_vector);
