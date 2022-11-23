%% Demo: Compare OPIT against FAPI -- the best power-based subspace tracker

clear; clc; close all;
addpath(genpath('Subspace_Trackers'))

%% Setup
N_exp        = 5;              % number of independent experiments
N_run        = 2;              % number of independent runs / experiment
N            = 100;            % data dimesion
r            = 10;             % rank
T            = 2000;           % number of observations
sparsity     = 0.0;            % sparse density
SNR          = 1e-2;           % noise level (dB)
time_varying = 1e-2*ones(1,T); % time-varying factors
time_varying(500)  = 1;        % time-varying factors
time_varying(1500) = 1;        % time-varying factors
beta         = 0.97;           % forgetting factor
%%
eta0_OPIT = zeros(length(N),T);  rho0_OPIT = zeros(length(N),T);
eta0_FAPI = zeros(length(N),T);  rho0_FAPI = zeros(length(N),T);

for jj = 1 : length(N)
    n = N(jj);
   
    eta_OPIT = zeros(1,T);   rho_OPIT =  zeros(1,T);
    eta_FAPI = zeros(1,T);   rho_FAPI =  zeros(1,T);

    for n_exp = 1 : N_exp
        fprintf('\n + run %d/%d: ',n_exp,N_exp)
        
        %% Data Generation
        [X_stream,U_stream] = online_data_generate(n,T,r,sparsity,time_varying,SNR);
               
        eta_FAPI_i = zeros(1,T); 
        rho_FAPI_i = zeros(1,T);
        eta_OPIT_i = zeros(1,T); 
        rho_OPIT_i = zeros(1,T);
        
       
        OPTS_OPIT.method = 'normalization'; 
        OPTS_OPIT.lambda = beta;
        
        for ii = 1 : N_run
            [~,eta_OPIT_ii,rho_OPIT_ii]  = OPIT(X_stream,OPTS_OPIT,U_stream);
            [~,eta_FAPI_ii,rho_FAPI_ii]  = FAPI_Tracking(X_stream,beta,U_stream);
            
            eta_FAPI_i   = eta_FAPI_i + eta_FAPI_ii;
            rho_FAPI_i   = rho_FAPI_i + rho_FAPI_ii;
            
            eta_OPIT_i   = eta_OPIT_i + eta_OPIT_ii;
            rho_OPIT_i   = rho_OPIT_i + rho_OPIT_ii;
        end
        eta_OPIT   = eta_OPIT + eta_OPIT_i/N_run;
        rho_OPIT   = rho_OPIT + rho_OPIT_i/N_run;
         
        eta_FAPI   = eta_FAPI + eta_FAPI_i/N_run;
        rho_FAPI   = rho_FAPI + rho_FAPI_i/N_run;
  
    end
    
    eta0_FAPI(jj,:) = eta_FAPI/N_exp;
    rho0_FAPI(jj,:) = rho_FAPI/N_exp;
    eta0_OPIT(jj,:) = eta_OPIT/N_exp;
    rho0_OPIT(jj,:) = rho_OPIT/N_exp;
end

%% PLOT RESULTS
makerSize   = 14;
numbMarkers = 50;
LineWidth   = 2;
set(0, 'defaultTextInterpreter', 'latex');
color   = get(groot,'DefaultAxesColorOrder');
red_o   = [1,0,0];
blue_o  = [0, 0, 1];
gree_o  = [0, 0.5, 0];
black_o = [0.25, 0.25, 0.25];
blue_n  = color(1,:);
oran_n  = color(2,:);
yell_n  = color(3,:);
viol_n  = color(4,:);
gree_n  = color(5,:);
lblu_n  = color(6,:);
brow_n  = color(7,:);
lbrow_n = [0.5350    0.580    0.2840];

%% 
fig = figure; 

hold on; 
semilogy(rho0_OPIT,'LineWidth',2,'Color','r');
semilogy(rho0_FAPI,'LineWidth',2,'Color','b');
% semilogy(1:T,rho0(7,1:T)','LineWidth',2.5,'Color','b');
% semilogy(1:T,rho0(8,1:T)','LineWidth',2.5,'Color','r');

xlabel('Time Index - $t$','interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel('$\sin \theta(\mathbf{A}_{t},\mathbf{U}_{t}) $','interpreter','latex','FontSize',13,'FontName','Times New Roman');

lgd = legend('OPIT','FAPI');
set(lgd,'Interpreter','latex',...
     'FontSize',22,'NumColumns',3);

 
h=gca;
set(gca, 'YScale', 'log')
set(h,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(h,'Xtick',0:round(T/5):T,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');
set(h,'FontSize', 30);
axis([0 T 5e-5 5e0]);
grid on;
box on;
set(fig, 'units', 'inches', 'position', [0.5 0.5 8 7]);
