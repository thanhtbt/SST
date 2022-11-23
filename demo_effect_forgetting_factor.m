%% Demo: Effect of the forgetting factor

clear; clc; close all;

%% Setup
N_exp        = 5;             % number of independent experiments
N_run        = 1;             % number of independent runs / experiment
n            = 50;            % data dimesion
r            = 5;             % rank
T            = 1000;          % number of observations
sparsity     = 0.5;           % sparse density
SNR          = 1e-3;          % noise level 
time_varying = 0e-3*ones(1,T);% time-varying factors
BETA         = [linspace(0.1,0.9,5) 0.95 0.98 1];

%%
eta0 = zeros(length(BETA),T);  rho0 = zeros(length(BETA),T);

%% Processing ...

for jj = 1 : length(BETA)
    beta = BETA(jj)
    fprintf('Forgetting factor =  %0.2f \n',beta)
    eta_OPIT = zeros(1,T);   rho_OPIT =  zeros(1,T);
   
    for n_exp = 1 : N_exp
        fprintf('\n + run %d/%d:\n ',n_exp,N_exp)
        
        %% Data Generation
        [X_stream,U_stream] = online_data_generate(n,T,r,sparsity,time_varying,SNR);
       
        %% Our Method
        eta_OPIT_i = zeros(1,T); 
        rho_OPIT_i = zeros(1,T);
        OPTS_OPIT.method = 'normalization';  %'orthonormalization'; %
        OPTS_OPIT.lambda = beta;
        for ii = 1 : N_run
            t1 = tic;
            [~,eta_ii,rho_ii]  = OPIT(X_stream,OPTS_OPIT,U_stream);
            t1_end = toc(t1);
            eta_OPIT_i   = eta_OPIT_i + eta_ii;
            rho_OPIT_i   = rho_OPIT_i + rho_ii;
        end
        eta_OPIT   = eta_OPIT + eta_OPIT_i/N_run;
        rho_OPIT   = rho_OPIT + rho_OPIT_i/N_run;
               
    end
    
    eta0(jj,:) = eta_OPIT/N_exp;
    rho0(jj,:) = rho_OPIT/N_exp;
   
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

 
fig = figure; hold on; 
semilogy(rho0(1:5,:)','LineWidth',2);
semilogy(rho0(6,1:T)','LineWidth',2.5,'Color','g');
semilogy(1:T,rho0(7,1:T)','LineWidth',2.5,'Color','b');
semilogy(1:T,rho0(8,1:T)','LineWidth',2.5,'Color','r');

xlabel('Time Index - $t$','interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel('$\sin \theta(\mathbf{A}_{t},\mathbf{U}_{t}) $','interpreter','latex','FontSize',13,'FontName','Times New Roman');
lgd = legend('$\beta=0.1$','$\beta=0.3$','$\beta=0.5$','$\beta=0.7$','$\beta=0.9$',...
             '$\beta=0.95$','$\beta=0.98$','$\beta=1$' );
set(lgd,'Interpreter','latex','FontSize',22,'NumColumns',3);

h=gca;
set(gca, 'YScale', 'log')
set(h,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(h,'Xtick',0:200:T,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');
set(h,'FontSize', 30);
axis([0 T 3e-5 1e-1]);
grid on;
box on;
set(fig, 'units', 'inches', 'position', [0.5 0.5 8 7]);
