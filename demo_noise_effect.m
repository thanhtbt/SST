%% Demo: Performance of OPIT in noisy environments


clear; clc; close all;

%% Setup 
N_exp        = 10;              % number of independent experiments
n_run        = 1;               % number of independent runs
n            = 100;             % data dimesion
r            = 5;               % rank
T            = 1000;            % number of observations
sparsity     = 0.9;             % sparse density
SNR_L        = [1e-1 1e-2 1e-3];% noise levels
time_varying = 0*ones(1,T);     % time-varying factors
time_varying(500) = 1;
beta = 0.9;

%% Plot
k  = 5;
TS = 100;
DT = 200;

eta0 = zeros(length(SNR_L),T);  rho0 = zeros(length(SNR_L),T);
eta1 = zeros(length(SNR_L),T);  rho1 = zeros(length(SNR_L),T);
eta2 = zeros(length(SNR_L),T);  rho2 = zeros(length(SNR_L),T);
eta3 = zeros(length(SNR_L),T);  rho3 = zeros(length(SNR_L),T);
eta4 = zeros(length(SNR_L),T);  rho4 = zeros(length(SNR_L),T);
eta5 = zeros(length(SNR_L),T);  rho5 = zeros(length(SNR_L),T);

%% Processing ... 
for jj = 1 : length(SNR_L)
    SNR = SNR_L(jj);
    fprintf('\n Noise level = %0.3f\n', SNR)
    eta_OPIT = zeros(1,T);   rho_OPIT  = zeros(1,T);
    eta_ours = zeros(1,T);   rho_ours =  zeros(1,T);
    
    for N_exp = 1 : N_exp
        fprintf(' + Run %d \n',N_exp)
        
        %% Data Generation     
        [X_stream,U_stream] = online_data_generate(n,T,r,sparsity,time_varying,SNR);
        %% Tracking ...
        eta_OPIT_i = zeros(1,T);
        rho_OPIT_i = zeros(1,T);
        OPTS_OPIT.lambda = beta;
        n_run = 1;
        for ii = 1 : n_run
            t1 = tic;
            [~,eta_ii,rho_ii]  = OPIT(X_stream,OPTS_OPIT,U_stream);
            t1_end = toc(t1);
            eta_OPIT_i   = eta_OPIT_i + eta_ii;
            rho_OPIT_i   = rho_OPIT_i + rho_ii;
        end
        eta_OPIT   = eta_OPIT + eta_OPIT_i/n_run;
        rho_OPIT   = rho_OPIT + rho_OPIT_i/n_run;      
    end
    
    eta0(jj,:) = eta_OPIT/N_exp;
    rho0(jj,:) = rho_OPIT/N_exp;
end

%%


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


% round 1
%%
d11 = semilogy(1:k:T,rho0(1,1:k:end),...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);
d111 = plot(1:TS:T,rho0(1,1:TS:end),...
    'marker','p','markersize',makerSize,...
    'linestyle','none','color',red_o,'LineWidth',LineWidth);
d112 = semilogy(1:1,rho0(1,1:1),...
    'marker','p','markersize',makerSize,...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);

d12 = semilogy(1:k:T,rho0(2,1:k:end),...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);
d121 = plot(1:TS:T,rho0(2,1:TS:end),...
    'marker','d','markersize',makerSize,...
    'linestyle','none','color',gree_o,'LineWidth',LineWidth);
d122 = semilogy(1:1,rho0(2,1:1),...
    'marker','d','markersize',makerSize,...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);

d13 = semilogy(1:k:T,rho0(3,1:k:end),...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);
d131 = plot(1:TS:T,rho0(3,1:TS:end),...
    'marker','s','markersize',makerSize,...
    'linestyle','none','color',blue_o,'LineWidth',LineWidth);
d132 = semilogy(1:1,rho0(3,1:1),...
    'marker','s','markersize',makerSize,...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);

% d14 = semilogy(1:k:T,rho0(4,1:k:end),...
%     'linestyle','-','color',viol_n,'LineWidth',LineWidth);
% d141 = plot(1:TS:T,rho0(4,1:TS:end),...
%  'marker','h','markersize',makerSize,...
%    'linestyle','none','color',viol_n,'LineWidth',LineWidth);
% d142 = semilogy(1:1,rho0(4,1:1),...
%     'marker','h','markersize',makerSize,...
%     'linestyle','-','color',viol_n,'LineWidth',LineWidth);
%


% lgd = legend([d112  d122 d132 d142],'$\sigma_n = 10^{-1}$','$\sigma_n = 10^{-2}$','$\sigma_n = 10^{-3}$','$\sigma_n = 10^{-4}$');
lgd = legend([d112  d122 d132],'$\sigma_n = 10^{-1}$','$\sigma_n = 10^{-2}$','$\sigma_n = 10^{-3}$');
lgd.FontSize = 22;
set(lgd, 'Interpreter', 'latex', 'Color', [0.95, 0.95, 0.95]);


xlabel('Time Index - $t$','interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel('$\sin \theta(\mathbf{A}_{t},\mathbf{U}_{t}) $','interpreter','latex','FontSize',13,'FontName','Times New Roman');

h=gca;
set(gca, 'YScale', 'log')
set(h,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(h,'Xtick',0:DT:T,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');
set(h,'FontSize', 30);
axis([0 T 1e-5 5]);
grid on;
box on;
set(fig, 'units', 'inches', 'position', [0.5 0.5 8 7]);


