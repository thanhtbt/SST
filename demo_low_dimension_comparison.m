%% Demo: Performance Comparison between OPIT and state-of-the-art algorithms
%%       in the classical regime (n is small)

clear; clc; close all;
addpath(genpath('Subspace_Trackers/'))

%% Setup
N_exp        = 3;             % number of independent experiments
N_run        = 1;             % number of independent runs
n            = 100;           % data dimesion
r            = 10;            % rank
T            = 1000;          % number of observations
sparsity     = 0.9;           % sparse density
SNR_L        = 1e-1;          % noise level (dB)
time_varying = 1e-3*ones(1,T);% time-varying factors
beta         = 0.97;

%% Evaluation metrics
eta0 = zeros(length(SNR_L),T);  rho0 = zeros(length(SNR_L),T);
eta1 = zeros(length(SNR_L),T);  rho1 = zeros(length(SNR_L),T);
eta2 = zeros(length(SNR_L),T);  rho2 = zeros(length(SNR_L),T);
eta3 = zeros(length(SNR_L),T);  rho3 = zeros(length(SNR_L),T);
eta4 = zeros(length(SNR_L),T);  rho4 = zeros(length(SNR_L),T);
eta5 = zeros(length(SNR_L),T);  rho5 = zeros(length(SNR_L),T);


%% Processing ... 
for snr_i = 1 : length(SNR_L)
    SNR = SNR_L(snr_i);
    
    eta_OPIT  = zeros(1,T);   rho_OPIT  =  zeros(1,T);
    eta_L1    = zeros(1,T);   rho_L1    = zeros(1,T);
    eta_GSS   = zeros(1,T);   rho_GSS   = zeros(1,T);
    eta_SS    = zeros(1,T);   rho_SS    = zeros(1,T);
    eta_PAST  = zeros(1,T);   rho_PAST  = zeros(1,T);
    eta_SSPCA = zeros(1,T);   rho_SSPCA = zeros(1,T);
    eta_Oja   = zeros(1,T);   rho_Oja   = zeros(1,T);
    eta_OPAST = zeros(1,T);   rho_OPAST = zeros(1,T);
    eta_FAPI  = zeros(1,T);   rho_FAPI  = zeros(1,T);
    eta_GYAST = zeros(1,T);   rho_GYAST = zeros(1,T);
    eta_LORAF = zeros(1,T);   rho_LORAF = zeros(1,T);

    for n_exp = 1 : N_exp
        fprintf('\n Run %d/%d:\n \n ',n_exp,N_exp)
        
        %% Data Generation
        [X_stream,U_stream] = online_data_generate(n,T,r,sparsity,time_varying,SNR);
       
        %% The Proposed  Method
        fprintf('+ OPIT (Proposed): \n');
        eta_OPIT_i = zeros(1,T); 
        rho_OPIT_i = zeros(1,T);
        OPTS_OPIT.lambda = beta;
        OPTS_OPIT.window = 1; 
        for ii = 1 : N_run
            t1 = tic;
            [~,eta_ii,rho_ii]  = OPIT(X_stream,OPTS_OPIT,U_stream);
            t1_end = toc(t1);
            eta_OPIT_i   = eta_OPIT_i + eta_ii;
            rho_OPIT_i   = rho_OPIT_i + rho_ii;
        end
        eta_OPIT   = eta_OPIT + eta_OPIT_i/N_run;
        rho_OPIT   = rho_OPIT + rho_OPIT_i/N_run;
        
        %% Classical Power-based Subspace Tracking Algorithms
        %% OPAST
        fprintf(' + OPAST:  \n');
        eta_OPAST_i = zeros(1,T);
        rho_OPAST_i = zeros(1,T);
        for ii = 1 : N_run
            t1 = tic;
            [~,eta_ii,rho_ii]  = OPAST_Tracking(X_stream,beta,U_stream);
            t1_end = toc(t1);
            eta_OPAST_i   = eta_OPAST_i + eta_ii;
            rho_OPAST_i   = rho_OPAST_i + rho_ii;
        end
        eta_OPAST  = eta_OPAST  + eta_OPAST_i/N_run;
        rho_OPAST  = rho_OPAST  + rho_OPAST_i/N_run;
        
        
        %% FAPI
        fprintf(' + FAPI:  \n');
        eta_FAPI_i = zeros(1,T);
        rho_FAPI_i = zeros(1,T);
        for ii = 1 : N_run
            t1 = tic;
            [~,eta_ii,rho_ii]  = FAPI_Tracking(X_stream,beta,U_stream);
            t1_end = toc(t1);
            eta_FAPI_i   = eta_FAPI_i + eta_ii;
            rho_FAPI_i   = rho_FAPI_i + rho_ii;
        end
        eta_FAPI  = eta_FAPI  + eta_FAPI_i/N_run;
        rho_FAPI  = rho_FAPI  + rho_FAPI_i/N_run;
        
        %% YAST
        fprintf(' + YAST:  \n');
        eta_GYAST_i = zeros(1,T);
        rho_GYAST_i = zeros(1,T);
        for ii = 1 : N_run
            t1 = tic;
            [~,eta_ii,rho_ii]  = GYAST_Tracking(X_stream,beta,U_stream);
            t1_end = toc(t1);
            eta_GYAST_i   = eta_GYAST_i + eta_ii;
            rho_GYAST_i   = rho_GYAST_i + rho_ii;
        end
        eta_GYAST  = eta_GYAST  + eta_GYAST_i/N_run;
        rho_GYAST  = rho_GYAST  + rho_GYAST_i/N_run;
        
        
        %% LORAF
        fprintf(' + LORAF:  \n');
        eta_LORAF_i = zeros(1,T);
        rho_LORAF_i = zeros(1,T);
        for ii = 1 : N_run
            t1 = tic;
            [~,eta_ii,rho_ii]  = LORAF_Tracking(X_stream,beta,U_stream);
            t1_end = toc(t1);
            eta_LORAF_i   = eta_LORAF_i + eta_ii;
            rho_LORAF_i   = rho_LORAF_i + rho_ii;
        end
        eta_LORAF  = eta_LORAF  + eta_LORAF_i/N_run;
        rho_LORAF  = rho_LORAF  + rho_LORAF_i/N_run;
                
        %% Sparse Subspace Tracking Algorithms
        %% L1-PAST
        fprintf(' + L1-PAST:  \n');
        eta_L1_i = zeros(1,T);
        rho_L1_i = zeros(1,T);
        for ii = 1 : N_run
            OPTS_l1_PAST.lambda = beta;
            t1 = tic;
            [~,eta_ii,rho_ii] = l1_PAST_ST(X_stream,OPTS_l1_PAST,U_stream);
            t1_end = toc(t1);
            eta_L1_i = eta_L1_i + eta_ii;
            rho_L1_i = rho_L1_i + rho_ii;      
        end
        eta_L1  = eta_L1  + eta_L1_i/N_run;
        rho_L1  = rho_L1  + rho_L1_i/N_run;
        
        %% SSFAPI 
        fprintf(' + SS-FAPI:  \n');
        eta_SS_i = zeros(1,T);
        rho_SS_i = zeros(1,T);
        OPTS_SS_FAPI.lambda = beta;
        for ii = 1 : N_run
            t1 = tic;
            [~,eta_ii,rho_ii]   = SS_FAPI_ST(X_stream,OPTS_SS_FAPI,U_stream);
            t1_end = toc(t1);
            eta_SS_i  = eta_SS_i  + eta_ii;
            rho_SS_i  = rho_SS_i  + rho_ii;
        end
      	eta_SS  = eta_SS  + eta_SS_i/N_run;
        rho_SS  = rho_SS  + rho_SS_i/N_run;
          
        %% Oja 
        fprintf(' + Oja: \n');
        eta_Oja_i = zeros(1,T); 
        rho_Oja_i = zeros(1,T); 
        for ii = 1 : N_run
            OPTS_OJA = [];
            t1 = tic;
            [~,eta_ii,rho_ii]   = Oja(X_stream,U_stream);
            t1_end   = toc(t1);
            eta_Oja_i  = eta_Oja_i  + eta_ii;
            rho_Oja_i  = rho_Oja_i  + rho_ii;
        end
        eta_Oja  = eta_Oja  + eta_Oja_i/N_run;
        rho_Oja  = rho_Oja  + rho_Oja_i/N_run;
        
        
        %% SSPCA  
        fprintf(' + SSPCA: \n');
        eta_SSPCA_i = zeros(1,T); 
        rho_SSPCA_i = zeros(1,T);
        OPTS_SSPCA = [];
        OPTS_SSPCA.sparsity = sparsity;
        for ii = 1 : N_run
            t1 = tic;
            [~,eta_i,rho_i]   = SSPCA(X_stream,OPTS_SSPCA,U_stream);
            t1_end   = toc(t1);
            eta_SSPCA_i  = eta_SSPCA_i  + eta_i;
            rho_SSPCA_i  = rho_SSPCA_i  + rho_i;
        end
        eta_SSPCA  = eta_SSPCA  + eta_SSPCA_i/N_run;
        rho_SSPCA  = rho_SSPCA  + rho_SSPCA_i/N_run;
        
        
    end
    
    
    eta0(snr_i,:) = eta_OPIT/N_exp;
    rho0(snr_i,:) = rho_OPIT/N_exp;
    eta1(snr_i,:) = eta_L1/N_exp;
    rho1(snr_i,:) = rho_L1/N_exp;
    eta2(snr_i,:) = eta_Oja/N_exp;
    rho2(snr_i,:) = rho_Oja/N_exp;
    eta3(snr_i,:) = eta_SS/N_exp;
    rho3(snr_i,:) = rho_SS/N_exp;
    eta4(snr_i,:) = eta_SSPCA/N_exp;
    rho4(snr_i,:) = rho_SSPCA/N_exp;
   
    eta5(snr_i,:) = eta_OPAST/N_exp;
    rho5(snr_i,:) = rho_OPAST/N_exp;
    eta6(snr_i,:) = eta_FAPI/N_exp;
    rho6(snr_i,:) = rho_FAPI/N_exp;
    eta7(snr_i,:) = eta_GYAST/N_exp;
    rho7(snr_i,:) = rho_GYAST/N_exp;
    eta8(snr_i,:) = eta_LORAF/N_exp;
    rho8(snr_i,:) = rho_LORAF/N_exp;
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
k = 5;
TS = round(T/5);
% round 1
d2 = semilogy(1:k:T,rho1(1,1:k:end),...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);
d21 = plot(1:TS:T,rho1(1,1:TS:end),...
 'marker','o','markersize',makerSize,...
   'linestyle','none','color',blue_o,'LineWidth',LineWidth);
d22 = semilogy(1:1,rho1(1,1:1),...
    'marker','o','markersize',makerSize,...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);

d3 = semilogy(1:k:T,rho2(1,1:k:end),...
    'linestyle','-','color',gree_n,'LineWidth',LineWidth);
d31 = plot(1:TS:T,rho2(1,1:TS:end),...
 'marker','d','markersize',makerSize,...
   'linestyle','none','color',gree_n,'LineWidth',LineWidth);
d32 = semilogy(1:1,rho2(1,1:1),...
    'marker','d','markersize',makerSize,...
    'linestyle','-','color',gree_n,'LineWidth',LineWidth);

d4 = semilogy(1:k:T,rho3(1,1:k:end),...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);
d41 = plot(1:TS:T,rho3(1,1:TS:end),...
 'marker','p','markersize',makerSize,...
   'linestyle','none','color',gree_o,'LineWidth',LineWidth);
d42 = semilogy(1:1,rho3(1,1:1),...
    'marker','p','markersize',makerSize,...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);

d44 = semilogy(1:k:T,rho4(1,1:k:end),...
    'linestyle','-','color',viol_n,'LineWidth',LineWidth);
d441 = plot(1:TS:T,rho4(1,1:TS:end),...
 'marker','^','markersize',makerSize,...
   'linestyle','none','color',viol_n,'LineWidth',LineWidth);
d442 = semilogy(1:1,rho4(1,1:1),...
    'marker','^','markersize',makerSize,...
    'linestyle','-','color',viol_n,'LineWidth',LineWidth);

d1 = semilogy(1:k:T,rho0(1,1:k:end),...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);
d11 = plot(1:TS:T,rho0(1,1:TS:end),...
 'marker','s','markersize',makerSize,...
   'linestyle','none','color',red_o,'LineWidth',LineWidth);
d12 = semilogy(1:1,rho0(1,1:1),...
    'marker','s','markersize',makerSize,...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);




d5 = semilogy(1:k:T,rho5(1,1:k:end),...
    'linestyle','-','color','k','LineWidth',LineWidth);
d51 = plot(1:TS:T,rho5(1,1:TS:end),...
 'marker','+    ','markersize',makerSize,...
   'linestyle','none','color','k','LineWidth',LineWidth);
d52 = semilogy(1:1,rho5(1,1:1),...
    'marker','+','markersize',makerSize,...
    'linestyle','-','color','k','LineWidth',LineWidth);

d6 = semilogy(1:k:T,rho6(1,1:k:end),...
    'linestyle','--','color',lbrow_n,'LineWidth',LineWidth);
d61 = plot(1:TS:T,rho6(1,1:TS:end),...
 'marker','>','markersize',makerSize,...
   'linestyle','none','color',lbrow_n,'LineWidth',LineWidth);
d62 = semilogy(1:1,rho6(1,1:1),...
    'marker','>','markersize',makerSize,...
    'linestyle','--','color',lbrow_n,'LineWidth',LineWidth);



d7 = semilogy(1:k:T,rho7(1,1:k:end),...
    'linestyle','-','color',yell_n,'LineWidth',LineWidth);
d71 = plot(1:TS:T,rho7(1,1:TS:end),...
 'marker','<','markersize',makerSize,...
   'linestyle','none','color',yell_n,'LineWidth',LineWidth);
d72 = semilogy(1:1,rho7(1,1:1),...
    'marker','<','markersize',makerSize,...
    'linestyle','-','color',yell_n,'LineWidth',LineWidth);

d8 = semilogy(1:k:T,rho8(1,1:k:end),...
    'linestyle','-','color',lblu_n,'LineWidth',LineWidth);
d81 = plot(1:TS:T,rho8(1,1:TS:end),...
 'marker','x','markersize',makerSize,...
   'linestyle','none','color',lblu_n,'LineWidth',LineWidth);
d82 = semilogy(1:1,rho8(1,1:1),...
    'marker','x','markersize',makerSize,...
    'linestyle','-','color',lblu_n,'LineWidth',LineWidth);




lgd = legend([d52 d22 d82 d42 d62 d72 d32  d442 d12],...
      '\texttt{OPAST}','\texttt{L1-PAST}','\texttt{LORAF}',...
      '\texttt{FAPI}','\texttt{SSFAPI}','\texttt{GYAST}','\texttt{Oja}','\texttt{SSPCA}','\texttt{OPIT}');
set(lgd,'Interpreter','latex',...
     'FontSize',22,'NumColumns',3);
    % 'FontSize',22,'NumColumns',3,'EdgeColor', 'none');
xlabel('Time Index - $t$','interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel('$\sin \theta(\mathbf{A}_{t},\mathbf{U}_{t}) $','interpreter','latex','FontSize',13,'FontName','Times New Roman');

h=gca;
set(gca, 'YScale', 'log')
set(h,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(h,'Xtick',0:round(T/5):T,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');
set(h,'FontSize', 30);
axis([0 T 0.09*SNR_L 10]);
grid on;
box on;
set(fig, 'units', 'inches', 'position', [0.5 0.5 8 7]);


