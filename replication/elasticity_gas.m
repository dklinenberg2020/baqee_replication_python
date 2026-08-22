clear all; clc; close all;

% Set latex as default interpreter
set(groot, 'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

pE = 1; pX = 1;
%alpha = 0.025;
%alpha = 0.04; %total GNE share of gas+oil+coal
%alpha = alpha*0.34;%gas = 905/(1077 + 905 + 606) = 34% of gas+oil+coal
alpha = 0.01; 

E_drop = 0.3;
%E_drop = 0.15;
%E_drop = 0.5;
E_after_drop = 1 - E_drop;

%for Leontief case need E/alpha = X/(1-alpha). Want E=1 to be the starting
%point. So let's set the endowment of the other factor to X = (1-alpha)/alpha.

X = (1-alpha)/alpha;

Emin = E_after_drop-0.05; Emax = 1.05; N=501;
E = linspace(Emin,Emax,N)';

N_sigma = 3;
sigmas = [0.1,0.2,0.999];
Y = zeros(N,N_sigma);
MPE = zeros(N,N_sigma);
MPX = zeros(N,N_sigma);
P = zeros(N,N_sigma);
exp_share = zeros(N,N_sigma);

MPE_CD = zeros(N,1);
MPX_CD = zeros(N,1);
MPE_Leontief = zeros(N,1);
MPX_Leontief = zeros(N,1);
exp_share_Leontief = zeros(N,1);

for i=1:N_sigma
    sigma = sigmas(i);
    Y(:,i) = (alpha^(1/sigma)*E.^((sigma-1)/sigma)  + (1-alpha)^(1/sigma)*X^((sigma-1)/sigma)).^(sigma/(sigma-1));
    Ybench(i) = (alpha^(1/sigma)  + (1-alpha)^(1/sigma)*X^((sigma-1)/sigma)).^(sigma/(sigma-1));
    opt_exp_share(i) = alpha*pE^(1-sigma)/(alpha*pE^(1-sigma) + (1-alpha)*pX^(1-sigma));
    MPE(:,i) = (alpha^(1/sigma)*E.^((sigma-1)/sigma)  + (1-alpha)^(1/sigma)*X^((sigma-1)/sigma)).^(sigma/(sigma-1)-1)*alpha^(1/sigma).*E.^(-1/sigma);
    MPX(:,i) = (alpha^(1/sigma)*E.^((sigma-1)/sigma)  + (1-alpha)^(1/sigma)*X^((sigma-1)/sigma)).^(sigma/(sigma-1)-1)*(1-alpha)^(1/sigma).*X.^(-1/sigma);
    P(:,i) = (alpha*MPE(:,i).^(1-sigma) + (1-alpha)*pX^(1-sigma)).^(1/(1-sigma));
    exp_share(:,i) = (MPE(:,i).*E)./(P(:,i).*Y(:,i));
end

Y_CD = E.^alpha.*X.^(1-alpha);
Ybench_CD = 1.^alpha.*X.^(1-alpha);
MPE_CD = alpha.*E.^(alpha-1).*X.^(1-alpha);
MPX_CD = (1-alpha).*E.^(alpha).*X.^(-alpha);
exp_share_CD = alpha*ones(N,1);

Y_Leontief = min(E/alpha,X/(1-alpha));
Ybench_Leontief = min(1/alpha,X/(1-alpha));
MPE_Leontief(E<1) = 1/alpha;
MPE_Leontief(E==1) = 1;
MPE_Leontief(E>1) = 0;
MPX_Leontief(E<1) = 0;
MPX_Leontief(E==1) = 1;
MPX_Leontief(E>1) = 1/(1-alpha);
exp_share_Leontief(E<1) = 1;
exp_share_Leontief(E==1) = alpha;
exp_share_Leontief(E>1) = 0;



figure(1)
plot(E,Y_Leontief./Ybench_Leontief,'--','LineWidth',2)
hold on
plot(E,Y(:,1)./Ybench(1),'-.','LineWidth',2)
hold on
plot(E,Y(:,2)./Ybench(2),'-.','LineWidth',2)
hold on
% plot(E,Y(:,3)./Ybench(3),'LineWidth',2)
% hold on
plot(E,Y_CD./Ybench_CD,'LineWidth',2)
line([E_after_drop E_after_drop],[Emin 1.02],'Color','k','LineWidth',1.5)
legend('$\sigma=0$ (Leontief)','$\sigma=0.1$','$\sigma=0.16$','$\sigma=1$ (Cobb-Douglas)','Fontsize',16,'Location','SouthEast')
grid
xlim([Emin Emax])
ylim([Emin 1.02])
set(gca,'FontSize',16)
xlabel('Energy, $E$','Fontsize',20)
ylabel('Production, $Y$','Fontsize',20)
print -dpdf elasticity_fig_gas.pdf
print -depsc elasticity_fig_gas.eps

[value index_bench] = min(abs(E-1));

figure(2)
plot(E,MPE_Leontief./MPE_Leontief(index_bench),'--','LineWidth',2)
hold on
plot(E,MPE(:,1)./MPE(index_bench,1),'-.','LineWidth',2)
hold on
plot(E,MPE(:,2)./MPE(index_bench,1),'-.','LineWidth',2)
% hold on
% plot(E,MPE(:,3)./MPE(index_bench,1),'LineWidth',2)
hold on
plot(E,MPE_CD./MPE_CD(index_bench),'LineWidth',2)
line([E_after_drop E_after_drop],[0.5 10],'Color','k','LineWidth',1.5)
legend('$\sigma=0$ (Leontief)','$\sigma=0.1$','$\sigma=0.16$','$\sigma=1$ (Cobb-Douglas)','Fontsize',16,'Location','NorthEast')
grid
xlim([Emin Emax])
ylim([0.5 10])
set(gca,'FontSize',16)
xlabel('Energy, $E$','Fontsize',20)
ylabel('Price $p_E$ ($=MPE$) rel. to baseline','Fontsize',20)
print -dpdf MPE_fig_gas.pdf
print -depsc MPE_fig_gas.eps

figure(3)
plot(E,MPX_Leontief./MPX_Leontief(index_bench),'--','LineWidth',2)
hold on
plot(E,MPX(:,1)./MPX(index_bench,1),'-.','LineWidth',2)
hold on
plot(E,MPX(:,2)./MPX(index_bench,1),'-.','LineWidth',2)
% hold on
% plot(E,MPE(:,3)./MPE(index_bench,1),'LineWidth',2)
hold on
plot(E,MPX_CD./MPX_CD(index_bench),'LineWidth',2)
line([E_after_drop E_after_drop],[0 1.05],'Color','k','LineWidth',1.5)
legend('$\sigma=0$ (Leontief)','$\sigma=0.1$','$\sigma=0.16$','$\sigma=1$ (Cobb-Douglas)','Fontsize',16,'Location','SouthEast')
grid
xlim([Emin Emax])
ylim([0 1.05])
set(gca,'FontSize',16)
xlabel('Energy, $E$','Fontsize',20)
ylabel('Price $p_X$ ($=MPX$) rel. to baseline','Fontsize',20)
print -dpdf MPX_fig_gas.pdf
print -depsc MPX_fig_gas.eps


figure(4)
plot(E,exp_share_Leontief,'--','LineWidth',2)
hold on
plot(E,exp_share(:,1),'-.','LineWidth',2)
hold on
plot(E,exp_share(:,2),'-.','LineWidth',2)
hold on
% plot(E,exp_share(:,3),'LineWidth',2)
% hold on
plot(E,exp_share_CD,'LineWidth',2)
line([E_after_drop E_after_drop],[0 1],'Color','k','LineWidth',1.5)
legend('$\sigma=0$ (Leontief)','$\sigma=0.1$','$\sigma=0.16$','$\sigma=1$ (Cobb-Douglas)','Fontsize',16,'Location','NorthEast')
grid
xlim([Emin Emax])
ylim([0 1])
set(gca,'FontSize',16)
xlabel('Energy, $E$','Fontsize',20)
ylabel('Expenditure share on energy $E$','Fontsize',20)
print -dpdf exp_share_fig_gas.pdf
print -depsc exp_share_fig_gas.eps


[value index] = min(abs(E-E_after_drop));

disp('Elasticity = 0 (Leontief)')
disp('Output loss, pE^new/pE^old, pX^new/pX^old, new energy share')
disp([Y_Leontief(index)./Ybench_Leontief-1,Inf,0,1])

disp('Elasticity = 0.1')
disp('Output loss, pE^new/pE^old, pX^new/pX^old, new energy share')
disp([Y(index,1)./Ybench(1)-1,MPE(index,1)./MPE(index_bench,1),MPX(index,1)./MPX(index_bench,1),exp_share(index,1)])

disp('Elasticity = 0.2')
disp('Output loss, pE^new/pE^old, pX^new/pX^old, new energy share')
disp([Y(index,2)./Ybench(2)-1,MPE(index,2)./MPE(index_bench,2),MPX(index,2)./MPX(index_bench,2),exp_share(index,2)])

disp('Elasticity = 1 (Cobb-Douglas)')
disp('Output loss, pE^new/pE^old, pX^new/pX^old, new energy share')
disp([Y_CD(index)./Ybench_CD-1,MPE_CD(index)./MPE_CD(index_bench),MPX_CD(index)./MPX_CD(index_bench),alpha])
