
function [dlambda_result, dlambda_F_star_result, dchi_std, dlogP_Vec, dlogP_Mat, dOmega_total, dPsi_total, dOmega_total_tilde, dPsi_total_tilde] = Nested_CES_linear_result_final(dlambda_F_all, data, shock)
%% Implement Linearization model of Nested CES
% Code by Sihwan Yang (PhD student in UCLA Economics), 10/21/2020

%% Replication Code Description
% 1. main_load_data: First part of main code. Calculate inputs from data
% 2. main_dlogW: Second part of main code. Solve system of linear equations to derive change in welfare
% 3. AES_func: Calculate Allen Elasticity of Substitution
% 4. Nested_CES_linear_final: Solve system of linear equations by inverting the system
% **5. Nested_CES_linear_result_final: Calculate derivatives that are used to iterate or derive welfare change

%% Input

C = data.C;
N = data.N;
F = data.F;
CN = data.CN;
CF= data.CF;
L = data.L;

Omega_total = data.Omega_total;
Psi_total = data.Psi_total;
Omega_total_tilde = data.Omega_total_tilde;
Psi_total_tilde = data.Psi_total_tilde;
lambda_CN = data.lambda_CN;
lambda_F = data.lambda_F;
chi_std = data.chi_std;

Phi_F = data.Phi_F;
Phi_T = data.Phi_T;

AES_C_Mat = data.AES_C_Mat;
AES_N_Mat = data.AES_N_Mat;
AES_F_Mat = data.AES_F_Mat;

dlogt = shock.dlogt;
dlogtau = shock.dlogtau;
init_t = shock.init_t;
dX = shock.dX;

dlambda_F = dlambda_F_all(1:CF,1); 
dlambda_F_star = dlambda_F_all(CF+1:end,1);

%% dPsi_total = f(dlambda_F)

dlogP_F = dlambda_F./lambda_F;

% Equation (4): dlogP_F -> dlogP_CN
dlogP_CN = Psi_total_tilde(1:C+CN,1:C+CN)*dX + Psi_total_tilde(1:C+CN,C+CN+1:end)*dlogP_F;
dlogP_Vec = [dlogP_CN ; dlogP_F]; % For HH and Producer: Marginal cost, For factor: Factor price

dlogP_Mat = zeros(L,L); % bilateral price - Matrix of size (C+CN+CF,C+CN+CF)
dlogP_Mat(1:C+CN,C+1:C+CN) = dlogt(:,1:CN) + dlogtau(:,1:CN) + repmat(dlogP_CN(C+1:end,1)',C+CN,1); 
for i=1:CN
    ic = floor((i-1)/N)+1;
    dlogP_Mat(C+i,C+CN+(ic-1)*F+1:C+CN+ic*F) = dlogt(C+i,CN+(ic-1)*F+1:CN+ic*F) + dlogtau(C+i,CN+(ic-1)*F+1:CN+ic*F) + dlogP_F((ic-1)*F+1:ic*F,1)';
end

% Equation (3): dlogP -> dOmega_tilde
dOmega_CN_tilde = zeros(C+CN,CN+CF);
for c=1:C
    temp = sum(repmat(Omega_total_tilde(c,C+1:C+CN),CN,1).*(AES_C_Mat(:,:,c)-1).*repmat(dlogP_Mat(c,C+1:C+CN),CN,1),2);
    dOmega_CN_tilde(c,1:CN) = Omega_total_tilde(c,C+1:C+CN).*dlogP_Mat(c,C+1:C+CN) + Omega_total_tilde(c,C+1:C+CN).*temp';
end

for kc=1:CN
    temp1 = sum(repmat(Omega_total_tilde(C+kc,C+1:end),CN,1).*(AES_N_Mat(:,:,kc)-1).*repmat(dlogP_Mat(C+kc,C+1:end),CN,1),2);
    temp2 = sum(repmat(Omega_total_tilde(C+kc,C+1:end),CF,1).*(AES_F_Mat(:,:,kc)-1).*repmat(dlogP_Mat(C+kc,C+1:end),CF,1),2);
    temp = [temp1 ; temp2];
    dOmega_CN_tilde(C+kc,:) = Omega_total_tilde(C+kc,C+1:end).*dlogP_Mat(C+kc,C+1:end) + Omega_total_tilde(C+kc,C+1:end).*temp';
end
dOmega_total_tilde = [zeros(C+CN,C), dOmega_CN_tilde ; 
                zeros(CF,L)];

% Equation (2): dOmega_tilde -> dOmega
dOmega_CN = Omega_total(1:C+CN,C+1:end)./Omega_total_tilde(1:C+CN,C+1:end).*dOmega_CN_tilde - Omega_total(1:C+CN,C+1:end).*dlogt;
dOmega_CN(isnan(dOmega_CN)) = 0;

% Equation (1): dOmega -> dPsi
dOmega_total = [zeros(C+CN,C), dOmega_CN ; 
                zeros(CF,L)];
dPsi_total = Psi_total*dOmega_total*Psi_total;
dPsi_total_tilde = Psi_total_tilde*dOmega_total_tilde*Psi_total_tilde;

%% dchi_std = g(dlambda_F, dlambda_F_star)
% Equation (5): dlambda_F, dlambda_F_star -> dchi
dt = init_t.*dlogt; % Matrix of Size (C+CN,CN)

temp1 = (lambda_CN'*squeeze(sum(Phi_T.*repmat(Omega_total(1:C+CN,C+1:end),[1 1 C]).*repmat(dt,[1 1 C]),2)))';
temp2 = (lambda_CN'*squeeze(sum(Phi_T.*repmat(dOmega_total(1:C+CN,C+1:end),[1 1 C]).*repmat(init_t-1,[1 1 C]),2)))';
dchi = Phi_F*dlambda_F + temp1 + temp2 + dlambda_F_star;

dchi_std = [dchi ; zeros(CN+CF,1)];

%% Output

dlambda_result = (dchi_std'*Psi_total + chi_std'*dPsi_total)';
dlambda_CN_result = dlambda_result(1:C+CN,1);

dlambda_F_star_result = (dlambda_CN_result'*squeeze(sum(Phi_T.*repmat(Omega_total(1:C+CN,C+1:end),[1 1 C]).*repmat(init_t-1,[1 1 C]),2)))';

end