
function [AES_C_Mat, AES_N_Mat, AES_F_Mat] = AES_func(data)
%% Derive Allen Elasticity of Substitution with Nested CES Model
% Code by Sihwan Yang (PhD student in UCLA Economics), 10/21/2020

%% Replication Code Description
% 1. main_load_data: First part of main code. Calculate inputs from data
% 2. main_dlogW: Second part of main code. Solve system of linear equations to derive change in welfare
% **3. AES_func: Calculate Allen Elasticity of Substitution
% 4. Nested_CES_linear_final: Solve system of linear equations by inverting the system
% 5. Nested_CES_linear_result_final: Calculate derivatives that are used to iterate or derive welfare change

%% Input
% C : The number of countries (In WIOD 2008, C = 41)
% N : The number of sectors per country (In WIOD 2008, N = 30)
% F : The number of factors per country (In WIOD 2008, F = 4)
% CN : The number of total sectors in the economy (CN = C*N)
% CF : The number of total factors in the economy (CF = C*F)

% Data
% alpha (alpha_ic) : How much sector i from country c uses value-added
% beta_s (b_ic) : How much HH from country c consumes sector i good
% Omega_s (omega_j^ic) : How much sector i from country c uses sector j
% Omega_total_C (Omega_jm^0c) : Hom much c buys from jm - Matrix of size (C,CN)
% Omega_total_N (Omega_jm^ic) : How much ic buys from jm (j can be good/factor) - Matrix of size (CN,CN+CF)

% Elasticities
% sigma : Elasticity of substitution across consumption goods
% theta : Elasticity of substitution across Intermediate and VA
% gamma : Elasticity of subsititution across factors
% epsilon :  Elasticity of subsititution across intermediate goods
% trade_elast : Elasticity of subsititution across countries engaging in trade

%% Output
% Below variables in paranthesis follow notation in Appendix C
% Indexes are ordered a country-sector (e.g. (c1,s1),(c1,s2),...(c2,s1),(c2,s2),...)
% AES_C_Mat (theta_0c(ic',jm)) - Matrix of size (CN,CN,C)
% (ic',jm,c) element : AES of HH in country c that substitues good i from country c' to good j from country m

% AES_N_Mat (theta_kc(ic',jm)) - Matrix of size (CN,CN+CF,CN)
% (ic',jm,kc) element : AES of producer k in country c that substitutes good i from country c' to good/factor j from country m

% AES_F_Mat (theta_kc(fc,jm)) - Matrix of size (CF,CN+CF,CN)
% (fc,jm,kc) element : AES of producer k in country c that substitues factor f in country c to good/factor j from country m

%% Load Data

C = data.C;
N = data.N;
F = data.F;
CN = data.CN;
CF= data.CF;
L = data.L;

sigma = data.sigma;
theta = data.theta;
gamma = data.gamma;
epsilon = data.epsilon;

trade_elast = data.trade_elast;
alpha = data.alpha;
beta_s = data.beta_s;
Omega_s = data.Omega_s;

Omega_total_C = data.Omega_total_C;
Omega_total_N = data.Omega_total_N;

%% Allen Elasticity of substitution (AES)

AES_C_Mat = sigma*ones(CN,CN,C); % Consumer
AES_N_Mat = theta*ones(CN,CN+CF,CN); % Producer that substitute goods to goods/factors
AES_F_Mat = theta*ones(CF,CN+CF,CN); % Producer that substitute factors to goods/factors

for c=1:C
    for n=1:N
        AES_C_Mat(n:N:end,n:N:end,c) = (trade_elast(n)+1)/beta_s(c,n)+sigma*(1-1/beta_s(c,n)); 
    end
    
    for i=1:CN
        ic = floor((i-1)/N+1);
        in = mod(i,N);
        if mod(i,N)==0
            in = N;
        end

        AES_C_Mat(i,i,c) = -(trade_elast(in)+1)/Omega_total_C(c,i)+(trade_elast(in)+1)/beta_s(c,in)+sigma*(1-1/beta_s(c,in));
        if Omega_total_C(c,i) == 0
            AES_C_Mat(i,i,c) = 1;
        end
    end
    
end

for kc=1:CN

    AES_N_Mat(:,1:CN,kc) = epsilon/(1-alpha(kc,1))+theta*(1-1/(1-alpha(kc,1)));
    
    for n=1:N
        AES_N_Mat(n:N:end,n:N:CN,kc) = (trade_elast(n)+1)/((1-alpha(kc,1))*Omega_s(kc,n))+...
            epsilon*(1/(1-alpha(kc,1))-1/((1-alpha(kc,1))*Omega_s(kc,n)))+theta*(1-1/(1-alpha(kc,1)));
    end
    
    for i=1:CN
        ic = floor((i-1)/N+1);
        in = mod(i,N);
        if mod(i,N)==0
            in = N;
        end

        AES_N_Mat(i,i,kc) = -(trade_elast(in)+1)/Omega_total_N(kc,i)+(trade_elast(in)+1)/((1-alpha(kc,1))*Omega_s(kc,in))+...
            epsilon*(1/(1-alpha(kc,1))-1/((1-alpha(kc,1))*Omega_s(kc,in)))+theta*(1-1/(1-alpha(kc,1)));
        if Omega_total_N(kc,i) == 0
            AES_N_Mat(i,i,kc) = 1;
        end
    end
    AES_N_Mat(:,CN+1:end,kc) = theta;
    
end

for kc=1:CN
    
    kcc = floor((kc-1)/N)+1;
    
    AES_F_Mat(:,1:CN,kc) = theta;
    AES_F_Mat(:,CN+1:end,kc) = 1;
    
    AES_F_Mat((kcc-1)*F+1:kcc*F,CN+(kcc-1)*F+1:CN+kcc*F,kc) = gamma/alpha(kc,1)+theta*(1-1/alpha(kc,1));
    for f=1:F
        AES_F_Mat((kcc-1)*F+f,CN+(kcc-1)*F+f,kc) = -gamma/Omega_total_N(kc,CN+(kcc-1)*F+f)+gamma/alpha(kc,1)+theta*(1-1/alpha(kc,1));
        
        if Omega_total_N(kc,CN+(kcc-1)*F+f)==0
            AES_F_Mat((kcc-1)*F+f,CN+(kcc-1)*F+f,kc) = 1;
        end
    end

end

end
