
function [data, shock] = main_load_data_rev(keep_c, initial_tariff_index, factor_index)
%% Solve System of Linear equations with Allen Elasticities
% Code by Sihwan Yang (PhD student in UCLA Economics), 02/08/2022

%% Replication Code Description
% **1. main_load_data_rev: First part of main code. Calculate inputs from data
% 2. main_dlogW_rev: Second part of main code. Solve system of linear equations to derive change in welfare
% 3. AES_func: Calculate Allen Elasticity of Substitution
% 4. Nested_CES_linear_final_rev: Solve system of linear equations by inverting the system
% 5. Nested_CES_linear_result_final: Calculate derivatives that are used to iterate or derive welfare change

%% Input
% C : The number of countries (In WIOD 2008, C = 41)
% N : The number of sectors per country (In WIOD 2008, N = 30)
% F : The number of factors per country (In WIOD 2008, F = 4)
% CN : The number of total sectors in the economy (CN = C*N)
% CF : The number of total factors in the economy (CF = C*F)

% Data
% Omega_tilde : Input-Output Matrix -  Matrix of size (CN,CN)
% beta : Household expenditure share - Matrix of size (CN,C)
% alpha : Value-added share - Vector of size (CN,1)
% alpha_VA : Factor share - Matrix of size (CN,F)
% GNE_weights: GNE share - Vector of size (C,1)
% trade_elast: Elasticity of subsititution across countries engaging in trade

% Initial tariff (optional)
% Tariffs_cons_matrix_new : Tariff matrix when household (column) buys goods (row) - Matrix of size (CN,C)
% Tariff_matrix_new : Tariff matrix when a sector (row) buys goods (column) - Matrix of size (CN,CN)
% init_t : Initial tariff - Matrix of size (C+CN,CN+CF)

%% Output
% data : Data struct
% shock : Shock struct

%% 1. Load Input-Output Data and choose initial tariff option
% Input Mat files (wiod2008data_without_tariff.mat, wiod2008data_initial_tariff.mat)
% come from WIOD of year 2008 with some data cleaning and reordering

% wiod2008data_without_tariff.mat -- Data without initial tariff
% (i,j) element of Omega: expenditure share of sector i on sector j
% (i,c) element of beta: expenditure share of household c on sector i
% (i,f) element of alpha_VA: expenditure share of sector i on factor f ouf of factor usage
% (i,1) element of alpha: value-added share of sector i
% (c,1) element of GNE_weights: GNE share of country c in world GNE
% Check: Each row of Omega (=Omega_tilde) should sum up to 1

% wiod2008data_initial_tariff.mat -- Data with initial tariff
% Under initial tariff, one should calculate Omega_tilde matrix first by
% multiplying expenditure share (IO_new) with initial tariff (Tariff_matrix_new) and normalizing
% Similarly, alpha considers initial tariff information for after tariff expenditure
% and beta, GNE_weights consider it for after tariff consumption
% Element-wise description of beta, alpha_VA, alpha, GNE_weights is same as without initial tariff
% Check: Each row of Omega_tilde should sum up to 1

% initial_tariff_index = 1; % without initial tariff -- 1, with tariff -- 2

if initial_tariff_index == 1
    [Omega, beta, alpha_VA, alpha, trade_elast, GDP_weights] = IO_reorder(keep_c);
    
    Omega_tilde = Omega;
    GNE_weights = GDP_weights;
elseif initial_tariff_index == 2

    [Omega, beta, alpha_VA, alpha, trade_elast, GDP_weights, Tariff_matrix_new, Tariffs_cons_matrix_new] = IO_reorder_init_tariff(keep_c);

    Omega_tilde = Omega;
    GNE_weights = GDP_weights;
end

alpha_VA_temp = alpha_VA;
alpha_VA(:,3) = alpha_VA_temp(:,4);
alpha_VA(:,4) = alpha_VA_temp(:,3);

C = length(keep_c)+1; % countries in keep_c + ROW
N = 30; % sectors in WIOD 2008

if factor_index == 1
    F = 4; % 1st factor: capital, 2nd factor: low-skilled, 3rd factor: medium-skilled, 4th factor: high-skilled
elseif factor_index == 2
    F = N; % Country-sector-specific labor
end

CN = C*N;
CF = C*F;

% (i,c) element of Tariffs_cons_matrix_new: tariff rate of household c (destination) on sector i (origin)
% (i,j) element of Tariff_matrix_new: tariff rate of sector i (destination) on sector j (origin)
% (i,j) element of init_t: initial tariff when household/sector i buys from sector/factor j

if initial_tariff_index == 1
    init_t = [ones(C+CN,CN), ones(C+CN,CF)];
elseif initial_tariff_index == 2
    init_t = [[Tariffs_cons_matrix_new' ; Tariff_matrix_new], ones(C+CN,CF)];
end

Omega_tilde = Omega_tilde(1:CN,1:CN);
Omega_tilde(isnan(Omega_tilde)) = 0;
Omega_tilde(isinf(Omega_tilde)) = 1;
beta = beta(1:CN,1:C);
beta(isnan(beta)) = 0;
beta(isinf(beta)) = 1;
alpha = alpha(1:CN,1);
alpha(isnan(alpha)) = 0;
alpha(isinf(alpha)) = 1;
alpha_VA = alpha_VA(1:CN,:);
alpha_VA(isnan(alpha_VA)) = 0;
alpha_VA(isinf(alpha_VA)) = 1;
GNE_weights = GNE_weights(1:C,1);

%% 2. Transform Raw Data into Input share
% Below variables in paranthesis follow notation in Appendix C
% beta_s (b_ic) : How much HH from country c consumes sector i good
% beta_disagg (delta_im^0c) : How much HH from country c consumes sector i from country m
% alpha (alpha_ic) : How much sector i from country c uses value-added
% alpha_VA (alpha_f^ic) : How much sector i from country c uses factor f
% Omega_s (omega_j^ic) : How much sector i from country c uses sector j
% Omega_disagg (delta_jm^ic) : How much sector i from country c uses sector j from country m

Omega_s = zeros(CN,N);
beta_s = zeros(C,N);

for i=1:C
   beta_s = beta_s + beta((i-1)*N+1:i*N,:)';
   Omega_s = Omega_s + Omega_tilde(:,(i-1)*N+1:i*N);
end

Omega_disagg = cell(CN,N);
beta_disagg = cell(C,N);

beta_disagg_temp = beta'./repmat(beta_s, 1, C);
Omega_disagg_temp = Omega_tilde./repmat(Omega_s, 1, C);

for i=1:C
    for j=1:N
        for k=1:C
            beta_disagg{i,j}(1,k) = beta_disagg_temp(i,(k-1)*N+j);
        end
    end
end

for i=1:CN
    for j=1:N
        for k=1:C
           Omega_disagg{i,j}(1,k) = Omega_disagg_temp(i,(k-1)*N+j);
        end
    end
end

% Omega_total_C (Omega_jm^0c) : Hom much c buys from jm - Matrix of size (C,CN)
% Omega_total_N (Omega_jm^ic) : How much ic buys from jm (j can be good/factor) - Matrix of size (CN,CN+CF)

Omega_total_C = zeros(C,CN); % same as beta'
Omega_total_N = zeros(CN,CN+CF);

for c=1:C
    for i=1:CN
        ic = floor((i-1)/N)+1;
        in = mod(i,N);
        if mod(i,N)==0
            in = N;
        end
        Omega_total_C(c,i) = beta_s(c,in)*beta_disagg{c,in}(1,ic);
    end
end

for i=1:CN
    ic = floor((i-1)/N)+1;
    
    for j=1:CN
        jc = floor((j-1)/N)+1;
        jn = mod(j,N);
        if mod(j,N)==0
            jn = N;
        end
        Omega_total_N(i,j) = (1-alpha(i,1))*Omega_s(i,jn)*Omega_disagg{i,jn}(1,jc);
    end
    if factor_index == 1
        Omega_total_N(i,CN+(ic-1)*F+1:CN+ic*F) = alpha(i,1)*alpha_VA(i,:);
    elseif factor_index == 2
        Omega_total_N(i,CN+i) = alpha(i,1);
    end
end

%% 3. Define Standard Form of Omega, lambda, and chi
L = C+CN+CF;
% Omega_total_tilde: Standard Form of Cost-based IO matrix - Matrix of size (L,L)
% (c,i) element of First C row of Omega_total_tilde: How much household c buys from sector i
% (i,j) element of C+1:C+CN row of Omega_total_tilde: How much sector i buys from sector/factor j
% Since factors are endowed, last CF rows of Omega_total_tilde are zeros
% Check: Each row of Omega_total_tilde should sum up to 1

Omega_total_tilde = [ zeros(C,C), Omega_total_C, zeros(C,CF) ;
                zeros(CN,C), Omega_total_N ; 
                zeros(CF,L)];

% Omega_total: Standard Form of Revenue-based IO matrix - Matrix of size (L,L)
% Conceptually, Omega_total = Omega_total_tilde/initial_tariff;

if initial_tariff_index == 1
    Omega_total = Omega_total_tilde; 
elseif initial_tariff_index == 2
    Omega_total = Omega_total_tilde;
    Omega_total(1:C,C+1:C+CN) = Omega_total_tilde(1:C,C+1:C+CN)./Tariffs_cons_matrix_new';
    Omega_total(C+1:C+CN,C+1:C+CN) = Omega_total_tilde(C+1:C+CN,C+1:C+CN)./Tariff_matrix_new;
end

% Psi_total : Standard form of Leontief inverse - Matrix of size )L,L)
Psi_total_tilde = linsolve(eye(L)-Omega_total_tilde,eye(L));
Psi_total = linsolve(eye(L)-Omega_total,eye(L));

% chi_std : Standard form of Household income - Vector of size (L,1)
% lambda_std : Standard form of Sales share/Factor income share - Vector of size (L,1)
chi_std = [GNE_weights ; zeros(CN+CF,1)];
lambda_std = (chi_std'*Psi_total)';

lambda_CN = lambda_std(1:C+CN,1);
lambda_F = lambda_std(C+CN+1:end,1);

%% 4. Store data

data.C = C;
data.N = N;
data.F = F;
data.CN = CN;
data.CF= CF;
data.L = L;

data.trade_elast = trade_elast;
%data.trade_elast = 5*ones(N,1);
data.alpha = alpha;
data.beta_s = beta_s;
data.Omega_s = Omega_s;
data.beta_disagg = beta_disagg;
data.Omega_disagg = Omega_disagg;
data.Omega_total_C = Omega_total_C;
data.Omega_total_N = Omega_total_N;

data.Omega_total_tilde = Omega_total_tilde;
data.Omega_total = Omega_total;
data.Psi_total_tilde = Psi_total_tilde;
data.Psi_total = Psi_total;
data.lambda_std = lambda_std;
data.lambda_CN = lambda_CN;
data.lambda_F = lambda_F;
data.chi_std = chi_std;

shock.init_t = init_t;
end