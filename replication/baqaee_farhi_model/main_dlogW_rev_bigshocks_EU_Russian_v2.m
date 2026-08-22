%% Solve System of Linear equations with Allen Elasticities
% Code by Sihwan Yang (PhD student in UCLA Economics), 03/05/2022
clear
clc
%% Replication Code Description (Revised version)
% 1. main_load_data_rev : First part of main code. Calculate inputs from data
% **2. main_dlogW_rev: Second part of main code. Solve system of linear equations to derive change in welfare
% 3. AES_func: Calculate Allen Elasticity of Substitution
% 4. Nested_CES_linear_final_rev: Solve system of linear equations by inverting the system
% 5. Nested_CES_linear_result_final: Calculate derivatives that are used to iterate or derive welfare change

%% Output
% dlogW : Change in real income of each country for each iteration of discretized shocks
% dlogW_sum : Change in real income of each country from linearized system
% dlogW_world : Change in real income of world
% dlogR : Reallocation term of each country for each iteration of discretized shocks
% dlogR_sum : Reallocation term of each country from linearized system
% dlogY_2nd : Change in world GDP to a 2nd order
tic
%% User Input

% (1) Initial tariff environment : initial_tariff_index
% initial_tariff_index = 1 -- without initial tariff
% initial_tariff_index = 2 -- with initial tariff 
initial_tariff_index = 1;

% (2) Shock type : shock_index, ngrid, intensity
% shock_index = 1 -- universal iceberg cost shock
% shock_index = 2 -- tariff shock
% shock_index = 3 -- EU iceberg cost on Russia
% shock_index = 4 -- EU tariff on Russian energy & transport
shock_index = 3;

% ngrid : number of grids used in for-loop
ngrid = 20; 

% e.g. 10 -- 10% shock
intensity = 150; 

% (3) Choice of # of countries : keep_c
% keep_c = [28 7 41]; % MEX (28), CHN (7), USA (40->41), ROW (35)
% keep_c = [6 7 10 15 16 23 41]; % CAN (6), CHN (7), DEU (10), FRA (15), GBR (16), JPN (23), USA (40->41), ROW (35)
keep_c = (1:41); % All countries
keep_c(35) = []; % Always exclude 35 for ROW

% (4) Choce of factors : factor_index
% factor_index = 1 -- country-specific factors (4 factors per country)
% factor_index = 2 -- country-sector-specific factors (30 factors per country)
factor_index = 2;

% (5) Countries/sectors of interest
EU = [2	3	4	8	9	10	11	12	13	14	15	17	18	21	22	25	26	27	29	30	31	32	33	35	36	37]; % EU : 26 countries
RUS = 34; % Russia : 34th country in WIOD
% DEU = 10; % Germany : 10th country in WIOD
Energy = 7; % Coke, Refined Petroleum and Nuclear Fuel : 7th sector in WIOD
Gas = 15; 
Transport = 20; % Inland transport : 20th sector

% coi = EU; % Country of Interest : EU
coi = EU;
% soi = [Energy, Transport]; % Sector of Interest : Energy, Transport
soi = 1:30; % Sector of Interest : All sectors

% AB_analysis = 0 -- does not conduct analysis, AB_analysis = 1 -- conduct
AB_analysis = 0;

% data, shock: struct from main_load_data.m
[data, shock] = main_load_data_rev(keep_c, initial_tariff_index, factor_index);

C = data.C;
N = data.N;
F = data.F;
CN = data.CN;
CF= data.CF;
L = data.L;

Omega_total = data.Omega_total;
lambda_std = data.lambda_std;
lambda_CN = data.lambda_CN;
chi_std = data.chi_std;
GNE_weights = data.chi_std(1:C,1);

sigma = 0.9; % Elasticity of substitution (1) across Consumption
theta = 0.05; % Elasticity of substitution (2) across Composite Value-added and Intermediates
gamma = 0.5; % Elasticity of substitution (3) across Primary Factors
epsilon = 0.05; % Elasticity of substitution (4) across Intermediate Inputs


data.sigma = sigma;
data.theta = theta;
data.gamma = gamma;
data.epsilon = epsilon;

% Phi_F : Ownership structure of factor - Matrix of size (C,CF)
% Phi_T : Ownership structure of tariff revenue - Matrix of size (C+CN,CN+CF,C)
% (c,f) element of Phi_F: factor income share of factor f owned by household c
% (i,j,c) element of Phi_T: tariff revenue share owned by household c when household/sector i buys from sector/factor j

Phi_F = zeros(C,CF);
Phi_T = zeros(C+CN,CN+CF,C);

for c=1:C
    Phi_F(c,(c-1)*F+1:c*F) = 1; % A country owns domestic factors
    Phi_T(c,1:CN,c) = 1;
    Phi_T(C+(c-1)*N+1:C+c*N,1:CN,c) = 1; % HH collects tariff revenue when imported
end

data.Phi_F = Phi_F;
data.Phi_T = Phi_T;

data_org = data;
shock_org = shock;

% Initial share : Trade of interest
trade_of_interest = zeros(C+CN,CN+CF);
trade_of_interest(EU,(RUS-1)*N+1:RUS*N) = 1; % EU put iceberg trade cost on Russian consumption good
for e=1:length(EU)
    trade_of_interest(C+(EU(e)-1)*N+1:C+EU(e)*N,(RUS-1)*N+1:RUS*N) = 1; % EU put iceberg trade cost on Russain intermediate input
end
% trade_share_EU = sum(data_org.lambda_std(1:C+CN)'*(data_org.Omega_total(1:C+CN,C+1:end).*(trade_of_interest==1)))/sum(GNE_weights(EU))
% trade_share_RUS = sum(data_org.lambda_std(1:C+CN)'*(data_org.Omega_total(1:C+CN,C+1:end).*(trade_of_interest==1)))/GNE_weights(RUS)

%% For-loop to solve nonlinear system of equations
dlogW = zeros(C,ngrid); % (Main output) Change in real income across countries
dlogW_world = zeros(1,ngrid); % Change in real income of world
dlogR = zeros(C,ngrid); % Reallocation term 
dlogP_s = zeros(C,N,ngrid); % Change in sector-level prices

for i = 1:ngrid
%% 1. Input shocks
% dlogtau : change in iceberg trade cost in log change - Matrix of size (C+CN,CN+CF)
% dlogt : change in tariff in log change - Matrix of size (C+CN,CN+CF)
% (i,j) element of dlogtau: log change in iceberg trade cost when household/sector i buys from sector/factor j
% (i,j) element of dlogt: log change in tariff when household/sector i buys from sector/factor j

dlogt = zeros(C+CN,CN+CF);
dlogtau = zeros(C+CN,CN+CF);

% intensity_grid = (log(1+(i/ngrid)*intensity/100)-log(1+((i-1)/ngrid)*intensity/100));
intensity_grid = log(1+intensity/100)/ngrid;

if shock_index == 1 % universal change in iceberg trade cost

    dlogtau(:,1:CN) =intensity_grid*ones(C+CN,CN);
    for c=1:C
        dlogtau(c,(c-1)*N+1:c*N) = 0;
        dlogtau(C+(c-1)*N+1:C+c*N,(c-1)*N+1:c*N) = 0;
    end
    
elseif shock_index == 2 % universal change in tariff
    
%     dlogt(:,1:CN) = (log(1+(i/ngrid)*intensity/100)-log(1+((i-1)/ngrid)*intensity/100))*ones(C+CN,CN);
    dlogt(:,1:CN) =intensity_grid*ones(C+CN,CN);

    for c=1:C
        dlogt(c,(c-1)*N+1:c*N) = 0;
        dlogt(C+(c-1)*N+1:C+c*N,(c-1)*N+1:c*N) = 0;
    end


elseif shock_index == 3 % EU iceberg trade cost on Russia

    dlogtau(EU,(RUS-1)*N+1:RUS*N) = intensity_grid; % EU put iceberg trade cost on Russian consumption good
    for e=1:length(EU)
        dlogtau(C+(EU(e)-1)*N+1:C+EU(e)*N,(RUS-1)*N+1:RUS*N) = intensity_grid; % EU put iceberg trade cost on Russain intermediate input
    end

elseif shock_index == 4 % German tariff on Russian energy % transport

    dlogt(DEU,(RUS-1)*N+Energy) = intensity_grid; % German consumer put tariff on Russian energy
    dlogt(C+(DEU-1)*N+1:DEU*N,(RUS-1)*N+Energy) = intensity_grid; % German producers put tariff on Russain energy
    dlogt(DEU,(RUS-1)*N+Transport) = intensity_grid;
    dlogt(C+(DEU-1)*N+1:DEU*N,(RUS-1)*N+Transport) = intensity_grid; 
    
end
    
    dX = sum(data.Omega_total_tilde(1:C+CN,C+1:end).*(dlogt+dlogtau),2); % Matrix of size (C+CN,1)
    
    shock.dlogt = dlogt;
    shock.dlogtau = dlogtau;
    shock.dX = dX;

%% 2. Evaluate Allen Elasticities

[AES_C_Mat, AES_N_Mat, AES_F_Mat] = AES_func(data);

data.AES_C_Mat = AES_C_Mat;
data.AES_N_Mat = AES_N_Mat;
data.AES_F_Mat = AES_F_Mat;

%% 3. Set system in explicit form: [dlambda_F ; dlambda_F_star] = A*[dlambda_F ; dlambda_F_star] + B

[A_orig, B_orig, A, B] = Nested_CES_linear_final_rev(data, shock, factor_index);

dlambda_F_all = inv(eye(C+CF)-A)*B;

dlambda_F = dlambda_F_all(1:CF,1);
dlambda_F_star = dlambda_F_all(CF+1:CF+C,1);

%% 4. Calculate derivatives to derive change in real income across countries

[dlambda_result, dlambda_F_star_result, dchi_std, dlogP_Vec, dlogP_Mat, dOmega_total, dPsi_total, dOmega_total_tilde, dPsi_total_tilde] =...
        Nested_CES_linear_result_final(dlambda_F_all, data, shock);

%% 5. Store change in real income and reallocation term
dlogW_std = dchi_std./(data.chi_std) - dlogP_Vec;
dlogW(:,i) = dlogW_std(1:C,1);
dlogW_world(1,i) = chi_std(1:C,1)'*dlogW(:,i);

% 10% iceberg cost shock without initial tariff
lambda_Fc = zeros(C,CF);
for c=1:C
   lambda_Fc(c,(c-1)*F+1:c*F) = data.lambda_F((c-1)*F+1:c*F,1)'/data.chi_std(c,1);
end
lambda_FWc = data.Psi_total_tilde(1:C,C+CN+1:end); 

for c=1:C
   dlogR(c,i) = sum((lambda_Fc(c,:) - lambda_FWc(c,:)).*(dlambda_F./data.lambda_F)');
end
%% 6. Update variables

% Update tariff matrix
shock.init_t = (shock.init_t).*exp(dlogt);
init_t_mat = [zeros(C+CN,C), shock.init_t ; zeros(CF,L)];

% Update Omega_total_tilde (renormalize so that all rows add up to 1)
Omega_total_tilde_org = data.Omega_total_tilde;
data.Omega_total_tilde = data.Omega_total_tilde + dOmega_total_tilde;
data.Omega_total_tilde(data.Omega_total_tilde<0)=0;
dOmega_total_tilde = data.Omega_total_tilde - Omega_total_tilde_org;
data.Omega_total_tilde = bsxfun(@rdivide, data.Omega_total_tilde, sum(data.Omega_total_tilde,2));
data.Omega_total_tilde(isnan(data.Omega_total_tilde)) = 0;

% All other updated variables depend on Omega_total_tilde
data.Omega_total = data.Omega_total_tilde./init_t_mat;
data.Omega_total(isnan(data.Omega_total)) = 0;
data.Psi_total_tilde = inv(eye(L)-data.Omega_total_tilde);
data.Psi_total = inv(eye(L)-data.Omega_total);
data.chi_std = data.chi_std + dchi_std;
lambda_total = ((data.chi_std)'*data.Psi_total)';
data.lambda_CN = lambda_total(1:C+CN,1);
data.lambda_F = lambda_total(C+CN+1:end,1);
data.lambda_std = [data.lambda_CN ; data.lambda_F];

% Update other input shares to calculate sector-level prices P_s
% beta_s (b_ic) : How much HH from country c consumes sector i good
% beta_disagg (delta_im^0c) : How much HH from country c consumes sector i from country m
% Omega_s (omega_j^ic) : How much sector i from country c uses sector j
% Omega_disagg (delta_jm^ic) : How much sector i from country c uses sector j from country m

dbeta_temp = dOmega_total_tilde(1:C,C+1:C+CN);
dOmega_temp = dOmega_total_tilde(C+1:C+CN,C+1:end);

dOmega_temp2 = dOmega_total_tilde(C+1:C+CN,C+1:C+CN).*repmat(1./(1-data.alpha),1,CN);

dbeta_s = zeros(C,N);
dOmega_s = zeros(CN,N);

for c=1:C
   dbeta_s = dbeta_s + dbeta_temp(:,(c-1)*N+1:c*N);
   dOmega_s = dOmega_s + dOmega_temp2(:,(c-1)*N+1:c*N);
end

data.beta_s = data.beta_s + dbeta_s;
data.Omega_s = data.Omega_s + dOmega_s;

data.Omega_total_C = data.Omega_total_C + dbeta_temp;
data.Omega_total_N = data.Omega_total_N + dOmega_temp;

beta_disagg_update = cell(C,N);
beta_disagg_temp_update = (data.Omega_total_C)./repmat(data.beta_s, 1, C);

for c=1:C
    for j=1:N
        for k=1:C
            beta_disagg_update{c,j}(1,k) = beta_disagg_temp_update(c,(k-1)*N+j);
        end
    end
end


for c=1:C
    for j=1:N
        index = j:N:C*N;
        dlogP_s(c,j,i) = 1/2*(data.beta_disagg{c,j}+beta_disagg_update{c,j})*dlogP_Mat(c,index)';
        data.beta_disagg{c,j} = beta_disagg_update{c,j};
    end
end

% data.Omega_total_tilde = data.Omega_total_tilde + dOmega_total_tilde;
% data.Omega_total_tilde = bsxfun(@rdivide, data.Omega_total_tilde, sum(data.Omega_total_tilde,2));
% data.Omega_total_tilde(isnan(data.Omega_total_tilde)) = 0;
% data.Omega_total = data.Omega_total + dOmega_total;
% data.Psi_total_tilde = data.Psi_total_tilde + dPsi_total_tilde;
% data.Psi_total = data.Psi_total + dPsi_total;
% data.lambda_CN = data.lambda_CN + dlambda_result(1:C+CN,1);
% data.lambda_F = data.lambda_F + dlambda_result(C+CN+1:end,1);
% data.chi_std = data.chi_std + dchi_std;

end

%% Output
% dlogW_sum : Change in real income of each country from linearized system
% dlogR_sum : Reallocation term of each country from linearized system
dlogW_sum = sum(dlogW,2);
dlogR_sum = sum(dlogR,2);
dlogW_sum_world = GNE_weights'*dlogW_sum;

Country_list = [keep_c 35 42];
countries = {
    'AUS'
    'AUT'
    'BEL'
    'BGR'
    'BRA'
    'CAN'
    'CHN'
    'CYP'
    'CZE'
    'DEU'
    'DNK'
    'ESP'
    'EST'
    'FIN'
    'FRA'
    'GBR'
    'GRC'
    'HUN'
    'IDN'
    'IND'
    'IRL'
    'ITA'
    'JPN'
    'KOR'
    'LTU'
    'LUX'
    'LVA'
    'MEX'
    'MLT'
    'NLD'
    'POL'
    'PRT'
    'ROU'
    'RUS'
    'SVK'
    'SVN'
    'SWE'
    'TUR'
    'TWN'
    'USA'
    'ROW'
    'World'};
countries_reorder = cell(42,1);
countries_reorder(1:34) = countries(1:34);
countries_reorder(35) = countries(41);
countries_reorder(36:41) = countries(35:40);
countries_reorder(end) = countries(end);
Tab1 = table(Country_list', countries_reorder(Country_list'), [round(dlogW_sum,4) ; round(dlogW_sum_world,4)])
Volume_of_Trade = sum(data.lambda_std(1:C+CN)'*(data.Omega_total(1:C+CN,C+1:end).*(dlogtau~=0)))
Expenditure_change = [sum(data_org.lambda_std(1:C+CN)'*(data_org.Omega_total(1:C+CN,C+1:end).*(trade_of_interest~=0))), sum(data.lambda_std(1:C+CN)'*(data.Omega_total(1:C+CN,C+1:end).*(trade_of_interest~=0)))]

sectors = {
    'Agriculture, Hunting, Forestry and Fishing'
    'Mining and Quarrying'
    'Food, Beverages and Tobacco'
    'Textiles and Textile Products AND Leather, Leather and Footwear'
    'Wood and Products of Wood and Cork'
    'Pulp, Paper, Paper , Printing and Publishing'
    'Coke, Refined Petroleum and Nuclear Fuel'
    'Chemicals and Chemical Products AND Rubber and Plastics'
    'Other Non-Metallic Mineral'
    'Basic Metals and Fabricated Metal'
    'Machinery, Nec'
    'Electrical and Optical Equipment'
    'Transport Equipment'
    'Manufacturing, Nec; Recycling'
    'Electricity, Gas and Water Supply'
    'Construction'
    'Sale, Maintenance and Repair of Motor Vehicles and Motorcycles; Retail Sale of Fuel AND Wholesale Trade and Commission Trade, Except of Motor Vehicles and Motorcycles'
    'Retail Trade, Except of Motor Vehicles and Motorcycles; Repair of Household Goods'
    'Hotels and Restaurants'
    'Inland Transport'
    'Water Transport'
    'Air Transport'
    'Other Supporting and Auxiliary Transport Activities; Activities of Travel Agencies'
    'Post and Telecommunications'
    'Financial Intermediation'
    'Real Estate Activities'
    'Renting of M&Eq and Other Business Activities'
    'Public Admin and Defence; Compulsory Social Security'
    'Education'
    'Health and Social Work AND Other Community, Social and Personal Services AND Private Households with Employed Persons'
};
%% HH Expenditure quantity change 
% Output: logQ_s_change
% Row : coi (country of interest)
% Column : soi (sector of interest)
beta_s_initial = data_org.beta_s(coi,soi); % initial consumption expenditure share of coi countries on soi goods
beta_s_final = data.beta_s(coi,soi); % final consumption expenditure share of coi countries on soi goods

logP_s_initial = zeros(length(coi),length(soi)); % initial log price level (normalized to 1)
logP_s_final = sum(dlogP_s(coi,soi,:),3); % final log price level

logQ_s_change = (log(beta_s_final.*repmat(data.chi_std(coi),1,length(soi)))-logP_s_final) - (log(beta_s_initial.*repmat(data_org.chi_std(coi),1,length(soi)))-logP_s_initial);

Tab2 = table(EU', countries(EU'), round(logQ_s_change,4))
% disp('Column 1: Coke, Refined Petroleum and Nuclear Fuel')
% disp('Column 2: Inland Transport')

%% Change in world GDP to a 2nd order
dlogOmega = dOmega_total./Omega_total;
dloglambda = dlambda_result./lambda_std;
dlogOmega(isnan(dlogOmega)) = 0;
dloglambda(isnan(dloglambda)) = 0;

dlogX = dlogOmega + repmat(dloglambda,1,L) - repmat(dlogP_Vec',L,1);
dlogY_2nd = 1/2*lambda_CN'*sum(Omega_total(1:C+CN,C+1:end).*dlogX(1:C+CN,C+1:end).*dlogt,2);

%% Analysis on A,B

if AB_analysis == 1

dlambda_F_true = inv(eye(C+CF)-A)*B;
[V,D] = eig(A);

% I. Distribution of eigenvalue of A : 4 eigenvalues > 1
% hist(abs(diag(D)))

% II. Contribution on each eigenvalue of A : Always true!
n_artif = [1 5:17]; % Range : 5 - 20
c_artif = rand(1,length(n_artif));
B_artif = zeros(C+CF,1);
for i=1:length(n_artif)
    B_artif = B_artif + c_artif(i)*V(:,n_artif(i)); 
end
eig_artif = diag(D(n_artif,n_artif));

dlambda_F_artif1 = inv(eye(C+CF)-A)*B_artif;
dlambda_F_artif2 = zeros(C+CF,1);
for i=1:length(n_artif)
    dlambda_F_artif2 = dlambda_F_artif2 + c_artif(i)/(1-eig_artif(i))*V(:,n_artif(i));
end
[dlambda_F_artif1, dlambda_F_artif2]

% III. dlogW comparison : Not working
N_approx = 10;
dlambda_F_approx = cell(N_approx+1,1);
dlambda_F_approx{1,1} = B;
for i=2:N_approx
    dlambda_F_approx{i,1} = A*dlambda_F_approx{i-1,1} + B;
end

dlambda_F_approx{end,1} = dlambda_F_true;

% 1. Input shocks
shock_org = input_shock(data_org,shock_org,shock_index,1,ngrid,intensity);

% 2. Evaluate Allen Elasticities
[AES_C_Mat, AES_N_Mat, AES_F_Mat] = AES_func(data_org);

data_org.AES_C_Mat = AES_C_Mat;
data_org.AES_N_Mat = AES_N_Mat;
data_org.AES_F_Mat = AES_F_Mat;

dlogW_approx = zeros(C,N_approx+1);

for i=1:N_approx+1
% 4. Calculate derivatives to derive change in real income across countries
[dlambda_result, dlambda_F_star_result, dchi_std, dlogP_Vec, dlogP_Mat, dOmega_total, dPsi_total, dOmega_total_tilde, dPsi_total_tilde] =...
        Nested_CES_linear_result_final(dlambda_F_approx{i,1}, data_org, shock_org);

% 5. Store change in real income and reallocation term
dlogW_std = dchi_std./((chi_std + data_org.chi_std)/2) - dlogP_Vec;
dlogW_approx(:,i) = dlogW_std(1:C,1);
end

dlogW_approx

end

toc