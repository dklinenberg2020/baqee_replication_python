
function [Omega, beta, alpha_VA, alpha, trade_elast, GDP_weights] = IO_reorder(keep_c)
%% Reorder IO matrix for arbitrary countries
% Code by Sihwan Yang (PhD student in UCLA Economics), 01/20/2022

%% Input
% keep_c : vector of country numbers
% set of countries to keep, in the ascending order
% #35 should always be omitted 

% Case 1: All countries
% keep_c = (1:41); keep_c(35) = []; 

% Case 2: 3 countries - MEX (28), CHN (7), USA (41)
% keep_c = [28 7 41]

%% Import Data & Initialize a IO Matrix

load wiott2008.mat;
load trade_elast_2008.mat;
load wiodsea2008.mat;

N = 31;
N_orig = 31;
J = 41;
year = 2008;
years = 15;
temp=2;

%% Data cleaning
%wiod = reshape(wiodprocess, J*N+3,  years, J*(N+3)+2);
%IO = reshape(IO, J*N+3, J*(N+3)+2);
IO = wiott2008;
%IO = mean(wiod, 2);
%IO(:, 1:2) = [];
IO(IO < 0) = 0;

%Ind = IO(:,1); %store industry codes

% N+4 -- Private Households with Employed Persons
% N+3 -- Other Community, Social and Personal Services
% N -- Public Admin and Defence; Compulsory Social Security
% 4 -- Textiles and Textile Products
% 5 -- Leather, Leather and Footwear
% 19 -- Sale, Maintenance and Repair of Motor Vehicles and Motorcycles; Retail Sale of Fuel
% 20 -- Wholesale Trade and Commission Trade, Except of Motor Vehicles and Motorcycles

IO(N:(N+4):(J-1)*(N+4)+N, :) = IO(N:(N+4):(J-1)*(N+4)+N, :)+IO(N+3:(N+4):(J-1)*(N+4)+N+3, :)+IO(N+4:(N+4):(J-1)*(N+4)+N+4, :);
IO(:, N:(N+9):(J-1)*(N+9)+N) = IO(:, N:(N+9):(J-1)*(N+9)+N)+IO(:, N+3:(N+9):(J-1)*(N+9)+N+3)+IO(:, N+4:(N+9):(J-1)*(N+9)+N+4);
IO(4:(N+4):(J-1)*(N+4)+4, :) = IO(4:(N+4):(J-1)*(N+4)+4, :)+IO(5:(N+4):(J-1)*(N+4)+5, :);
IO(:, 4:(N+9):(J-1)*(N+9)+4) = IO(:, 4:(N+9):(J-1)*(N+9)+4)+IO(:, 5:(N+9):(J-1)*(N+9)+5);
IO(19:(N+4):(J-1)*(N+4)+19, :) = IO(19:(N+4):(J-1)*(N+4)+19, :)+IO(20:(N+4):(J-1)*(N+4)+20, :);
IO(:, 19:(N+9):(J-1)*(N+9)+19) = IO(:, 19:(N+9):(J-1)*(N+9)+19)+IO(:, 20:(N+9):(J-1)*(N+9)+20);
IO(9:(N+4):(J-1)*(N+4)+9, :) = IO(9:(N+4):(J-1)*(N+4)+9, :)+IO(8:(N+4):(J-1)*(N+4)+8, :);
IO(:, 9:(N+9):(J-1)*(N+9)+9) = IO(:, 9:(N+9):(J-1)*(N+9)+9)+IO(:, 8:(N+9):(J-1)*(N+9)+8);

IO([N+3:(N+4):(J-1)*(N+4)+N+3, N+4:(N+4):(J-1)*(N+4)+N+4, 5:(N+4):(J-1)*(N+4)+5, 20:(N+4):(J-1)*(N+4)+20, 8:(N+4):(J-1)*(N+4)+8], :) = [];
IO(:, [N+3:(N+9):(J-1)*(N+9)+N+3, N+4:(N+9):(J-1)*(N+9)+N+4, 5:(N+9):(J-1)*(N+9)+5, 20:(N+9):(J-1)*(N+9)+20, 8:(N+9):(J-1)*(N+9)+8]) = [];

IO = IO';

N = 30;
N_orig = 30;

% should drop taxes less subsidies here (& gross output)
H_cons_h = IO(N+1:N+5:(J-1)*(N+5)+N+1, 1:J*N)'+IO(N+4:N+5:(J-1)*(N+5)+N+4, 1:J*N)'+IO(N+5:N+5:(J-1)*(N+5)+N+5, 1:J*N)';
NP_cons_h = IO(N+2:N+5:(J-1)*(N+5)+N+2, 1:J*N)';
Gov_cons_h = IO(N+3:N+5:(J-1)*(N+5)+N+3, 1:J*N)'; 
Tot_cons_h = (H_cons_h + NP_cons_h + Gov_cons_h);

IO([N+5:N+5:(J-1)*(N+5)+N+5, N+4:N+5:(J-1)*(N+5)+N+4, N+3:N+5:(J-1)*(N+5)+N+3, N+2:N+5:(J-1)*(N+5)+N+2, N+1:N+5:(J-1)*(N+5)+N+1], :) = [];
IO(1:J*N, end) = sum(IO(:, 1:J*N), 1)'+sum(Tot_cons_h, 2) - sum(IO(1:J*N, 1:J*N),  2);

VA_all = IO(1:N*J, end);
alpha = VA_all;
alpha(isnan(alpha)) = 0;
alpha = max(0, alpha);
IO_normed = IO;
IO_normed(isnan(IO_normed)) = 0;

Omega = IO_normed(1:J*N, 1:J*N);
Omega(isnan(Omega)) = 0;
Omega(isinf(Omega)) = 1;

%Total consumption
Tot_cons_h = (H_cons_h + NP_cons_h + Gov_cons_h);

beta = Tot_cons_h(1:J*N, :);
sea = wiodsea2008;
sea(isnan(sea)) = 0;

N = 31;
N_orig = 31;
VA = horzcat(sea(1:(J-1)*(N+4), 1), sea((J-1)*(N+4)+1:2*(J-1)*(N+4), 1).*sea(3*(J-1)*(N+4)+1:4*(J-1)*(N+4), 1), sea((J-1)*(N+4)+1:2*(J-1)*(N+4), 1).*sea(2*(J-1)*(N+4)+1:3*(J-1)*(N+4), 1), sea(1*(J-1)*(N+4)+1:2*(J-1)*(N+4), 1).*sea(4*(J-1)*(N+4)+1:5*(J-1)*(N+4), 1));

% ROW -- 35
VA_new = zeros(J*(N+4), 4);
VA_new(1:34*(N+4), :) = VA(1:34*(N+4), :);
VA_new(34*(N+4)+1:35*(N+4), :) = 250*ones(size(VA_new(33*(N+4)+1:34*(N+4), :)));
VA_new(35*(N+4)+1:end, :) = VA(34*(N+4)+1:end, :);
VA = VA_new;

VA(4:(N+4):(J-1)*(N+4)+4, :) = VA(4:(N+4):(J-1)*(N+4)+4, :)+VA(5:(N+4):(J-1)*(N+4)+5, :);
VA(19:(N+4):(J-1)*(N+4)+19, :) = VA(19:(N+4):(J-1)*(N+4)+19, :)+VA(20:(N+4):(J-1)*(N+4)+20, :);
VA(9:(N+4):(J-1)*(N+4)+9, :) = VA(9:(N+4):(J-1)*(N+4)+9, :)+VA(8:(N+4):(J-1)*(N+4)+8, :);
VA([N+3:(N+4):(J-1)*(N+4)+N+3, N+4:(N+4):(J-1)*(N+4)+N+4, 5:(N+4):(J-1)*(N+4)+5,   20:(N+4):(J-1)*(N+4)+20, 8:(N+4):(J-1)*(N+4)+8], :) = [];

VA = max(VA, 0);
alpha_VA = VA;
alpha_VA(isnan(alpha_VA)) = 0;
alpha_VA = max(0, alpha_VA);

N = 30;
N_orig = 30;

for i = 1:J*N_orig
    vect = sum(alpha_VA, 2);
    if vect(i) == 0
        alpha(i) = 0;
    end
end

%% Transforming the dataset
% IO_new : Dollar value of IO matrix within 'keep_c' countries
% beta_new : Dollar value of Final consumption within 'keep_c' countries
% alpha_VA_new : Dollar value of Factor usage of 'keep_c' countries
% alpha_new : Dollar value of Value-added usage of 'keep_c' countries

M = 4; 
J_new = length(keep_c);
Omega_new = zeros(N*J_new, N*J_new);
beta_new = zeros(J_new, N*J_new)';
alpha_VA_new = zeros(N*J_new, M);
alpha_new = zeros(N*J_new, 1);
Sum_tot = sum(Tot_cons_h, 1)';
sum_cons_new = zeros(J_new, 1);

for i = 1:J
    for j = 1:J
        if (ismember(i, keep_c) == 1) && (ismember(j, keep_c) == 1)
            k = find(keep_c == i);
            l = find(keep_c == j);
            sum_cons_new(k, 1) = Sum_tot(i, 1);
            IO_new((k-1)*N+1:k*N, (l-1)*N+1:l*N) = IO((i-1)*N+1:i*N, (j-1)*N+1:j*N);
            beta_new((k-1)*N+1:k*N, l) = beta((i-1)*N+1:i*N, j);
            alpha_VA_new((k-1)*N+1:k*N, :) = alpha_VA((i-1)*N+1:i*N, :);
            alpha_new((k-1)*N+1:k*N, 1) = alpha((i-1)*N+1:i*N, 1);
        end
    end
end

%% Transforming the dataset
% IO_hh : Dollar value of IO matrix within 'keep_c' countries
% IO_ff : Dollar value of IO matrix within ROW countries
% IO_hf : Dollar value of IO matrix that ROW buys from 'keep_c' countries
% IO_fh : Dollar value of IO matrix that 'keep_c' countries buys from ROW

J_new = length(keep_c);

IO_hh = IO_new;
IO_ff = zeros(N, N);
IO_ff_mod = IO(1:J*N, 1:J*N);
for j = 1:41
    if ismember(j, keep_c)
        IO_ff_mod(1:J*N, N*(j-1)+1:N*j) = 0;
        IO_ff_mod(N*(j-1)+1:N*j, 1:J*N) = 0;
    end
end

for k = 1:J
    for l = 1:J
        IO_ff = IO_ff+IO_ff_mod(N*(k-1)+1:N*k, N*(l-1)+1:N*l);
    end
end

IO_mod_fh = IO(1:J*N, 1:J*N);
% IO_fh -- home usage of the foreign goods
IO_fh = zeros(J_new*N, N);
IO_mod_hf = IO(1:J*N, 1:J*N);
% IO_hf -- foreign usage of the home goods 
IO_hf = zeros(N, J_new*N);

for j = 1:41
    if ismember(j, keep_c)
        IO_mod_fh(N*(j-1)+1:N*j, N*(j-1)+1:N*j) = 0;
        IO_mod_fh(1:J*N, N*(j-1)+1:N*j) = 0;
        IO_mod_hf(N*(j-1)+1:N*j, N*(j-1)+1:N*j) = 0;
        IO_mod_hf(N*(j-1)+1:N*j, 1:J*N) = 0;
    else
        IO_mod_fh(N*(j-1)+1:N*j, 1:J*N) = 0;
        IO_mod_hf(1:J*N, N*(j-1)+1:N*j) = 0;
    end
end
            
for k = 1:length(keep_c)
    for l= 1:41
        j = keep_c(k);
        
        IO_fh(N*(k-1)+1:N*k, 1:N) = IO_fh(N*(k-1)+1:N*k, 1:N)+IO_mod_fh(N*(j-1)+1:N*j, (l-1)*N+1:l*N);
        IO_hf(1:N, N*(k-1)+1:N*k) = IO_hf(1:N, N*(k-1)+1:N*k)+IO_mod_hf((l-1)*N+1:l*N, N*(j-1)+1:N*j);
    end
end

%% Final consumption
% Home consumption
beta_hh = beta_new;
beta_fh = zeros(N, J_new);
beta_ff = zeros(N, 1);
beta_hf = zeros(N*J_new, 1);
    
    for j = 1:41
        for l =1:length(keep_c)
        if not(ismember(j, keep_c))
            k = keep_c(l);
        beta_fh(1:N, l) = beta_fh(1:N, l)+beta(N*(j-1)+1:N*j, k);
        beta_hf(N*(l-1)+1:N*l, 1) = beta_hf(N*(l-1)+1:N*l, 1)+beta(N*(k-1)+1:N*k, j);
        end
        
        end
        for k = 1:41
            if not(ismember(k, keep_c)) && not(ismember(j, keep_c))
            beta_ff(1:N, 1) = beta_ff(1:N, 1)+beta(N*(j-1)+1:N*j, k);
            end
        end
    end
    
%% Value added
% Home value added 
alpha_h = alpha_new;
% Foreign value added
alpha_f = zeros(N, 1);

%% Value Added Factors
alpha_VA_h = alpha_VA_new;
alpha_VA_f = zeros(N, M); 

for k = 1:41
    if not(ismember(k, keep_c))
    alpha_f = alpha_f + alpha(N*(k-1)+1:N*k, :);
    alpha_VA_f = alpha_VA_f + alpha_VA(N*(k-1)+1:N*k, :);
    end
end

%% Output
IO_new = [IO_hh, IO_fh; IO_hf, IO_ff];
alpha_new = [alpha_h; alpha_f];
alpha_VA_new = [alpha_VA_h; alpha_VA_f];
beta_new = [beta_hh, beta_hf; beta_fh, beta_ff];

J = J_new;
alpha_VA = bsxfun(@ldivide, sum(alpha_VA_new, 2), alpha_VA_new);

Omega = bsxfun(@ldivide, sum(IO_new, 2), IO_new);
beta = bsxfun(@ldivide, sum(beta_new, 1), beta_new);

alpha = alpha_new./(alpha_new+sum(IO_new, 2));

vect1 = sum(IO_new, 2)+alpha_new;
vect2 = sum(beta_new, 2)+sum(IO_new, 1)';

shares_1 = sum((beta')*inv(eye((J+1)*N) - diag(1 - alpha)*Omega), 1)';
matr1 = (beta')*inv(eye((J+1)*N) - diag(1 - alpha)*Omega);

GDP_weights = sum(beta_new, 1)';
GDP_weights = GDP_weights./sum(GDP_weights);

Omega(isnan(Omega)) = 0;
Omega(isinf(Omega)) = 1;
beta(isnan(beta)) = 0;
beta(isinf(beta)) = 1;
%% Set Elasticities
J = J_new+1; 
M = 4;

epsilon = 0.5*ones(J*N, 1); % elasticity across VA and Intermediates
theta = 0.2*ones(J*N, 1); % elasticity across intermediates
eta = 0.5*ones(J*N, 1); % between VA types
trade_elast = 3*ones(N, 1) ;
trade_elast = [trade_elast_2008(1:8); trade_elast_2008(10:end)];

end