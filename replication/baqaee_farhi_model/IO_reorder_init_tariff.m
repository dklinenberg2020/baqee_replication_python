
function [Omega, beta, alpha_VA, alpha, trade_elast, GDP_weights, Tariff_matrix_new, Tariffs_cons_matrix_new] = IO_reorder_init_tariff(keep_c)
%% Reorder IO matrix for arbitrary countries
% Code by Sihwan Yang (PhD student in UCLA Economics), 01/24/2022

%% Input
% keep_c : vector of country numbers
% set of countries to keep, in the ascending order
% #35 should always be omitted 

% Case 1: All countries
% keep_c = (1:41); keep_c(35) = []; 

% Case 2: 3 countries - MEX (28), CHN (7), USA (41)
% keep_c = [28 7 41]

%% Import Data & Initialize a IO Matrix
Initial_tariffs_mode = 2;

year = 2008;

trade_flow_file_name = strcat('wiot',  num2str(year), '_row_apr12.xlsx');

load('trade_elast_2008.mat');

load('WIOD_SEA_14.mat')
sea_data = WIODSEA14;
sea_data = sea_data(:, year-1994);
sea_data(isnan(sea_data))=0;

load('wiot2008_row_apr12.mat')
trade_flow_data = trade_flow_data(5:end, 5:end);

load('ahs_all.mat')
tariff_data = ahs_all;
tariff_data(find(tariff_data(:, 4) ~=year), :) = [];
tariff_data(isnan(tariff_data)) = 0;

vect = [];
for i= 1:length(tariff_data)
    if tariff_data(i, 4)== 0 || isnan(tariff_data(i, 4))
        vect = vertcat(vect, i);
    end 
end 

tariff_data(vect, :) = [];

N = 31;
N_orig = 31;
J = 41;
year = 2008;
years = 15;
temp=2;

wiot_data = zeros(J*(N+4)+1, J*(N+9));

list1 = [1:1:41];
list2 = [1:1:34, 41, 35:1:40];

for j = 1:J
    data = horzcat(trade_flow_data(1:J*(N+4), (list2(j)-1)*(N+4)+1:(list2(j))*(N+4)), trade_flow_data(1:J*(N+4), J*(N+4)+(list2(j)-1)*(5)+1:J*(N+4)+(list2(j))*(5)));
    data_va = horzcat(trade_flow_data(1441,(list2(j)-1)*(N+4)+1:(list2(j))*(N+4)), zeros(1, 5));
    wiot_data(1:J*(N+4)+1, (list1(j)-1)*(N+9)+1:(list1(j))*(N+9)) = vertcat(data, data_va);
end

row_data = wiot_data((41-1)*(N+4)+1:(41)*(N+4), :);
wiot_data((36-1)*(N+4)+1:(41)*(N+4), :) = wiot_data((35-1)*(N+4)+1:(40)*(N+4), :);
wiot_data((35-1)*(N+4)+1:(35)*(N+4), :) = row_data;

%% Data cleaning
%wiod = reshape(wiodprocess, J*N+3,  years, J*(N+3)+2);
%IO = reshape(IO, J*N+3, J*(N+3)+2);
IO = wiot_data;
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


IO_normed = IO;
IO_normed(isnan(IO_normed)) = 0;

Omega = IO_normed(1:J*N, 1:J*N);
Omega(isnan(Omega)) = 0;
Omega(isinf(Omega)) = 1;

%Total consumption
Tot_cons_h = (H_cons_h + NP_cons_h + Gov_cons_h);

beta = Tot_cons_h(1:J*N, :);

N = 31;
N_orig = 31;
% order: capital, low skill, medium skill, high skill
VA = horzcat(sea_data(1:(J-1)*(N+4), 1), sea_data((J-1)*(N+4)+1:2*(J-1)*(N+4), 1).*sea_data(3*(J-1)*(N+4)+1:4*(J-1)*(N+4), 1), sea_data((J-1)*(N+4)+1:2*(J-1)*(N+4), 1).*sea_data(4*(J-1)*(N+4)+1:5*(J-1)*(N+4), 1), sea_data(1*(J-1)*(N+4)+1:2*(J-1)*(N+4), 1).*sea_data(2*(J-1)*(N+4)+1:3*(J-1)*(N+4), 1));

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

%% Adding tariffs to the IO data

tau_int_old = ones(J*N_orig, J*N_orig);
tau_cons_old = ones(J, J*N_orig);

Tau_c_old = ones(J*N_orig+J, J*N_orig+J, J);

tariff_rates = tariff_data(:, 5);

c_origin = tariff_data(:, 1);

c_origin(c_origin == 41) = 42;
c_origin(c_origin == 40) = 41;
c_origin(c_origin == 39) = 40;
c_origin(c_origin == 38) = 39;
c_origin(c_origin == 37) = 38;
c_origin(c_origin == 36) = 37;
c_origin(c_origin == 35) = 36;
c_origin(c_origin == 42) = 35;

ind_origin = tariff_data(:, 3);
origin = N_orig*(c_origin-1)+ind_origin;
c_destination = tariff_data(:, 2);
c_destination(c_destination == 41) = 42;
c_destination(c_destination == 40) = 41;
c_destination(c_destination == 39) = 40;
c_destination(c_destination == 38) = 39;
c_destination(c_destination == 37) = 38;
c_destination(c_destination == 36) = 37;
c_destination(c_destination == 35) = 36;
c_destination(c_destination == 42) = 35;

destination = zeros(length(c_destination), N_orig);
for i = 1:length(c_destination)
destination(i, :) = [N_orig*(c_destination(i)-1)+1:N_orig*(c_destination(i))];
end
destination_c = [c_destination, J+destination];
rates = tariff_rates/100+1;

if Initial_tariffs_mode == 1
    rates = ones(size(rates));
end

for i = 1:length(origin)
    tau_cons_old(c_destination(i), origin(i)) = tau_cons_old(c_destination(i), origin(i))*rates(i);
    for j = 1:N_orig
        
        tau_int_old(destination(i, j), origin(i)) = tau_int_old(destination(i, j), origin(i))*rates(i);
        
    end 
end

for i = 1:length(origin)
    for j = 1:N_orig+1
        Tau_c_old(J+origin(i), destination_c(i, j), c_destination(i)) = Tau_c_old(J+origin(i), destination_c(i, j), c_destination(i))*rates(i);       
    end 
end

%% Value addded

VA_all = sum(IO(:, 1:N*J), 1)'+sum(Tot_cons_h, 2) - sum(IO(1:N*J, 1:N*J), 2) - sum(IO(1:N*J, 1:N*J).*(tau_int_old-1), 2);
alpha = VA_all;
alpha(isnan(alpha)) = 0;
alpha = max(0, alpha);

for i = 1:J*N_orig
    vect = sum(alpha_VA, 2);
    if vect(i) == 0
        alpha(i) = 0;
    end
end

%% Transforming the dataset


M = 4; 
J_new = length(keep_c);
Omega_new = zeros(N*J_new, N*J_new);
beta_new = zeros(J_new, N*J_new)';
alpha_VA_new = zeros(N*J_new, M);
alpha_new = zeros(N*J_new, 1);
Sum_tot = sum(Tot_cons_h, 1)';
sum_cons_new = zeros(J_new, 1);
Tariffs_new = zeros(N*J_new, N*J_new);
Tariffs_cons_new = zeros(J_new, N*J_new)';
Tau_c_new = ones(J_new*N_orig+J_new, J_new*N_orig+J_new, J);
tau_cons_old_ = tau_cons_old';

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
            Tariffs_new((k-1)*N+1:k*N, (l-1)*N+1:l*N) = tau_int_old((i-1)*N+1:i*N, (j-1)*N+1:j*N);
            Tariffs_cons_new((k-1)*N+1:k*N, l) = tau_cons_old_((i-1)*N+1:i*N, j);
            Tau_c_new(J_new+(l-1)*N+1:J_new+l*N, J_new+(k-1)*N+1:J_new+k*N, :) = Tau_c_old(J+(j-1)*N+1:J+j*N, J+(i-1)*N+1:J+i*N, :)-1;
            Tau_c_new(J_new+(k-1)*N+1:J_new+k*N, l, :) = Tau_c_old(J+(i-1)*N+1:J+i*N, j, :)-1;
        end
    end
end

%% Transforming the dataset
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
    
%% Transforming the tariff matrices
Tariffs_hh = Tariffs_new.*IO_new;
Tau_c_hh = Tau_c_new(J_new+1:end, J_new+1:end, :).*repmat(IO_new', 1, 1, J);
Tariffs_ff = zeros(N, N);
Tau_c_ff = zeros(N, N, J);
Tariffs_ff_mod = tau_int_old(1:J*N, 1:J*N).*IO(1:J*N, 1:J*N);
Tau_c_ff_mod = (Tau_c_old(J+1:end, J+1:end, :)-1).*repmat(IO(1:J*N, 1:J*N)', 1, 1, J);
for j = 1:41
    if ismember(j, keep_c)
        Tariffs_ff_mod(1:J*N, N*(j-1)+1:N*j) = 0;
        Tariffs_ff_mod(N*(j-1)+1:N*j, 1:J*N) = 0;
        Tau_c_ff_mod(N*(j-1)+1:N*j, 1:J*N, :) = 0;
        Tau_c_ff_mod(1:J*N, N*(j-1)+1:N*j, :) = 0;
    end
end

for k = 1:J
    for l = 1:J
        Tariffs_ff = Tariffs_ff+Tariffs_ff_mod(N*(k-1)+1:N*k, N*(l-1)+1:N*l);
        Tau_c_ff = Tau_c_ff+Tau_c_ff_mod(N*(l-1)+1:N*l, N*(k-1)+1:N*k, :);
    end
end

Tariffs_mod_fh = tau_int_old(1:J*N, 1:J*N).*IO(1:J*N, 1:J*N);
Tau_c_mod_fh = (Tau_c_old(J+1:end, J+1:end, :)-1).*repmat(IO(1:J*N, 1:J*N)', 1, 1, J);
% IO_fh -- home usage of the foreign goods
Tariffs_fh = zeros(J_new*N, N);
Tau_c_fh = zeros(N, J_new*N, J);
Tariffs_mod_hf = tau_int_old(1:J*N, 1:J*N).*IO(1:J*N, 1:J*N);
Tau_c_mod_hf = (Tau_c_old(J+1:end, J+1:end, :)-1).*repmat(IO(1:J*N, 1:J*N)', 1, 1, J);
% IO_hf -- foreign usage of the home goods 
Tariffs_hf = zeros(N, J_new*N);
Tau_c_hf = zeros(J_new*N, N, J);

for j = 1:41
    if ismember(j, keep_c)
        Tariffs_mod_fh(N*(j-1)+1:N*j, N*(j-1)+1:N*j) = 0;
        Tariffs_mod_fh(1:J*N, N*(j-1)+1:N*j) = 0;
        Tariffs_mod_hf(N*(j-1)+1:N*j, N*(j-1)+1:N*j) = 0;
        Tariffs_mod_hf(N*(j-1)+1:N*j, 1:J*N) = 0;
        
        Tau_c_mod_fh(N*(j-1)+1:N*j, N*(j-1)+1:N*j, :) = 0;
        Tau_c_mod_fh(N*(j-1)+1:N*j, 1:J*N, :) = 0;
        Tau_c_mod_hf(N*(j-1)+1:N*j, N*(j-1)+1:N*j, :) = 0;
        Tau_c_mod_hf(1:J*N, N*(j-1)+1:N*j, :) = 0;
    else
        Tariffs_mod_fh(N*(j-1)+1:N*j, 1:J*N) = 0;
        Tariffs_mod_hf(1:J*N, N*(j-1)+1:N*j) = 0;
        
        Tau_c_mod_fh(1:J*N, N*(j-1)+1:N*j, :) = 0;
        Tau_c_mod_hf(N*(j-1)+1:N*j, 1:J*N, :) = 0;
    end
end
            
for k = 1:length(keep_c)
    for l= 1:41
        j = keep_c(k);
        
        Tariffs_fh(N*(k-1)+1:N*k, 1:N) = Tariffs_fh(N*(k-1)+1:N*k, 1:N)+Tariffs_mod_fh(N*(j-1)+1:N*j, (l-1)*N+1:l*N);
        Tariffs_hf(1:N, N*(k-1)+1:N*k) = Tariffs_hf(1:N, N*(k-1)+1:N*k)+Tariffs_mod_hf((l-1)*N+1:l*N, N*(j-1)+1:N*j);
        Tau_c_fh(1:N, N*(k-1)+1:N*k, :) = Tau_c_fh(1:N, N*(k-1)+1:N*k, :)+Tau_c_mod_fh((l-1)*N+1:l*N, N*(j-1)+1:N*j, :);
        Tau_c_hf(N*(k-1)+1:N*k, 1:N, :) = Tau_c_hf(N*(k-1)+1:N*k, 1:N, :)+Tau_c_mod_hf(N*(j-1)+1:N*j, (l-1)*N+1:l*N, :);
    end
end    

%% Final consumption
% Home consumption
Tariffs_cons_hh = beta_new.*Tariffs_cons_new;
Tariffs_cons_fh = zeros(N, J_new);
Tariffs_cons_ff = zeros(N, 1);
Tariffs_cons_hf = zeros(N*J_new, 1);

Tau_c_cons_hh = Tau_c_new(J_new+1:end, 1:J_new, :).*repmat(beta_new, 1, 1, J);
Tau_c_cons_fh = zeros(N, J_new, J);
Tau_c_cons_ff = zeros(N, 1, J);
Tau_c_cons_hf = zeros(N*J_new, 1, J);
    
    for j = 1:41
        for l =1:length(keep_c)
        if not(ismember(j, keep_c))
            k = keep_c(l);
        Tariffs_cons_fh(1:N, l) = Tariffs_cons_fh(1:N, l)+beta(N*(j-1)+1:N*j, k).*tau_cons_old_(N*(j-1)+1:N*j, k);
        Tariffs_cons_hf(N*(l-1)+1:N*l, 1) = Tariffs_cons_hf(N*(l-1)+1:N*l, 1)+beta(N*(k-1)+1:N*k, j).*tau_cons_old_(N*(k-1)+1:N*k, j);
        
        Tau_c_cons_fh(1:N, l, :) = Tau_c_cons_fh(1:N, l, :)+repmat(beta(N*(j-1)+1:N*j, k), 1, 1, J).*(Tau_c_old(J+N*(j-1)+1:J+N*j, k, :)-1);
        Tau_c_cons_hf(N*(l-1)+1:N*l, 1, :) = Tau_c_cons_hf(N*(l-1)+1:N*l, 1, :)+repmat(beta(N*(k-1)+1:N*k, j), 1, 1, J).*(Tau_c_old(J+N*(k-1)+1:J+N*k, j, :)-1);
        end
        
        end
        for k = 1:41
            if not(ismember(k, keep_c)) && not(ismember(j, keep_c))
            Tariffs_cons_ff(1:N, 1) = Tariffs_cons_ff(1:N, 1)+beta(N*(j-1)+1:N*j, k).*tau_cons_old_(N*(j-1)+1:N*j, k);
            Tau_c_cons_ff(1:N, 1, :) = Tau_c_cons_ff(1:N, 1, :)+repmat(beta(N*(j-1)+1:N*j, k), 1, 1, J).*(Tau_c_old(J+N*(j-1)+1:J+N*j, k, :)-1);
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

IO_new = [IO_hh, IO_fh; IO_hf, IO_ff];
alpha_new = [alpha_h; alpha_f];
alpha_VA_new = [alpha_VA_h; alpha_VA_f];
beta_new = [beta_hh, beta_hf; beta_fh, beta_ff];
alpha_VA = bsxfun(@ldivide, sum(alpha_VA_new, 2), alpha_VA_new);
alpha_VA = bsxfun(@ldivide, sum(alpha_VA, 2), alpha_VA);
alpha_VA(isnan(alpha_VA)) = 0;
alpha_VA(find(sum(alpha_VA, 2)==0), :) = 1/M;
alpha_VA(isinf(alpha_VA)) = 1;
Tariff_matrix_new = [Tariffs_hh, Tariffs_fh; Tariffs_hf, Tariffs_ff]./IO_new;
Tariff_matrix_new(isnan(Tariff_matrix_new)) = 1;
Tariffs_cons_matrix_new = [Tariffs_cons_hh, Tariffs_cons_hf; Tariffs_cons_fh, Tariffs_cons_ff]./beta_new;
Tariffs_cons_matrix_new(isnan(Tariffs_cons_matrix_new)) = 1;
Tau_c_new = [ones(J_new+1, J_new+1+(J_new+1)*N, J);[Tau_c_cons_hh, Tau_c_cons_hf; Tau_c_cons_fh, Tau_c_cons_ff], [Tau_c_hh, Tau_c_hf; Tau_c_fh, Tau_c_ff]];
Tau_c_new_h = Tau_c_new(:, :, keep_c);
Tau_c_new_f = (sum(Tau_c_new, 3) - sum(Tau_c_new(:, :, keep_c), 3));
Tau_c_new = cat(3, Tau_c_new_h, Tau_c_new_f);

Tau_c_new(J_new+2:end, J_new+2:end, :) = Tau_c_new(J_new+2:end, J_new+2:end, :)./repmat(IO_new', 1, 1, J_new+1)+1;

Tau_c_new(J_new+2:end, 1:J_new+1, :) = Tau_c_new(J_new+2:end, 1:J_new+1, :)./repmat(beta_new,1, 1, J_new+1)+1;

%Tau_c_new(Tau_c_new ==0) = 1;

Tau_c_new(1:J_new+1, :, end) = 1;
Tau_c_new(isnan(Tau_c_new)) = 1;

Sectoral_tariffs = (sum(Tariff_matrix_new.*IO_new, 1)'+sum(Tariffs_cons_matrix_new.*beta_new, 2))./(sum(IO_new, 1)'+sum(beta_new, 2)); 

%% Outputs
J = J_new;

alpha = alpha_new./(alpha_new+sum(IO_new.*Tariff_matrix_new, 2));

Omega = bsxfun(@ldivide, sum(IO_new.*Tariff_matrix_new, 2), IO_new.*Tariff_matrix_new);
beta = bsxfun(@ldivide, sum(beta_new.*Tariffs_cons_matrix_new, 1), beta_new.*Tariffs_cons_matrix_new);

vect1 = sum(IO_new, 2)+alpha_new;
vect2 = sum(beta_new, 2)+sum(IO_new, 1)';

shares_1 = sum((beta')*inv(eye((J+1)*N) - diag((1 - alpha)./Sectoral_tariffs)*Omega), 1)';
matr1 = (beta')*inv(eye((J+1)*N) - diag((1 - alpha)./Sectoral_tariffs)*Omega);

GDP_weights = sum(beta_new.*Tariffs_cons_matrix_new, 1)';
GDP_weights = GDP_weights./sum(GDP_weights);

Omega(isnan(Omega)) = 0;
Omega(isinf(Omega)) = 1;
beta(isnan(beta)) = 0;
beta(isinf(beta)) = 1;
alpha(isnan(alpha)) = 0;
alpha(isinf(alpha)) = 1;

Sectoral_tariffs(isinf(Sectoral_tariffs)) = 1;
Sectoral_tariffs(isnan(Sectoral_tariffs)) = 1;

%% Set Elasticities
J = J_new+1; 
M = 4;

epsilon = 0.5*ones(J*N, 1); % elasticity across VA and Intermediates
theta = 0.2*ones(J*N, 1); % elasticity across intermediates
eta = 0.5*ones(J*N, 1); % between VA types
trade_elast = 3*ones(N, 1) ;
trade_elast = [trade_elast_2008(1:8); trade_elast_2008(10:end)];
end