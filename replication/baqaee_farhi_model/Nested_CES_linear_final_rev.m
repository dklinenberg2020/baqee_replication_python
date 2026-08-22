
function [A_orig, B_orig, A_result, B_result] = Nested_CES_linear_final_rev(data, shock, factor_index)
%% Implement Linearization model of Nested CES
% Code by Sihwan Yang (PhD student in UCLA Economics), 02/08/2022

%% Replication Code Description
% 1. main_load_data: First part of main code. Calculate inputs from data
% 2. main_dlogW: Second part of main code. Solve system of linear equations to derive change in welfare
% 3. AES_func: Calculate Allen Elasticity of Substitution
% **4. Nested_CES_linear_final: Solve system of linear equations by inverting the system
% 5. Nested_CES_linear_result_final: Calculate derivatives that are used to iterate or derive welfare change

%% Input

C = data.C;
N = data.N;
F = data.F;
CN = data.CN;
CF = data.CF;
L = data.L;

Omega_total = data.Omega_total;
Psi_total = data.Psi_total;
Omega_total_tilde = data.Omega_total_tilde;
Psi_total_tilde = data.Psi_total_tilde;
lambda_CN = data.lambda_CN;
lambda_F = data.lambda_F;

Phi_F = data.Phi_F;
Phi_T = data.Phi_T;

AES_C_Mat = data.AES_C_Mat;
AES_N_Mat = data.AES_N_Mat;
AES_F_Mat = data.AES_F_Mat;

dlogt = shock.dlogt;
dlogtau = shock.dlogtau;
init_t = shock.init_t;
dX = shock.dX;

%% dPsi_total = f(dlambda_F) - 2.5 min

% Equation (4): dlogP_F -> dlogP_CN
% dlogP_Vec = A_2*dlambda_F + B_2;
A_2 = Psi_total_tilde(C+1:end,:)*[zeros(C+CN,CF); eye(CF)]*diag(1./lambda_F);
B_2 = Psi_total_tilde(C+1:end,:)*[dX ; zeros(CF,1)];

A_2N = A_2(1:CN,:);
A_2F = A_2(CN+1:end,:);
B_2N = B_2(1:CN,1);
B_2F = B_2(CN+1:end,1);

% Equation (3): dlogP -> dOmega_tilde
% dOmega_CN_tilde_Vec = A_star*dlogP_Vec + B_star;

Omega_total_tilde_CN = Omega_total_tilde(1:C+CN,C+1:end);
Omega_total_tilde_CN_cell = num2cell(Omega_total_tilde_CN,1);
dlogttau = dlogt+dlogtau;

    if factor_index == 1 % country-specific factors (4 factors per country) : CN > CF

% A_star Ver.1 : 4 factors, slow
% Mat_C = permute(cat(2,bsxfun(@times, permute(repmat(Omega_total_tilde(1:C,C+1:C+CN),[1 1 CN]),[3 2 1]), AES_C_Mat-1), zeros(CN,CF,C)), [3 2 1]);
% Mat_N = permute(bsxfun(@times, permute(repmat(Omega_total_tilde(C+1:C+CN,C+1:end),[1 1 CN]),[3 2 1]), AES_N_Mat-1), [3 2 1]);
% Mat_F = permute(bsxfun(@times, permute(repmat(Omega_total_tilde(C+1:C+CN,C+1:end),[1 1 CF]),[3 2 1]), AES_F_Mat-1), [3 2 1]);
% Mat1 = cat(1,reshape(permute(cat(1,Mat_C,Mat_N),[2 1 3]),CN+CF,CN*(C+CN))',reshape(permute(cat(1,zeros(C,CN+CF,CF),Mat_F),[2 1 3]),CN+CF,CF*(C+CN))');

% A_star Ver.2 : 4 factors, fast
Mat = reshape(permute(repmat(Omega_total_tilde_CN,[1 1 CN+CF]).*cat(3,cat(1,permute(cat(2,AES_C_Mat-1,zeros(CN,CF,C)),[3 2 1]),permute(AES_N_Mat-1,[3 2 1])),cat(1,zeros(C,CN+CF,CF),permute(AES_F_Mat-1,[3 2 1]))), [2 1 3]), CN+CF, [])';

A_star = blkdiag(Omega_total_tilde_CN_cell{:})+repmat(reshape(Omega_total_tilde_CN,(C+CN)*(CN+CF),1),1,CN+CF).*Mat;
B_star = vec(Omega_total_tilde_CN.*dlogttau) + vec(Omega_total_tilde_CN).*sum(bsxfun(@times, Mat, repmat(dlogttau,CN+CF,1)),2);

    elseif factor_index == 2 % country-sector-specific factors (30 factors per country) : CN = CF

% A_star Ver.3 : Country-sector-specific factors, Seperated A_star and B_star, fast
Omega_total_tilde_CN_CN = Omega_total_tilde(1:C+CN,C+1:C+CN);
Omega_total_tilde_CN_CF_diag = diag(Omega_total_tilde(C+1:C+CN,C+CN+1:L));
Omega_total_tilde_CN_CN_cell = num2cell(Omega_total_tilde_CN_CN,1);

Mat_CN_CN = reshape(permute(repmat(Omega_total_tilde_CN_CN,[1 1 CN]).*cat(1,permute(AES_C_Mat-1,[3 2 1]),permute(AES_N_Mat(:,1:CN,:)-1,[3 2 1])), [2 1 3]), CN, [])';
Mat_CN_CF_temp = reshape(permute(AES_N_Mat(:,CN+1:end,:)-1,[3 2 1]),1,CN*CN,[]);
Mat_CN_CF = repmat(Omega_total_tilde_CN_CF_diag',CN,1).*reshape(Mat_CN_CF_temp(1,1:CN+1:end,:),CN,CN)';

% A_11 = A_star_CN_CN
A_star_CN_CN = blkdiag(Omega_total_tilde_CN_CN_cell{:})+repmat(reshape(Omega_total_tilde_CN_CN,(C+CN)*CN,1),1,CN).*Mat_CN_CN;

% A_12(Sparse) = A_star_CN_CF_sparse
A_star_CN_CF = Omega_total_tilde_CN_CN(C+1:end,:)'.*Mat_CN_CF;
A_star_CN_CF = A_star_CN_CF';
A_star_CN_CF_diag = full(sparse(1:numel(A_star_CN_CF), repmat(1:size(A_star_CN_CF,1),1,size(A_star_CN_CF,2)), A_star_CN_CF(:)));
A_star_CN_CF_sparse = sparse(reshape([zeros(C,CN*CF);reshape(A_star_CN_CF_diag,CN,[])],[],CF)); 
% A_star_CN_CF_stack = reshape([zeros(C,CN*CF);reshape(A_star_CN_CF_diag,CN,[])],[],CF);
% A_star_CN = [A_star_CN_CN, A_star_CN_CF_stack]; % A_star_up = [A_11 , A_12];

% B_1 = B_star_CN
B_star_CN_CN = vec(Omega_total_tilde_CN_CN.*dlogttau(:,1:CN)) + vec(Omega_total_tilde_CN_CN).*sum(bsxfun(@times, Mat_CN_CN, repmat(dlogttau(:,1:CN),CN,1)),2);
B_star_CN_CF = Omega_total_tilde_CN_CN(C+1:end,:)'.*Mat_CN_CF.*repmat(diag(dlogttau(C+1:C+CN,CN+1:end))',CN,1);
B_star_CN_CF = B_star_CN_CF';
B_star_CN_CF_stack = reshape([zeros(C,CN);B_star_CN_CF],[],1);
B_star_CN =  B_star_CN_CN + B_star_CN_CF_stack; % B_star_up = B_11 + B_12;

Mat_CF_CN = zeros(CN,CN);
Mat_CF_CF = zeros(CF,1);
for i=1:CF
    Mat_CF_CF(i,1) = Omega_total_tilde_CN_CF_diag(i,1)*(AES_F_Mat(i,CN+i,i)-1);
    Mat_CF_CN(i,:) = Omega_total_tilde_CN_CN(C+i,:).*(AES_F_Mat(i,1:CN,i)-1);
end

A_star_CF_CN = repmat(Omega_total_tilde_CN_CF_diag,1,CN).*Mat_CF_CN;
A_star_CF_CF = diag(Omega_total_tilde_CN_CF_diag+Omega_total_tilde_CN_CF_diag.*Mat_CF_CF);
A_star_CF_CN_stack = sparse([zeros(C,CN) ; reshape([reshape(A_star_CF_CN,1,[]);zeros(C+CN,(CN)*CN)],[],CN)]);
A_star_CF_CF_stack = sparse([zeros(C,CF) ; reshape([reshape(A_star_CF_CF,1,[]);zeros(C+CN,(CF)*CN)],[],CF)]);

% A_21(Sparse) = A_star_CF_CN_sparse
A_star_CF_CN_sparse = sparse(A_star_CF_CN_stack(1:(C+CN)*CF,:));

% A_22(Sparse) = A_star_CF_CF_sparse
A_star_CF_CF_sparse = sparse(A_star_CF_CF_stack(1:(C+CN)*CF,:));

B_star_CF_CN = Omega_total_tilde_CN_CF_diag.*sum(Mat_CF_CN.*dlogttau(C+1:C+CN,1:CN),2);
B_star_CF_CF = (Omega_total_tilde_CN_CF_diag+Omega_total_tilde_CN_CF_diag.*Mat_CF_CF).*diag(dlogttau(C+1:C+CN,CN+1:end));
B_star_CF = B_star_CF_CN + B_star_CF_CF;
B_star_CF_stack = [zeros(C,1) ; reshape([reshape(B_star_CF,1,[]);zeros(C+CN,CN)],[],1)];

% B_2(Sparse) = B_star_CF_sparse 
B_star_CF_sparse = sparse(B_star_CF_stack(1:(C+CN)*CF,:));

% A_star = [A_star_CN ; A_star_CF_stack(1:(C+CN)*CF,:)];
% B_star = [B_star_CN ; B_star_CF_stack(1:(C+CN)*CF,:)];

    end

% Equation (2): dOmega_tilde -> dOmega
% dOmega_CN_Vec = A_1*dOmega_CN_tilde_Vec + B_1
% dOmega_CN_Vec = A_Omega*dlambda_F + B_Omega;

    if factor_index == 1

rA_1 = vec(Omega_total(1:C+CN,C+1:end)./Omega_total_tilde(1:C+CN,C+1:end));
rA_1(isnan(rA_1))=0;
B_1 = -vec(Omega_total(1:C+CN,C+1:end)).*vec(dlogt);

A_Omega = (repmat(rA_1,1,size(A_star,2)).*A_star)*A_2;
B_Omega = B_1 + repmat(rA_1,1,size(B_star,2)).*B_star + (repmat(rA_1,1,size(A_star,2)).*A_star)*B_2;

    elseif factor_index == 2

rA_1N = vec(Omega_total(1:C+CN,C+1:C+CN)./Omega_total_tilde(1:C+CN,C+1:C+CN));
rA_1F = sparse(vec(Omega_total(1:C+CN,C+CN+1:end)./Omega_total_tilde(1:C+CN,C+CN+1:end)));
rA_1N(isnan(rA_1N))=0;
rA_1F(isnan(rA_1F))=0;
B_1N = -vec(Omega_total(1:C+CN,C+1:C+CN)).*vec(dlogt(:,1:CN));
B_1F = sparse(-vec(Omega_total(1:C+CN,C+CN+1:end)).*vec(dlogt(:,CN+1:end)));

clearvars AES_C_Mat AES_N_Mat AES_F_Mat Mat_CN_CN Mat_CN_CF Mat_CN_CF_temp
clearvars A_star_CN_CF_diag A_star_CF_CN_stack A_star_CF_CF_stack

A_Omega = [(repmat(rA_1N,1,size(A_star_CN_CN,2)).*A_star_CN_CN)*A_2N+(repmat(rA_1N,1,size(A_star_CN_CF_sparse,2)).*A_star_CN_CF_sparse)*A_2F ;
                 (repmat(rA_1F,1,size(A_star_CF_CN_sparse,2)).*A_star_CF_CN_sparse)*A_2N+(repmat(rA_1F,1,size(A_star_CF_CF_sparse,2)).*A_star_CF_CF_sparse)*A_2F ];
B_Omega = [B_1N + repmat(rA_1N,1,size(B_star_CN,2)).*B_star_CN + (repmat(rA_1N,1,size(A_star_CN_CN,2)).*A_star_CN_CN)*B_2N+(repmat(rA_1N,1,size(A_star_CN_CF_sparse,2)).*A_star_CN_CF_sparse)*B_2F ;
                 B_1F + repmat(rA_1F,1,size(B_star_CF_sparse,2)).*B_star_CF_sparse + (repmat(rA_1F,1,size(A_star_CF_CN_sparse,2)).*A_star_CF_CN_sparse)*B_2N+(repmat(rA_1F,1,size(A_star_CF_CF_sparse,2)).*A_star_CF_CF_sparse)*B_2F ];

    end

% A_Omega_N =(repmat(rA_1N,1,size(A_star_CN_CN,2)).*A_star_CN_CN)*A_2N+(repmat(rA_1N,1,size(A_star_CN_CF_sparse,2)).*A_star_CN_CF_sparse)*A_2F;
% A_Omega_F = sparse((repmat(rA_1F,1,size(A_star_CF_CN_sparse,2)).*A_star_CF_CN_sparse)*A_2N+(repmat(rA_1F,1,size(A_star_CF_CF_sparse,2)).*A_star_CF_CF_sparse)*A_2F);
% Equation (1): dOmega -> dPsi
% dPsi_total_Vec = A_3*dOmega_CN_Vec

% A_3 = kron(Psi_total',Psi_total)*kron([eye(CN+CF); zeros(C,CN+CF)],[eye(C+CN); zeros(CF,C+CN)]);

%% dchi_std = g(dlambda_F, dlambda_F_star) - 2 min
% Equation (5): dlambda_F, dlambda_F_star -> dchi
% ksi = G*dlambda_F + H;
% dchi = A_5*dlambda_F + dlambda_F_star + B_5

dt = init_t.*dlogt; % Matrix of Size (C+CN,CN+CF)

H = (lambda_CN'*squeeze(sum(Phi_T.*repmat(reshape(B_Omega,C+CN,CN+CF),[1 1 C]).*repmat(init_t-1,[1 1 C]),2)))';
G = zeros(C,CF);

% G Ver.1 : Obtain each column of G (Very slow when factors become sector-specific)
% for i=1:CF
%     temp = zeros(CF,1);
%     temp(i) = 1;
%     dOmega_CN_temp = reshape(A_Omega*temp + B_Omega,C+CN,CN+CF);
%     GH_temp = (lambda_CN'*squeeze(sum(Phi_T.*repmat(dOmega_CN_temp,[1 1 C]).*repmat(init_t-1,[1 1 C]),2)))';
%     G(:,i) = GH_temp - H;
% end

% G Ver.2 : Obtain each row of G
for j=1:C
    G(j,:) = sum(bsxfun(@times,A_Omega,vec(Phi_T(:,:,j).*(init_t-1).*repmat(lambda_CN,1,CN+CF))),1);
end

% G = squeeze(sum(bsxfun(@times,repmat(A_Omega,[1 1 C]),reshape(Phi_T,(C+CN)*(CN+CF),1,C).*repmat(vec(init_t-1),[1 1 C]).*repmat(repmat(lambda_CN,CN+CF,1),[1 1 C])),1))';

temp1 = (lambda_CN'*squeeze(sum(Phi_T.*repmat(Omega_total(1:C+CN,C+1:end),[1 1 C]).*repmat(dt,[1 1 C]),2)))';

A_5 = Phi_F + G;
B_5 = temp1 + H;

%% Output - 0.5 min
% dlambda_F = A_up*[dlambda_F ; dlambda_F_star]+B_up
% dlambda_CN = A_CN*[dlambda_F ; dlambda_F_star]+B_CN 
% dlambda_F_star = A_down*[dlambda_F ; dlambda_F_star]+B_down

A_up = [Psi_total(:,C+CN+1:end)'*[eye(C); zeros(L-C,C)]*A_5 + kron(Psi_total(:,C+CN+1:end)'*[zeros(C,CN+CF); eye(CN+CF)],lambda_CN')*A_Omega, Psi_total(:,C+CN+1:end)'*[eye(C); zeros(L-C,C)]];
B_up = Psi_total(:,C+CN+1:end)'*[eye(C); zeros(L-C,C)]*B_5 + kron(Psi_total(:,C+CN+1:end)'*[zeros(C,CN+CF); eye(CN+CF)],lambda_CN')*B_Omega;

A_CN = [Psi_total(:,1:C+CN)'*[eye(C); zeros(L-C,C)]*A_5 + kron(Psi_total(:,1:C+CN)'*[zeros(C,CN+CF) ; eye(CN+CF)],lambda_CN')*A_Omega, Psi_total(:,1:C+CN)'*[eye(C); zeros(L-C,C)]];
B_CN = Psi_total(:,1:C+CN)'*[eye(C); zeros(L-C,C)]*B_5 + kron(Psi_total(:,1:C+CN)'*[zeros(C,CN+CF) ; eye(CN+CF)],lambda_CN')*B_Omega;

A_down = squeeze(sum(Phi_T.*repmat(Omega_total(1:C+CN,C+1:end),[1 1 C]).*repmat(init_t-1,[1 1 C]),2))'*A_CN;
B_down = squeeze(sum(Phi_T.*repmat(Omega_total(1:C+CN,C+1:end),[1 1 C]).*repmat(init_t-1,[1 1 C]),2))'*B_CN;

A_result = [A_up ; A_down];
B_result = [B_up ; B_down];

A_orig = A_result;
B_orig = B_result;

% Substitute first equation to \sum dchi = 0
CC = [sum(G,1), zeros(1,C)];
D = sum(H,1);

A_result(1,:) = -1/(1+CC(1,1))*[0, 1+CC(1,2:CF), ones(1,C)];
B_result(1) = -(D+sum(temp1))/(1+CC(1,1));
end