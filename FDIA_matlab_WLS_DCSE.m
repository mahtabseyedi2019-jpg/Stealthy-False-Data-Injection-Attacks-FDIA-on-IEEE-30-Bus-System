close all;
clear all;
clc;
format compact;
path2tools = 'PR_Toolbox\';
addpath(path2tools);

mpc = loadcase('case30');
n_loadbus = sum(mpc.bus(:,2) == 1);

loads = load_prep_data('D:\payan name\متلبی\nyiso_2013\*.csv', 61757);
%%
loads1 = loads(1:30700,:) 
total1 = get_state_vars_with_load( mpc, loads1 );

%%
loads2 = loads(10001:44600,:) 
total2 = get_state_vars_with_load(mpc, loads2) %include of H,R_inv,z,x_se,x_DCPF,diff(x_SE,x_DCPF)
total = total1;  
save('nyiso_load_statevars', 'total');

%%

load nyiso_load_statevars

%%
H = total.H(:,:,1); %for frist set of load
R_inv = total.R_Inv(:,:,1);
z = total.z(1,:)';
x_rad = total.x_se(1,:)'; 



c = zeros(size(H,2),1);
c(1) = 2;
c(5) = 2;
a = H*c;
z_a = z+a;
%%Weighted least square (WLS) 

x = (H'*R_inv*H)^-1 * H'*R_inv * z;
x_bad = (H'*R_inv*H)^-1 * H'*R_inv * z_a;
disp('residual of noraml data(WLS)')
norm(z-H*x)
disp('residual of attack data(WLS)')
norm(z_a-H*x_bad)


%% result check


define_constants;
total_load = sum(mpc.bus(:,PD));
total_gen = sum(mpc.gen(:,PG));
disp(['Total Load: ', num2str(total_load)]);
disp(['Total Gen: ', num2str(total_gen)]);
disp(['Mismatch (Load - Gen): ', num2str(total_load - total_gen)]);% power blance
disp('Min and Max of estimated angles (radians):');
disp([min(x), max(x)]);
density_H = nnz(H) / numel(H); 
disp(['Density of H: ', num2str(density_H)]);

tau = 3.0;  

R = pinv(R_inv); 

G = H' * R_inv * H;
L = chol(G, 'lower');  
n_meas = size(H,1); 
omega_diag = zeros(n_meas,1);

for m = 1:n_meas
    h_m = H(m,:)';  
    T_m = L' \ (L \ h_m);  %G^-1
    omega_diag(m) = R(m,m) - h_m' * T_m;
    
    if omega_diag(m) <= 0
        omega_diag(m) = 1e-6;
    end
end

r_normal = z - H * x;       
r_attacked = z_a - H * x_bad ; 

r_norm_normal = r_normal ./ sqrt(omega_diag);
r_norm_attacked = r_attacked ./ sqrt(omega_diag);

[max_rn, idx_rn] = max(abs(r_norm_normal));
[max_rn_a, idx_rn_a] = max(abs(r_norm_attacked));

disp('-----------------------------')
disp('LNR Test Result (Normal Data)')
disp(['Max normalized residual: ', num2str(max_rn)])
disp(['Index of max residual: ', num2str(idx_rn)])
if max_rn > tau
    disp('Detected bad data in normal measurements (False Alarm)');
else
    disp('No bad data detected in normal data');
end

disp('-----------------------------')
disp('LNR Test Result (Attacked Data)')
disp(['Max normalized residual: ', num2str(max_rn_a)])
disp(['Index of max residual: ', num2str(idx_rn_a)])
if max_rn_a > tau
    disp('Detected bad data in attacked measurements (Correct Detection)');
else
    disp('No bad data detected in attacked data (Missed Detection)');
end
disp('-----------------------------')



cond_G = cond(G);
disp(['Condition number of Gain matrix: ', num2str(cond_G)]);
if cond_G > 1e12
    warning('⚠Gain matrix is ill-conditioned. WLS solution may be inaccurate.');
end
g = H' * R_inv * z;
x = G \ g;

r = z - H * x;

J_normal = r' * R_inv * r;
J_normal = full(J_normal(1)); 

k = size(H,1) - size(H,2);  

confidence_level = 0.95;
threshold_chi2 = chi2inv(confidence_level, k);

disp('==============================');
disp('Chi-Squared Test Result for Normal Data:');
disp(['Threshold (t_J) at 95% confidence: ', num2str(threshold_chi2)]);
disp(['J value for normal data: ', num2str(J_normal)]);
if (J_normal > threshold_chi2) && (J_normal > 1e-6)
    disp('Bad data detected in normal measurements (False Alarm)!');
else
    disp('No bad data detected in normal measurements.');
end
disp('==============================');


z_a = z_a(:);  
g_a = H' * R_inv * z_a;
x_a = G \ g_a;

r_a = z_a - H * x_a;

J_attacked = r_a' * R_inv * r_a;
J_attacked = full(J_attacked(1));  

disp('==============================');
disp('Chi-Squared Test Result for Attacked Data:');
disp(['Threshold (t_J) at 95% confidence: ', num2str(threshold_chi2)]);
disp(['J value for attacked data: ', num2str(J_attacked)]);
if (J_attacked > threshold_chi2) && (J_attacked > 1e-6)
    disp('Bad data detected in attacked measurements (Correct Detection)!');
else
    disp('No bad data detected in attacked measurements (Missed Detection).');
end
disp('==============================');

load nyiso_load_statevars

threshold_LNR = 3.0;
confidence_level = 0.95;
attack_prob = 0.05; 
timesteps = size(total.z, 1);

X = zeros(size(total.z));
Y = zeros(timesteps, 1);
is_attacked = zeros(timesteps, 1);
Detected_LNR = zeros(timesteps, 1);
Detected_ChiSq = zeros(timesteps, 1);
LNR_results = zeros(timesteps, 1);
ChiSq_results = zeros(timesteps, 1);

for t = 1:timesteps
    H = total.H(:,:,t);
    R_inv = total.R_Inv(:,:,t);
    z = total.z(t,:)';

    if rand < attack_prob
        c = 2 * ones(size(H,2),1);  
        a = H * c;
        z_a = z + a+0.1;
        is_attacked(t) = 1;
        Y(t) = 1;
    else
        z_a = z;
        Y(t) = 0;
    end
    X(t,:) = z_a';

    G = H' * R_inv * H;
    g = H' * R_inv * z_a;
    x_hat = G \ g;
    r = z_a - H * x_hat;

    R = pinv(R_inv);
    omega_diag = zeros(size(H,1),1);
    try
        L = chol(G, 'lower');
        for m = 1:size(H,1)
            h_m = H(m,:)';
            T_m = L' \ (L \ h_m);
            omega_diag(m) = R(m,m) - h_m' * T_m;
            if omega_diag(m) <= 0
                omega_diag(m) = 1e-6;
            end
        end
    catch
        omega_diag(:) = max(diag(R), 1e-6);
    end
    r_norm = r ./ sqrt(omega_diag);
    LNR_results(t) = max(abs(r_norm));
    Detected_LNR(t) = LNR_results(t) > threshold_LNR;

    J = r' * R_inv * r;
    ChiSq_results(t) = full(J);
    k = size(H,1) - size(H,2);
    threshold_chi2 = chi2inv(confidence_level, k);
    Detected_ChiSq(t) = (J > threshold_chi2) && (J > 1e-6);
end

disp('Randomized attacks applied at timesteps:');
disp(find(is_attacked));
disp('Detected by LNR:');
disp(find(Detected_LNR));
disp('Detected by Chi-Squared:');
disp(find(Detected_ChiSq));



%%
function [ total ] = get_state_vars_with_load( mpc, loads )

define_constants;

n_loadbus = sum(mpc.bus(:,2) == 1); 
idx_loadbus = find(mpc.bus(:,2) == 1); 
n_genbus = sum(mpc.bus(:,2) == 2);
idx_genbus = find(mpc.bus(:,2) == 2);


PD_original = mpc.bus(idx_loadbus, PD); 
max_load = max(loads(:));  
t0 = tic;
for set=1:size(loads,1) 
    for msmt=1:size(loads,2)
        normalized_factor = loads(set, msmt) / max_load; 
        new_PD = PD_original(msmt) * normalized_factor;   
        mpc.bus(idx_loadbus(msmt), PD) = new_PD;
       mpc.bus(idx_loadbus(msmt),QD) = 0;
       mpc.bus(idx_loadbus(msmt),GS) = 0;
       mpc.bus(idx_loadbus(msmt),BS) = 0;
    end

    current_load = sum(mpc.bus(:,3));
    current_gen = sum(mpc.gen(:,2));

    inc_per_gen = (current_load - current_gen) / n_genbus;

    for msmt=2:n_genbus+1
        mpc.gen(msmt,2) = mpc.gen(msmt,2) + inc_per_gen;
    end

    mpc.baseMVA = 100;
   
    
    results = rundcpf(mpc, mpoption('out.all',0));
    [ x, H, R_Inv, z ] = dc_state_est(mpc, results);
    
    disp(set);
    
    total.H(:,:,set) = full(H);
    total.R_Inv(:,:,set) = R_Inv;
    total.z(set,:) = z;
    total.x_se(set,:) = x;
    total.x_pf(set,:) = (pi/180) .* results.bus(2:end,9);
    total.x_diff(set,:) = total.x_se(set,:) - total.x_pf(set,:);
end

toc(t0)

end



%%
function loads = load_prep_data(folder_path, ptid_filter)


if nargin < 1 
    error('Folder path to CSV files is required.');
end

if nargin < 2 
    ptid_filter = 61757; 
end
files = dir(folder_path); 
if isempty(files) 
    error('No CSV files found in the specified folder.');
end

loads_all = [];

for k = 1:length(files) 
    file_name = fullfile(files(k).folder, files(k).name); 
    try
        data = readmatrix(file_name); 
    catch 
        warning('Could not read file: %s. Skipping...', file_name);
        continue;
    end

    [rows, cols] = size(data);
    if cols < 5
        warning('File %s has fewer than 5 columns. Skipping...', file_name);
        continue;
    end

    filtered_rows = data(data(:,4) == ptid_filter, :);

    loads_all = [loads_all; filtered_rows(:,5)];
end

loads = loads_all;
end



%%
function [x, H, R_inv, z] = dc_state_est(mpc, results) 
    define_constants; 

    Va = results.bus(:, VA); 
    Va_rad = deg2rad(Va);    
    
    ref_bus = find(mpc.bus(:, BUS_TYPE) == REF, 1); 
    if isempty(ref_bus) 
        warning('No reference bus found. Assuming bus 1 as reference.');
        ref_bus = 1;
    end
    
   
    n_bus = size(mpc.bus, 1); %number of rows=number of buss
    x =Va_rad(setdiff(1:n_bus, ref_bus)); %setdiff=hame anasor nbus ke dar refbus nist->delet voltage of ref bus


    branch = mpc.branch;
    from = branch(:, F_BUS);
    to = branch(:, T_BUS);
    n_branch = size(branch, 1);
    n_bus = size(mpc.bus, 1);

    %%  susptance cal (b = -1/X)
    b = -1 ./ branch(:, BR_X); % in DCSE: P_ij = -b_ij*(θ_i - θ_j)
    b(isinf(b)) = 0; % line with ampdance=0->isinf(b):if b=inf replace with 0

    %%  make H (H)
    H = zeros(n_branch, n_bus - 1);
    for k = 1:n_branch
        i = from(k);
        j = to(k);
        
        
        if i ~= ref_bus
            col_i = i - (i > ref_bus);
            H(k, col_i) = b(k);
        end
        if j ~= ref_bus
            col_j = j - (j > ref_bus);
            H(k, col_j) = -b(k);
        end
    end

    %% Power flow
    try
        z = results.branch(:, PF) / mpc.baseMVA; %p.u result of DCPF (just Pij)
    catch
       
        z = zeros(n_branch, 1);
        for k = 1:n_branch
            i = from(k);
            j = to(k);
            z(k) = b(k) * (Va_rad(i) - Va_rad(j)); %DCPF
        end
        z = z / mpc.baseMVA;
    end

    %% (R_inv)
    measurement_error = 0.01 * max(abs(z), 0.1); % min eror 1% 
    R_inv = diag(1 ./ (measurement_error.^2));  %diag make a qotri matrix(inv(cov(eror))
end

disp('Extracting variables from struct...');
H_3D = total.H;
z_2D = total.z;
x_se_2D = total.x_se;
R_Inv_3D = total.R_Inv;


disp('Saving 2D matrices (z, x_se) to CSV...');
writematrix(z_2D, 'z_data.csv');
writematrix(x_se_2D, 'x_se_data.csv');
disp('✅ z_data.csv and x_se_data.csv saved.');


disp('Reshaping and saving 3D matrices (H, R_Inv) to CSV...');

[rows_H, cols_H, num_timesteps_H] = size(H_3D);
H_reshaped_2D = reshape(permute(H_3D, [1, 3, 2]), [], cols_H);
writematrix(H_reshaped_2D, 'H_data_reshaped.csv');
disp('✅ H_data_reshaped.csv saved.');

[rows_R, cols_R, num_timesteps_R] = size(R_Inv_3D);
R_Inv_reshaped_2D = reshape(permute(R_Inv_3D, [1, 3, 2]), [], cols_R);
writematrix(R_Inv_reshaped_2D, 'R_Inv_data_reshaped.csv');
disp('✅ R_Inv_data_reshaped.csv saved.');

disp('--- All CSV files saved successfully! ---');