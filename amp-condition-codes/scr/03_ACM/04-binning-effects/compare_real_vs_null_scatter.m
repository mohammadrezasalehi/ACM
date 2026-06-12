% Compares Real Amplitude-Modulated results against Random Null Model.

clc; clear; close all;

%% 1. Configuration
project_root = 'F:\PhD Code\My_PhD_Project';
session = 'REST1_RL'; 

% Paths to Group Level directories
real_dir = fullfile(project_root, 'results', 'Modulated_FC', 'Group_Level', session);
null_dir = fullfile(project_root, 'results', 'Modulated_FC_RandomNull', 'Group_Level', session);

%% 2. Load Data Pairs
fprintf('Loading Group-Level matrices...\n');

% Variance (Consensus)
V_real = load(fullfile(real_dir, 'Group_Variance_Consensus.mat')); v_r = V_real.Consensus_Percent(:);
V_null = load(fullfile(null_dir, 'Group_Variance_Consensus.mat')); v_n = V_null.Consensus_Percent(:);

% Slope (T-map)
T_real = load(fullfile(real_dir, 'Group_Slope_Tmap.mat')); t_r = T_real.T_map(:);
T_null = load(fullfile(null_dir, 'Group_Slope_Tmap.mat')); t_n = T_null.T_map(:);

% ZCR (Median)
Z_real = load(fullfile(real_dir, 'Group_ZCR_Median.mat')); z_r = Z_real.Median(:);
Z_null = load(fullfile(null_dir, 'Group_ZCR_Median.mat')); z_n = Z_null.Median(:);

% Poly (Prevalence)
P_real = load(fullfile(real_dir, 'Group_Poly_Prevalence.mat')); p_r = P_real.Prevalence_Percent(:);
P_null = load(fullfile(null_dir, 'Group_Poly_Prevalence.mat')); p_n = P_null.Prevalence_Percent(:);

%% 3. Filter NaNs (Important for correlation and plotting)
v_valid = ~isnan(v_r) & ~isnan(v_n);
t_valid = ~isnan(t_r) & ~isnan(t_n);
z_valid = ~isnan(z_r) & ~isnan(z_n);
p_valid = ~isnan(p_r) & ~isnan(p_n);

%% 4. Visualization (2x2 Binscatter)
figure('Name', 'Real vs Random Null Comparison', 'Position', [100, 100, 1000, 900], 'Color', 'w');
sgtitle('Sensitivity to Amplitude: Real Modulator vs. Random Shuffling', 'FontSize', 16, 'FontWeight', 'bold');

% --- Subplot 1: Variance Consensus ---
subplot(2,2,1);
binscatter(v_n(v_valid), v_r(v_valid), [100 100]); % 100x100 bins for heat
colormap(jet); hold on;
plot(xlim, xlim, 'r--', 'LineWidth', 2); % Identity line (y=x)
r_val = corr(v_n(v_valid), v_r(v_valid));
title(sprintf('Variance Consensus Map (r = %.2f)', r_val));
xlabel('Random Null Model'); ylabel('Real Amplitude Model');

% --- Subplot 2: Slope T-Map ---
subplot(2,2,2);
binscatter(t_n(t_valid), t_r(t_valid), [100 100]);
colormap(jet); hold on;
plot([-10 10], [-10 10], 'r--', 'LineWidth', 2); % Identity line
xlim([-15 15]); ylim([-15 15]);
r_val = corr(t_n(t_valid), t_r(t_valid));
title(sprintf('Slope T-Value Map (r = %.2f)', r_val));
xlabel('Random Null Model'); ylabel('Real Amplitude Model');

% --- Subplot 3: ZCR Median ---
subplot(2,2,3);
binscatter(z_n(z_valid), z_r(z_valid), [100 100]);
colormap(jet); hold on;
plot(xlim, xlim, 'r--', 'LineWidth', 2);
r_val = corr(z_n(z_valid), z_r(z_valid));
title(sprintf('ZCR/MCR Median (r = %.2f)', r_val));
xlabel('Random Null Model'); ylabel('Real Amplitude Model');

% --- Subplot 4: Non-Linear Prevalence ---
subplot(2,2,4);
binscatter(p_n(p_valid), p_r(p_valid), [100 100]);
colormap(jet); hold on;
plot(xlim, xlim, 'r--', 'LineWidth', 2);
r_val = corr(p_n(p_valid), p_r(p_valid));
title(sprintf('Non-Linear Prevalence (r = %.2f)', r_val));
xlabel('Random Null Model'); ylabel('Real Amplitude Model');