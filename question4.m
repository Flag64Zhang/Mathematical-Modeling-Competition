% 投票选拔问题 - 决策变量n、k/m、α、m及m、n联合分布的灵敏度分析
% 文件名: sensitivity_analysis.m

clear; clc; close all;
%% 全局参数
N = 19;   % 专家总人数
M = 10000; % 蒙特卡洛模拟次数
%% ----------------- 灵敏度分析（n变化） -----------------
fprintf('\n=== 灵敏度分析 (n变化) ===\n');
m = 5; k = 3; % 固定为实例 A
n_range = 10:19; % 出席人数
P_sim_list = zeros(size(n_range));
P_theory_list = zeros(size(n_range));

for i = 1:length(n_range)
    n = n_range(i);
    t = ceil(2*n/3);
    p = k / m;
    P_theory = 1 - binocdf(t-1, n, p);
    P_theory_list(i) = P_theory;

    % 蒙特卡洛
    success_count = 0;
    for iter = 1:M
        votes = zeros(1,m);
        for expert = 1:n
            picks = randperm(m,k);
            votes(picks) = votes(picks) + 1;
        end
        success_count = success_count + (votes(1) >= t);
    end
    P_sim = success_count / M;
    P_sim_list(i) = P_sim;

    fprintf('n=%d: 理论=%.4f, 模拟=%.4f\n', n, P_theory, P_sim);
end

% 计算皮尔逊相关系数
R_n_theory = corr(n_range', P_theory_list', 'Type', 'Pearson');
R_n_sim = corr(n_range', P_sim_list', 'Type', 'Pearson');
fprintf('皮尔逊相关系数（理论）：R = %.4f\n', R_n_theory);
fprintf('皮尔逊相关系数（模拟）：R = %.4f\n', R_n_sim);

figure;
plot(n_range, P_theory_list, '-o', 'LineWidth',1.5); hold on;
plot(n_range, P_sim_list, '-s', 'LineWidth',1.5);
xlabel('实到专家人数 n');
ylabel('单候选人当选概率');
legend('理论(Binom)', '模拟(MC)');
title(sprintf('灵敏度分析：出席人数 n 对当选概率的影响 (R理论=%.2f, R模拟=%.2f)', R_n_theory, R_n_sim));
grid on;

%% ----------------- 灵敏度分析（k/m变化） -----------------
fprintf('\n=== 灵敏度分析 (k/m变化) ===\n');
n = 18; alpha = 2/3; m = 10;
km_range = 0.4:0.05:0.8;
P_sim_list_k = zeros(size(km_range));
P_theory_list_k = zeros(size(km_range));

for i = 1:length(km_range)
    k = round(km_range(i) * m);
    t = ceil(alpha * n);
    p = k / m;
    P_theory = 1 - binocdf(t-1, n, p);
    P_theory_list_k(i) = P_theory;

    % 蒙特卡洛
    success_count = 0;
    for iter = 1:M
        votes = zeros(1,m);
        for expert = 1:n
            picks = randperm(m,k);
            votes(picks) = votes(picks) + 1;
        end
        success_count = success_count + (votes(1) >= t);
    end
    P_sim = success_count / M;
    P_sim_list_k(i) = P_sim;

    fprintf('k=%d (k/m=%.2f): 理论=%.4f, 模拟=%.4f\n', k, k/m, P_theory, P_sim);
end

% 计算皮尔逊相关系数
R_km_theory = corr(km_range', P_theory_list_k', 'Type', 'Pearson');
R_km_sim = corr(km_range', P_sim_list_k', 'Type', 'Pearson');
fprintf('皮尔逊相关系数（理论）：R = %.4f\n', R_km_theory);
fprintf('皮尔逊相关系数（模拟）：R = %.4f\n', R_km_sim);

figure;
plot(km_range, P_theory_list_k, '-o', 'LineWidth',1.5); hold on;
plot(km_range, P_sim_list_k, '-s', 'LineWidth',1.5);
xlabel('投票比例 k/m');
ylabel('单候选人当选概率');
legend('理论(Binom)', '模拟(MC)');
title(sprintf('灵敏度分析：投票比例 k/m 对当选概率的影响 (n=18, m=10, α=2/3, R理论=%.2f, R模拟=%.2f)', R_km_theory, R_km_sim));
grid on;

%% ----------------- 灵敏度分析（阈值比例α变化） -----------------
fprintf('\n=== 灵敏度分析 (阈值比例α变化) ===\n');
n = 18; m = 5; k = 3; % 按照建模方案进行构造
alpha_range = 0.5:0.05:0.9; % 阈值比例
P_sim_list_alpha = zeros(size(alpha_range));
P_theory_list_alpha = zeros(size(alpha_range));

for i = 1:length(alpha_range)
    alpha = alpha_range(i);
    t = ceil(alpha * n);
    p = k / m;
    P_theory = 1 - binocdf(t-1, n, p);
    P_theory_list_alpha(i) = P_theory;

    % 蒙特卡洛
    success_count = 0;
    for iter = 1:M
        votes = zeros(1,m);
        for expert = 1:n
            picks = randperm(m,k);
            votes(picks) = votes(picks) + 1;
        end
        success_count = success_count + (votes(1) >= t);
    end
    P_sim = success_count / M;
    P_sim_list_alpha(i) = P_sim;

    fprintf('α=%.2f, t=%d: 理论=%.4f, 模拟=%.4f\n', alpha, t, P_theory, P_sim);
end

% 计算皮尔逊相关系数
R_alpha_theory = corr(alpha_range', P_theory_list_alpha', 'Type', 'Pearson');
R_alpha_sim = corr(alpha_range', P_sim_list_alpha', 'Type', 'Pearson');
fprintf('皮尔逊相关系数（理论）：R = %.4f\n', R_alpha_theory);
fprintf('皮尔逊相关系数（模拟）：R = %.4f\n', R_alpha_sim);

figure;
plot(alpha_range, P_theory_list_alpha, '-o', 'LineWidth',1.5); hold on;
plot(alpha_range, P_sim_list_alpha, '-s', 'LineWidth',1.5);
xlabel('阈值比例 α');
ylabel('单候选人当选概率');
legend('理论(Binom)', '模拟(MC)');
title(sprintf('灵敏度分析：阈值比例 α 对当选概率的影响 (R理论=%.2f, R模拟=%.2f)', R_alpha_theory, R_alpha_sim));
grid on;


%% ----------------- 灵敏度分析（候选人数m变化） -----------------
fprintf('\n=== 灵敏度分析 (候选人数m变化) ===\n');
m_samples = 3:12;
n = 15; alpha = 2/3;
k_samples = ceil(m_samples/2);
P_theory_samples = zeros(size(m_samples));
P_sim_samples = zeros(size(m_samples));

for i = 1:length(m_samples)
    m = m_samples(i);
    k = k_samples(i);
    t = ceil(alpha * n);
    p = k / m;
    P_theory = 1 - binocdf(t-1, n, p);
    P_theory_samples(i) = P_theory;

    % 蒙特卡洛
    success_count = 0;
    for iter = 1:M
        votes = zeros(1,m);
        for expert = 1:n
            picks = randperm(m,k);
            votes(picks) = votes(picks) + 1;
        end
        success_count = success_count + (votes(1) >= t);
    end
    P_sim = success_count / M;
    P_sim_samples(i) = P_sim;
     fprintf('α=%.2f, t=%d: 理论=%.4f, 模拟=%.4f\n', alpha, t, P_theory, P_sim);
end

% 计算皮尔逊相关系数
R_theory = corr(m_samples', P_theory_samples', 'Type', 'Pearson');
R_sim = corr(m_samples', P_sim_samples', 'Type', 'Pearson');
fprintf('皮尔逊相关系数（理论）：R = %.4f\n', R_theory);
fprintf('皮尔逊相关系数（模拟）：R = %.4f\n', R_sim);

figure;
plot(m_samples, P_theory_samples, '-o', 'LineWidth',1.5); hold on;
plot(m_samples, P_sim_samples, '-s', 'LineWidth',1.5);
xlabel('候选人数 m');
ylabel('单候选人当选概率');
legend('理论(Binom)', '模拟(MC)');
title(sprintf('全局敏感性分析：m对P的影响 (R理论=%.2f, R模拟=%.2f)', R_theory, R_sim));
grid on;


%% ----------------- 全局灵敏度分析（m和n联合分布） -----------------
fprintf('\n=== 全局灵敏度分析 (m和n联合分布对P的影响) ===\n');
m_samples_joint = 3:12;
n_samples_joint = 10:19;
alpha = 2/3;
P_theory_joint = zeros(length(m_samples_joint), length(n_samples_joint));
P_sim_joint = zeros(length(m_samples_joint), length(n_samples_joint));

for i = 1:length(m_samples_joint)
    for j = 1:length(n_samples_joint)
        m = m_samples_joint(i);
        n = n_samples_joint(j);
        k = ceil(m/2);
        t = ceil(alpha * n);
        p = k / m;
        P_theory = 1 - binocdf(t-1, n, p);
        P_theory_joint(i,j) = P_theory;

        % 蒙特卡洛
        success_count = 0;
        for iter = 1:M
            votes = zeros(1,m);
            for expert = 1:n
                picks = randperm(m,k);
                votes(picks) = votes(picks) + 1;
            end
            success_count = success_count + (votes(1) >= t);
        end
        P_sim = success_count / M;
        P_sim_joint(i,j) = P_sim;
    end
end

% 修正展开方式，保证与矩阵一致
[m_grid, n_grid] = meshgrid(m_samples_joint, n_samples_joint);
m_flat = m_grid(:);
n_flat = n_grid(:);
P_theory_flat = P_theory_joint'; P_theory_flat = P_theory_flat(:);
P_sim_flat = P_sim_joint'; P_sim_flat = P_sim_flat(:);

% 计算皮尔逊相关系数
R_m_theory = corr(m_flat, P_theory_flat, 'Type', 'Pearson');
R_n_theory = corr(n_flat, P_theory_flat, 'Type', 'Pearson');
R_m_sim = corr(m_flat, P_sim_flat, 'Type', 'Pearson');
R_n_sim = corr(n_flat, P_sim_flat, 'Type', 'Pearson');

fprintf('皮尔逊相关系数（理论）：R_m = %.4f, R_n = %.4f\n', R_m_theory, R_n_theory);
fprintf('皮尔逊相关系数（模拟）：R_m = %.4f, R_n = %.4f\n', R_m_sim, R_n_sim);

% 绘制3D散点图
figure;
scatter3(m_flat, n_flat, P_theory_flat, 50, P_theory_flat, 'filled');
xlabel('候选人数 m');
ylabel('实到专家人数 n');
zlabel('单候选人当选概率');
title(sprintf('全局敏感性分析：m和n联合对P的影响 (理论, R_m=%.2f, R_n=%.2f)', R_m_theory, R_n_theory));
colorbar;
grid on;

figure;
scatter3(m_flat, n_flat, P_sim_flat, 50, P_sim_flat, 'filled');
xlabel('候选人数 m');
ylabel('实到专家人数 n');
zlabel('单候选人当选概率');
title(sprintf('全局敏感性分析：m和n联合对P的影响 (模拟, R_m=%.2f, R_n=%.2f)', R_m_sim, R_n_sim));
colorbar;
grid on;