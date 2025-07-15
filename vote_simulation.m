% 投票选拔问题 - MATLAB完整示例（实例A+B+灵敏度分析）
% 文件名: vote_simulation.m

clear; clc; close all;

%% 全局参数
N = 19;   % 专家总人数
M = 10000; % 蒙特卡洛模拟次数

%% ----------------- 实例 A -----------------
% 参数
n_A = 15;          % 实到专家人数
m_A = 5;           % 候选人数
k_A = floor(m_A/2)+1; % 每位专家投票数
t_A = ceil(2*n_A/3);  % 当选票数阈值
s_A = floor(m_A/2);   % 推优名额

fprintf('=== 实例 A ===\n');
run_simulation(N, n_A, m_A, k_A, t_A, M, 'A');

% 绘制实例A中单候选人当选概率随n增长的折线图
fprintf('\n=== 实例A：单候选人当选概率随n增长 ===\n');
n_plot = 10:19;
P_n_theory = zeros(size(n_plot));
P_n_sim = zeros(size(n_plot));

for i = 1:length(n_plot)
    n = n_plot(i);
    t = ceil(2*n/3);
    p = k_A / m_A;
    P_theory = 1 - binocdf(t-1, n, p);
    P_n_theory(i) = P_theory;

    % 蒙特卡洛
    success_count = 0;
    for iter = 1:M
        votes = zeros(1,m_A);
        for expert = 1:n
            picks = randperm(m_A,k_A);
            votes(picks) = votes(picks) + 1;
        end
        success_count = success_count + (votes(1) >= t);
    end
    P_sim = success_count / M;
    P_n_sim(i) = P_sim;
end

figure;
plot(n_plot, P_n_theory, '-o', 'LineWidth',2, 'MarkerSize',8); hold on;
plot(n_plot, P_n_sim, '-s', 'LineWidth',2, 'MarkerSize',8);
xlabel('实到专家人数 n');
ylabel('单候选人当选概率');
legend('理论(Binom)', '模拟(MC)', 'Location', 'best');
title('实例A：单候选人当选概率随n增长 (m=5, k=3, α=2/3)');
grid on;

%% ----------------- 实例 B -----------------
% 参数
n_B = 15;          % 实到专家人数
m_B = 9;           % 候选人数
k_B = floor(m_B/2)+1; % 每位专家投票数
t_B = ceil(2*n_B/3);  % 当选票数阈值
s_B = floor(m_B/2);   % 推优名额

fprintf('\n=== 实例 B ===\n');
run_simulation(N, n_B, m_B, k_B, t_B, M, 'B');

% 绘制实例B中单候选人当选概率随n增长的折线图
fprintf('\n=== 实例B：单候选人当选概率随n增长 ===\n');
n_plot_B = 10:19;
P_n_theory_B = zeros(size(n_plot_B));
P_n_sim_B = zeros(size(n_plot_B));

for i = 1:length(n_plot_B)
    n = n_plot_B(i);
    t = ceil(2*n/3);
    p = k_B / m_B;
    P_theory = 1 - binocdf(t-1, n, p);
    P_n_theory_B(i) = P_theory;

    % 蒙特卡洛
    success_count = 0;
    for iter = 1:M
        votes = zeros(1,m_B);
        for expert = 1:n
            picks = randperm(m_B,k_B);
            votes(picks) = votes(picks) + 1;
        end
        success_count = success_count + (votes(1) >= t);
    end
    P_sim = success_count / M;
    P_n_sim_B(i) = P_sim;
end

figure;
plot(n_plot_B, P_n_theory_B, '-o', 'LineWidth',2, 'MarkerSize',8); hold on;
plot(n_plot_B, P_n_sim_B, '-s', 'LineWidth',2, 'MarkerSize',8);
xlabel('实到专家人数 n');
ylabel('单候选人当选概率');
legend('理论(Binom)', '模拟(MC)', 'Location', 'best');
title('实例B：单候选人当选概率随n增长 (m=9, k=5, α=2/3)');
grid on;

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

figure;
plot(n_range, P_theory_list, '-o', 'LineWidth',1.5); hold on;
plot(n_range, P_sim_list, '-s', 'LineWidth',1.5);
xlabel('实到专家人数 n');
ylabel('单候选人当选概率');
legend('理论(Binom)', '模拟(MC)');
title('灵敏度分析：出席人数 n 对当选概率的影响');
grid on;

%% ----------------- 灵敏度分析（k/m变化） -----------------
fprintf('\n=== 灵敏度分析 (k/m变化) ===\n');
n = 15; m = 5; % 固定为实例 A
k_range = 1:m; % 每位专家投票数
P_sim_list_k = zeros(size(k_range));
P_theory_list_k = zeros(size(k_range));

for i = 1:length(k_range)
    k = k_range(i);
    t = ceil(2*n/3);
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

figure;
plot(k_range/m, P_theory_list_k, '-o', 'LineWidth',1.5); hold on;
plot(k_range/m, P_sim_list_k, '-s', 'LineWidth',1.5);
xlabel('投票比例 k/m');
ylabel('单候选人当选概率');
legend('理论(Binom)', '模拟(MC)');
title('灵敏度分析：投票比例 k/m 对当选概率的影响');
grid on;

%% ----------------- 灵敏度分析（阈值比例α变化） -----------------
fprintf('\n=== 灵敏度分析 (阈值比例α变化) ===\n');
n = 15; m = 5; k = 3; % 固定为实例 A
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

figure;
plot(alpha_range, P_theory_list_alpha, '-o', 'LineWidth',1.5); hold on;
plot(alpha_range, P_sim_list_alpha, '-s', 'LineWidth',1.5);
xlabel('阈值比例 α');
ylabel('单候选人当选概率');
legend('理论(Binom)', '模拟(MC)');
title('灵敏度分析：阈值比例 α 对当选概率的影响');
grid on;


%% ----------------- 灵敏度分析（候选人数m变化） -----------------
fprintf('\n=== 灵敏度分析 (候选人数m变化) ===\n');
n = 15; alpha = 2/3; % 固定为实例 A
m_range = 3:12; % 候选人数
P_sim_list_m = zeros(size(m_range));
P_theory_list_m = zeros(size(m_range));

for i = 1:length(m_range)
    m = m_range(i);
    k = ceil(m/2); % 每位专家投票数，随m变化
    t = ceil(alpha * n);
    p = k / m;
    P_theory = 1 - binocdf(t-1, n, p);
    P_theory_list_m(i) = P_theory;

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
    P_sim_list_m(i) = P_sim;

    fprintf('m=%d, k=%d: 理论=%.4f, 模拟=%.4f\n', m, k, P_theory, P_sim);
end

figure;
plot(m_range, P_theory_list_m, '-o', 'LineWidth',1.5); hold on;
plot(m_range, P_sim_list_m, '-s', 'LineWidth',1.5);
xlabel('候选人数 m');
ylabel('单候选人当选概率');
legend('理论(Binom)', '模拟(MC)');
title('灵敏度分析：候选人数 m 对当选概率的影响');
grid on;

%% ----------------- 全局敏感性分析（m） -----------------
fprintf('\n=== 全局敏感性分析 (只分析m对P的影响，其他参数固定) ===\n');
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

%% ----------------- 全局敏感性分析（n） -----------------
fprintf('\n=== 全局敏感性分析 (只分析n对P的影响，其他参数固定) ===\n');
n_samples = 10:19;
m = 5; k = 3; alpha = 2/3;
P_theory_samples_n = zeros(size(n_samples));
P_sim_samples_n = zeros(size(n_samples));

for i = 1:length(n_samples)
    n = n_samples(i);
    t = ceil(alpha * n);
    p = k / m;
    P_theory = 1 - binocdf(t-1, n, p);
    P_theory_samples_n(i) = P_theory;

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
    P_sim_samples_n(i) = P_sim;
end

% 计算皮尔逊相关系数
R_theory_n = corr(n_samples', P_theory_samples_n', 'Type', 'Pearson');
R_sim_n = corr(n_samples', P_sim_samples_n', 'Type', 'Pearson');
fprintf('皮尔逊相关系数（理论）：R = %.4f\n', R_theory_n);
fprintf('皮尔逊相关系数（模拟）：R = %.4f\n', R_sim_n);

figure;
plot(n_samples, P_theory_samples_n, '-o', 'LineWidth',1.5); hold on;
plot(n_samples, P_sim_samples_n, '-s', 'LineWidth',1.5);
xlabel('实到专家人数 n');
ylabel('单候选人当选概率');
legend('理论(Binom)', '模拟(MC)');
title(sprintf('全局敏感性分析：n对P的影响 (R理论=%.2f, R模拟=%.2f)', R_theory_n, R_sim_n));
grid on;

%% ----------------- 全局敏感性分析（m和n联合分布） -----------------
fprintf('\n=== 全局敏感性分析 (m和n联合分布对P的影响) ===\n');
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

% 将2D矩阵展平为1D向量用于相关系数计算
m_flat = repmat(m_samples_joint', 1, length(n_samples_joint));
m_flat = m_flat(:);
n_flat = repmat(n_samples_joint, length(m_samples_joint), 1);
n_flat = n_flat(:);
P_theory_flat = P_theory_joint(:);
P_sim_flat = P_sim_joint(:);

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


%% ----------------- 通用子函数 -----------------
function run_simulation(N, n, m, k, t, M, tag)
    p = k / m;
    P_binom = 1 - binocdf(t-1, n, p);
    fprintf('单候选人理论当选概率（Binomial）：%.4f\n', P_binom);

    success_count = zeros(m,1);
    R_list = zeros(M,1);

    for iter = 1:M
        votes = zeros(1,m);
        for expert = 1:n
            picks = randperm(m,k);
            votes(picks) = votes(picks) + 1;
        end
        success_count = success_count + (votes' >= t);
        R_list(iter) = sum(votes >= t);
    end

    P_sim = success_count(1) / M;
    fprintf('单候选人蒙特卡洛当选概率：%.4f\n', P_sim);

    fprintf('当选人数分布 P(R=r):\n');
    for r = 0:m
        Pr = sum(R_list == r) / M;
        if Pr > 0
            fprintf('P(R=%d) = %.4f\n', r, Pr);
        end
    end

    figure;
    histogram(R_list, 'BinMethod','integers');
    xlabel('当选人数 R');
    ylabel('频数');
    title(sprintf('蒙特卡洛模拟：实例 %s 当选人数分布', tag));
    grid on;
end