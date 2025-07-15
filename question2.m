% 投票选拔问题 - 实例A+B蒙特卡洛模拟
% 文件名: question2.m

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
run_simulation(N, n_A, m_A, k_A, t_A, M, 'A', s_A);

% 绘制实例A中单候选人当选概率随n增长的折线图
fprintf('实例A：单候选人当选概率随n增长折线图绘制完毕\n');
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
title('实例A：单候选人当选概率随n增长变化情况 (m=5, k=3, α=2/3)');
grid on;

%% ----------------- 实例 B -----------------
% 参数
n_B = 15;          % 实到专家人数
m_B = 9;           % 候选人数
k_B = floor(m_B/2)+1; % 每位专家投票数
t_B = ceil(2*n_B/3);  % 当选票数阈值
s_B = floor(m_B/2);   % 推优名额

fprintf('\n=== 实例 B ===\n');
run_simulation(N, n_B, m_B, k_B, t_B, M, 'B', s_B);

% 绘制实例B中单候选人当选概率随n增长的折线图
fprintf('实例B：单候选人当选概率随n增长折线图绘制完毕\n');
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
title('实例B：单候选人当选概率随n增长变化情况 (m=9, k=5, α=2/3)');
grid on;


%% ----------------- 通用子函数 -----------------
function run_simulation(N, n, m, k, t, M, tag, s)
    p = k / m;
    P_binom = 1 - binocdf(t-1, n, p);
    fprintf('单候选人理论当选概率（Binomial）：%.4f\n', P_binom);

    success_count = zeros(m,1);
    R_list = zeros(M,1);
    final_success_count = zeros(m,1); % 统计最终被推优的次数

    for iter = 1:M
        votes = zeros(1,m);
        for expert = 1:n
            picks = randperm(m,k);
            votes(picks) = votes(picks) + 1;
        end
        success_count = success_count + (votes >= t);
        R_list(iter) = sum(votes >= t);

        % 推优名额限制
        winners = find(votes >= t);
        if length(winners) > s
            [~, idx] = sort(votes(winners), 'descend');
            winners = winners(idx(1:s));
        end
        final_success_count(winners) = final_success_count(winners) + 1;
    end

    P_sim = success_count(1) / M;
    P_final = final_success_count(1) / M;
    fprintf('单候选人蒙特卡洛当选概率（达标）：%.4f\n', P_sim);
    fprintf('单候选人蒙特卡洛最终推优概率（含名额限制）：%.4f\n', P_final);

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