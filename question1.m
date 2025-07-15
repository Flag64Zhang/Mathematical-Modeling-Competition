% 投票选拔问题概率模拟与可视化
% 文件名: question1.m

clear; clc; close all;

%% 1. 参数定义（根据题目设定）
params_A = struct('m',5,'k',3,'s',2);  % m:候选人数，k:每人可投票数，s:推优名额
params_B = struct('m',9,'k',6,'s',5);

n_list = 10:19;  % 专家到场人数
iterations = 10000;  % 蒙特卡洛模拟次数

%% 2. 概率计算函数（每位候选人）
function [prob_vec, threshold] = calc_prob_each(n, m, k, iter)
    threshold = ceil(2*n/3);
    count = zeros(m,1);
    for i = 1:iter
        votes = zeros(1, m);
        for j = 1:n
            idx = randperm(m, k);
            votes(idx) = votes(idx) + 1;
        end
        count = count + (votes >= threshold)';
    end
    prob_vec = count / iter;
end

%% 3. 蒙特卡洛模拟
P_A = zeros(params_A.m, length(n_list));
P_B = zeros(params_B.m, length(n_list));
thresholds = zeros(size(n_list));

for i = 1:length(n_list)
    n = n_list(i);
    [P_A(:,i), thresholds(i)] = calc_prob_each(n, params_A.m, params_A.k, iterations);
    [P_B(:,i), ~] = calc_prob_each(n, params_B.m, params_B.k, iterations);
end

%% 4. 可视化结果
figure('Name','实例A：每位候选人达标概率随n变化','Position',[100 100 800 600]);
plot(n_list, P_A', '-o', 'LineWidth', 1.5);
xlabel('实到专家人数 n');
ylabel('P(票数≥2/3n)');
title('实例A：每位候选人达标概率随n变化');
legend(arrayfun(@(x) sprintf('候选人%d',x), 1:params_A.m, 'UniformOutput', false));
grid on;

figure('Name','实例B：每位候选人达标概率随n变化','Position',[100 100 800 600]);
plot(n_list, P_B', '-o', 'LineWidth', 1.5);
xlabel('实到专家人数 n');
ylabel('P(票数≥2/3n)');
title('实例B：每位候选人达标概率随n变化');
legend(arrayfun(@(x) sprintf('候选人%d',x), 1:params_B.m, 'UniformOutput', false));
grid on;

%% 5 阈值变化图
figure('Name','阈值变化图','Position',[100 100 800 400]);
linear_threshold = (2/3)*n_list;  % 理论线性阈值（非取整）
plot(n_list, linear_threshold, 'r--', 'LineWidth',1.5, 'DisplayName','理论线性阈值（2/3n）'); hold on;
plot(n_list, thresholds, 'd-', 'Color',[0.5 0.5 0.5], 'LineWidth',1.5, 'MarkerSize',6, 'DisplayName','实际阈值（ceil(2/3n)）');
xlabel('实到专家人数 n', 'FontSize',12);
ylabel('票数阈值', 'FontSize',12);
title('阈值变化图', 'FontSize',14);
grid on; legend('Location','best');

%% 6. 分票效应分析（每位候选人）
m_list = 5:2:13;  % 5,7,9,11,13人候选
n_fixed = 14;
P_m = zeros(max(m_list), length(m_list));
for i = 1:length(m_list)
    m = m_list(i);
    k = floor(m/2) + 1;
    P_m(1:m,i) = calc_prob_each(n_fixed, m, k, iterations);
end
figure('Name','分票效应分析','Position',[200 200 800 600]);
for i = 1:length(m_list)
    plot(1:m_list(i), P_m(1:m_list(i),i), '-o', 'LineWidth', 1.5, 'DisplayName', sprintf('m=%d',m_list(i)));
    hold on;
end
xlabel('候选人编号');
ylabel('P(票数≥2/3n)');
title('固定到场14人时，不同候选人达标概率（分票效应）');
grid on; legend('show');

%% 7. 结果输出（关键数据）
disp('实例A和B在不同到场人数下的每位候选人达标概率：');
for i = 1:params_A.m
    fprintf('实例A 候选人%d: %s\n', i, mat2str(P_A(i,:),3));
end
for i = 1:params_B.m
    fprintf('实例B 候选人%d: %s\n', i, mat2str(P_B(i,:),3));
end

% 输出为表格形式，便于查阅
T_A = array2table(P_A, 'VariableNames', compose('n%d', n_list), 'RowNames', compose('A_候选人%d', 1:params_A.m));
T_B = array2table(P_B, 'VariableNames', compose('n%d', n_list), 'RowNames', compose('B_候选人%d', 1:params_B.m));
disp('实例A概率表：'); disp(T_A);
disp('实例B概率表：'); disp(T_B);

% 分票效应分析结果输出
fprintf('\n分票效应分析（n=%d时，不同m下每位候选人达标概率）：\n', n_fixed);
for i = 1:length(m_list)
    m = m_list(i);
    fprintf('m=%d: ', m);
    fprintf('%.3f ', P_m(1:m,i));
    fprintf('\n');
end