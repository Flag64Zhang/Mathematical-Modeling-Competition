% 投票选拔问题 - 不同到场专家数n及不同候选人数m对候选人当选的总人数R的影响
% 文件名: question3.m

clear; clc; close all;

n_range = 10:19;      % 专家到场人数
m_range = 3:9;        % 候选人数
alpha = 2/3;          % 阈值比例
k_func = @(m) ceil(m/2); % 每位专家投票数
M = 10000;            % 蒙特卡洛次数

max_R = max(m_range); % 最大可能当选人数
P_R = zeros(length(n_range), length(m_range), max_R+1); % 存储概率分布

for i = 1:length(n_range)
    n = n_range(i);
    t = ceil(alpha * n);
    for j = 1:length(m_range)
        m = m_range(j);
        k = k_func(m);
        R_list = zeros(M,1);
        for iter = 1:M
            votes = zeros(1,m);
            for expert = 1:n
                picks = randperm(m, k);
                votes(picks) = votes(picks) + 1;
            end
            R_list(iter) = sum(votes >= t);
        end
        % 统计R=0,1,...,m的概率
        for r = 0:m
            P_R(i,j,r+1) = sum(R_list == r) / M;
        end
    end
end

% 可视化：以n=15为例，画m随R分布的热力图
n_idx = find(n_range==15);
figure;
imagesc(m_range, 0:max_R, squeeze(P_R(n_idx,:,:)));
xlabel('候选人数 m');
ylabel('当选人数 R');
title('n=15时，不同m下R的概率分布');
colorbar;
set(gca,'YDir','normal');

% 可视化：以m=5为例，画n随R分布的热力图
m_idx = find(m_range==5);
figure;
imagesc(n_range, 0:max_R, squeeze(P_R(:,m_idx,:)));
xlabel('专家人数 n');
ylabel('当选人数 R');
title('m=5时，不同n下R的概率分布');
colorbar;
set(gca,'YDir','normal');

%% 分组柱状图展示当选总人数R的概率分布（draft.m功能）
m_list = [5, 7, 9];
n_list_bar = 10:2:18;
iterations = 10000;
prob_R = cell(length(m_list), length(n_list_bar));

function prob = calculate_R_prob(n, m, iter)
    k = floor(m/2) + 1;
    threshold = ceil(2 * n / 3);
    r_counts = zeros(1, m + 1);
    for i = 1:iter
        votes = zeros(n, m);
        for j = 1:n
            selected = randperm(m, k);
            votes(j, selected) = 1;
        end
        total_votes = sum(votes, 1);
        R = sum(total_votes >= threshold);
        r_counts(R + 1) = r_counts(R + 1) + 1;
    end
    prob = r_counts / iter;
end

for m_idx = 1:length(m_list)
    m = m_list(m_idx);
    for n_idx = 1:length(n_list_bar)
        n = n_list_bar(n_idx);
        prob_R{m_idx, n_idx} = calculate_R_prob(n, m, iterations);
        fprintf('已完成：候选人数m=%d，到场专家数n=%d\n', m, n);
    end
end

figure('Name','当选总人数R的概率分布','Position',[100 100 1200 700]);
bar_width = 0.15;
gap = 0.05;
for m_idx = 1:length(m_list)
    m = m_list(m_idx);
    for n_idx = 1:length(n_list_bar)
        n = n_list_bar(n_idx);
        prob = prob_R{m_idx, n_idx};
        max_R_bar = find(prob > 0, 1, 'last') - 1;
        if isempty(max_R_bar), max_R_bar = 0; end
        x = (0:max_R_bar) + (m_idx - 1)*(length(n_list_bar)*(bar_width + gap)) + (n_idx - 1)*bar_width;
        bar(x, prob(1:max_R_bar+1), bar_width, ...
            'DisplayName', sprintf('m=%d人, n=%d人', m, n));
        hold on;
    end
end
xlabel('当选总人数R', 'FontSize',12);
ylabel('概率', 'FontSize',12);
title('不同候选人数、到场专家数下，当选总人数R的概率分布', 'FontSize',14, 'FontWeight','bold');
xticks(0:max(m_list));
xticklabels(0:max(m_list));
legend('Location','eastoutside');
grid on;
box on;

%% 关键场景数据输出

disp('===== 关键场景下当选总人数R的概率（R=1,2,3,...）=====');
for m_idx = 1:length(m_list)
    m = m_list(m_idx);
    for n_idx = 1:length(n_list_bar)
        n = n_list_bar(n_idx);
        prob = prob_R{m_idx, n_idx};
        fprintf('\n候选人数m=%d，到场专家数n=%d：\n', m, n);
        for r = 1:min(5, m)
            fprintf('  R=%d的概率：%.4f\n', r, prob(r+1));
        end
        if m > 5
            fprintf('  R>5的概率：%.4f\n', sum(prob(7:end)));
        end
    end
end
