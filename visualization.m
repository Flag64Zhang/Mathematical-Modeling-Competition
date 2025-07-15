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
imagesc(m_range, 0:max_R, squeeze(P_R(n_idx,:,:))');
xlabel('候选人数 m');
ylabel('当选人数 R');
title('n=15时，不同m下R的概率分布');
colorbar;
set(gca,'YDir','normal');

% 也可以画n随R分布的热力图
m_idx = find(m_range==5);
figure;
imagesc(n_range, 0:max_R, squeeze(P_R(:,m_idx,:))');
xlabel('专家人数 n');
ylabel('当选人数 R');
title('m=5时，不同n下R的概率分布');
colorbar;
set(gca,'YDir','normal');




%% ----------------- n减小时，单候选人得票≥2/3n的概率通常降低可视化验证 -----------------
fprintf('\n=== n减小时，单候选人得票≥2/3n的概率) ===\n');
n_range = 6:18;
m = 5; k = 3; alpha = 2/3;
P_theory = zeros(size(n_range));
for i = 1:length(n_range)
    n = n_range(i);
    t = ceil(alpha * n);
    p = k / m;
    P_theory(i) = 1 - binocdf(t-1, n, p);
end
plot(n_range, P_theory, '-o');
xlabel('实到专家人数 n');
ylabel('P(票数≥2/3n)');
title('n减小时单候选人达标概率');
grid on;
