% 投票选拔问题 可调节参数 GUI
% 文件名: vote_gui.m

function vote_gui
    % 创建界面
    fig = uifigure('Name','投票选拔模拟 GUI','Position',[500 300 400 400]);

    % 标签 & 输入框
    uilabel(fig,'Position',[30 340 150 30],'Text','实到专家人数 n:');
    nEdit = uieditfield(fig,'numeric','Position',[200 345 100 22],'Value',15);

    uilabel(fig,'Position',[30 300 150 30],'Text','候选人数 m:');
    mEdit = uieditfield(fig,'numeric','Position',[200 305 100 22],'Value',5);

    uilabel(fig,'Position',[30 260 150 30],'Text','每位专家投票数 k:');
    kEdit = uieditfield(fig,'numeric','Position',[200 265 100 22],'Value',3);

    uilabel(fig,'Position',[30 220 150 30],'Text','蒙特卡洛次数 M:');
    MEdit = uieditfield(fig,'numeric','Position',[200 225 100 22],'Value',10000);

    uilabel(fig,'Position',[30 180 150 30],'Text','阈值比例 α (0~1):');
    alphaEdit = uieditfield(fig,'numeric','Position',[200 185 100 22],'Value',2/3);

    % 按钮
    btn = uibutton(fig,'Text','运行模拟',...
        'Position',[150 120 100 40],...
        'ButtonPushedFcn', @(btn,event) runSimulation);

    % 输出区
    txt = uitextarea(fig,'Position',[30 20 340 80],'Editable','off');

    % 内部函数
    function runSimulation
        n = nEdit.Value;
        m = mEdit.Value;
        k = kEdit.Value;
        M = MEdit.Value;
        alpha = alphaEdit.Value;

        t = ceil(alpha * n);
        p = k / m;

        % 理论值
        P_theory = 1 - binocdf(t-1, n, p);

        % 蒙特卡洛
        success = 0;
        R_list = zeros(M,1);
        for iter = 1:M
            votes = zeros(1,m);
            for expert = 1:n
                picks = randperm(m,k);
                votes(picks) = votes(picks) + 1;
            end
            success = success + (votes(1) >= t);
            R_list(iter) = sum(votes >= t);
        end
        P_sim = success / M;

        % 显示结果
        result = sprintf(['参数: n=%d, m=%d, k=%d, α=%.2f, t=%d\n' ...
                          '单候选人理论值：%.4f\n' ...
                          '单候选人模拟值：%.4f\n'], ...
                          n, m, k, alpha, t, P_theory, P_sim);
        txt.Value = result;

        % 绘图
        figure(99); clf;
        histogram(R_list, 'BinMethod','integers');
        xlabel('当选人数 R');
        ylabel('频数');
        title('当选人数分布');
        grid on;
    end
end
