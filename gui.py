import tkinter as tk
from tkinter import ttk, Text
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binom
import random
import matplotlib

# 设置matplotlib支持中文显示
plt.rcParams["font.family"] = ["SimHei", "WenQuanYi Micro Hei", "Heiti TC"]
plt.rcParams["axes.unicode_minus"] = False  # 解决负号显示问题

def runSimulation():
    # 获取输入框的值
    n = int(n_entry.get())
    m = int(m_entry.get())
    M = int(M_entry.get())
    alpha = float(alpha_entry.get())

    # 自动设置 k 和 s
    k = m // 2 + 1
    s = m // 2

    t = int(np.ceil(alpha * n))
    p = k / m

    # 理论值计算
    P_theory = 1 - binom.cdf(t - 1, n, p)

    # 蒙特卡洛模拟
    success = 0
    R_list = np.zeros(M)
    final_success = 0
    for iter in range(M):
        votes = np.zeros(m, dtype=int)
        for expert in range(n):
            picks = random.sample(range(m), k)
            votes[picks] += 1
        success += 1 if votes[0] >= t else 0
        R_list[iter] = np.sum(votes >= t)

        # 推优名额限制相关逻辑
        winners = np.where(votes >= t)[0]
        if len(winners) > s:
            # 按票数降序排序，取前 s 个
            winners = winners[np.argsort(votes[winners])[::-1][:s]]
        final_success += 1 if 0 in winners else 0  # 假设关注候选人 1（索引 0 ）

    P_sim = success / M
    P_final = final_success / M

    # 显示结果到文本框
    result = (
        f"参数: n={n}, m={m}, k={k}, s={s}, α={alpha:.2f}, t={t}, M={M}\n"
        f"单候选人理论值：{P_theory:.4f}\n"
        f"单候选人模拟值（达标）：{P_sim:.4f}\n"
        f"单候选人最终推优概率：{P_final:.4f}\n"
    )
    txt.delete(1.0, tk.END)
    txt.insert(tk.END, result)

    # 绘图
    plt.figure(99)
    plt.clf()
    plt.hist(R_list, bins=np.arange(int(R_list.min()), int(R_list.max()) + 2), edgecolor='black')
    plt.xlabel('当选人数 R')
    plt.ylabel('频数')
    plt.title('当选人数分布')
    plt.grid(True)
    plt.show()

# 创建主窗口
root = tk.Tk()
root.title("投票选拔模拟 GUI")
root.geometry("500x600")

# 设置ttk控件的字体样式
style = ttk.Style()
style.configure("TLabel", font=("SimHei", 10))
style.configure("TButton", font=("SimHei", 10))
style.configure("TEntry", font=("SimHei", 10))

# 创建标签和输入框
ttk.Label(root, text="实到专家人数 n:").grid(row=0, column=0, padx=10, pady=10)
n_entry = ttk.Entry(root)
n_entry.grid(row=0, column=1, padx=10, pady=10)
n_entry.insert(0, "15")

ttk.Label(root, text="候选人数 m:").grid(row=1, column=0, padx=10, pady=10)
m_entry = ttk.Entry(root)
m_entry.grid(row=1, column=1, padx=10, pady=10)
m_entry.insert(0, "5")

ttk.Label(root, text="蒙特卡洛次数 M:").grid(row=2, column=0, padx=10, pady=10)
M_entry = ttk.Entry(root)
M_entry.grid(row=2, column=1, padx=10, pady=10)
M_entry.insert(0, "10000")

ttk.Label(root, text="阈值比例 α (0~1):").grid(row=3, column=0, padx=10, pady=10)
alpha_entry = ttk.Entry(root)
alpha_entry.grid(row=3, column=1, padx=10, pady=10)
alpha_entry.insert(0, "0.67")  # 对应 2/3

# 创建运行按钮
run_btn = ttk.Button(root, text="运行模拟", command=runSimulation)
run_btn.grid(row=4, column=0, columnspan=2, pady=20)

# 创建文本框用于显示结果，设置支持中文的字体
txt = Text(root, width=50, height=10, font=("SimHei", 10))
txt.grid(row=5, column=0, columnspan=2, padx=10, pady=10)

# 启动主循环
root.mainloop()