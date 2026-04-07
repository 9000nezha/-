import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

plt.rcParams['font.family'] = ['sans-serif']
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'DejaVu Sans']

script_dir = os.path.dirname(os.path.abspath(__file__))

df = pd.read_csv(os.path.join(script_dir, 'results_table.csv'))
print("数据列:", df.columns.tolist())
print("数据行数:", len(df))
print("\n前10行数据:")
print(df.head(10))

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

ax1 = axes[0, 0]
ax1.plot(df['t_ms'], df['p_MPa'], 'b-', linewidth=1.5)
ax1.set_xlabel('t / ms', fontsize=12)
ax1.set_ylabel('p / MPa', fontsize=12)
ax1.set_title('图2-8 预测的57mm高射炮 p-t 曲线', fontsize=14)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, 7)
ax1.set_ylim(0, 100)

ax2 = axes[0, 1]
ax2.plot(df['l_dm'], df['p_MPa'], 'b-', linewidth=1.5)
ax2.set_xlabel('l / dm', fontsize=12)
ax2.set_ylabel('p / MPa', fontsize=12)
ax2.set_title('图2-9 预测的57mm高射炮 p-l 曲线', fontsize=14)
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0, 36)
ax2.set_ylim(0, 100)

ax3 = axes[1, 0]
ax3.plot(df['t_ms'], df['v_m_s'], 'r-', linewidth=1.5)
ax3.set_xlabel('t / ms', fontsize=12)
ax3.set_ylabel(r'v / (m·s$^{-1}$)', fontsize=12)
ax3.set_title('图2-10 预测的57mm高射炮 v-t 曲线', fontsize=14)
ax3.grid(True, alpha=0.3)
ax3.set_xlim(0, 7)
ax3.set_ylim(0, 400)

ax4 = axes[1, 1]
ax4.plot(df['l_dm'], df['v_m_s'], 'r-', linewidth=1.5)
ax4.set_xlabel('l / dm', fontsize=12)
ax4.set_ylabel(r'v / (m·s$^{-1}$)', fontsize=12)
ax4.set_title('图2-11 预测的57mm高射炮 v-l 曲线', fontsize=14)
ax4.grid(True, alpha=0.3)
ax4.set_xlim(0, 36)
ax4.set_ylim(0, 600)

plt.tight_layout()

save_path = os.path.join(script_dir, 'curves_all.png')
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"曲线图已保存到: {save_path}")

plt.show()