import matplotlib.pyplot as plt
import numpy as np

# ----- 設定 -----
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['figure.dpi'] = 300

# ==========================================
# 図2: G2ルート系とフレーバー起源 (G2 Roots & Flavor Origin)
# ==========================================
def plot_g2_flavor_origin():
    fig, ax = plt.subplots(figsize=(8, 8))

    # --- G2ルート定義 ---
    # Short roots (長さ 1) -> Down-type
    angles_short = np.arange(0, 360, 60) * np.pi / 180

    # Long roots (長さ sqrt(3)) -> Up-type
    angles_long = np.arange(30, 390, 60) * np.pi / 180
    r_long = np.sqrt(3)

    # --- プロット関数 ---
    def plot_vector(angle, length, color, style, width, label=None):
        x = length * np.cos(angle)
        y = length * np.sin(angle)
        ax.arrow(0, 0, x, y, head_width=0.08, head_length=0.1,
                 fc=color, ec=color, linestyle=style, linewidth=width, length_includes_head=True, alpha=0.9)
        if label:
            # 凡例用にダミーの点を打つ
            ax.plot([], [], '-', color=color, linewidth=width, label=label)

    # --- 描画実行 ---

    # 1. Short Roots (Down-type)
    for i, ang in enumerate(angles_short):
        lbl = 'Short Roots ($|v|^2=1$)\n$\\rightarrow$ Down-type Quarks ($d,s,b$)' if i == 0 else None
        # 少し細い線、色はマゼンタ系（ダウンっぽいイメージ）
        plot_vector(ang, 1.0, '#C71585', '-', 1.5, lbl)

    # 2. Long Roots (Up-type)
    for i, ang in enumerate(angles_long):
        lbl = 'Long Roots ($|v|^2=3$)\n$\\rightarrow$ Up-type Quarks ($u,c,t$)' if i == 0 else None
        # 太い線、色は濃い赤（アップ/トップの重いイメージ）
        plot_vector(ang, r_long, '#8B0000', '-', 2.5, lbl)

    # 3. Center (Leptons)
    ax.plot(0, 0, 'o', color='blue', markersize=12, label='Singlet (Center)\n$\\rightarrow$ Leptons ($e,\mu,\\tau$)', zorder=10)

    # --- 解説テキスト (注釈) ---
    # Leptons
    ax.annotate('Leptons\n(No Friction)', xy=(0, 0), xytext=(-0.8, 0.2),
                arrowprops=dict(facecolor='blue', shrink=0.05, width=1, headwidth=5),
                fontsize=11, color='blue', fontweight='bold')

    # Top Quark suggestion
    ax.text(1.3, 1.3, r'Top Quark Direction?', fontsize=10, color='#8B0000', style='italic')

    # --- グラフの体裁 ---
    ax.set_xlim(-2.2, 2.2)
    ax.set_ylim(-2.2, 2.2)
    ax.set_aspect('equal')
    ax.grid(True, linestyle=':', alpha=0.6)

    # タイトルと凡例
    ax.set_title(r'Vacuum Texture: $G_2$ Root System & Flavor Split', fontsize=16, pad=15)

    # 凡例を枠外に出すか、邪魔にならない位置に
    ax.legend(loc='upper left', fontsize=10, framealpha=0.9, edgecolor='gray')

    plt.tight_layout()
    plt.savefig('RQT_G2_Flavor_Roots.png')
    plt.show()

plot_g2_flavor_origin()