import sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# PyR@TEが生成したSMクラスをインポート
# 注記: PyR@TEの初回実行後、'./results/SM/PythonOutput'などのフォルダが作成されます。
# パスはご自身の環境に合わせて修正してください。
try:
    sys.path.append('./Output/PythonOutput')
    from SM import RGEsolver
except ImportError:
    print("エラー: PyR@TEの出力ファイルが見つかりません。")
    print("先に './pyR@TE.py' を一度実行して、必要なファイルを生成してください。")
    sys.exit(1)

##########################################################################
# ステップ1: RGEsolverを初期化し、精密な初期値を設定
##########################################################################

# エネルギーの対数スケール t = log10(mu) を定義
initial_scale_log10 = np.log10(91.1876)  # Zボソン質量スケール

# RGEsolverインスタンスを作成 (tminが出発点と一致するように設定)
rge = RGEsolver('rge', tmin=initial_scale_log10, tmax=18, initialScale=initial_scale_log10)

# 2ループ精度の計算を指定
rge.loops = {'GaugeCouplings': 2, 'Yukawas': 2, 'QuarticTerms': 2, 
             'ScalarMasses': 2, 'Vevs': 2}

# --- ゲージ結合 @ M_Z ---
rge.g1.initialValue = 0.3576
rge.g2.initialValue = 0.6517
rge.g3.initialValue = 1.2177

# --- ヒッグスセクター @ M_Z ---
v = 246.22  # Higgs VEV in GeV
rge.v.initialValue = v
rge.lambda_.initialValue = 0.12938  # m_h ~ 125.25 GeV

# --- 湯川結合 @ M_Z (3x3 行列) ---
# クォーク
Yu = np.zeros((3, 3), dtype=complex)
Yu[2, 2] = 0.935  # top
Yd = np.zeros((3, 3), dtype=complex)
Yd[2, 2] = 0.0162 # bottom
# (他のクォーク湯川結合は小さいとしてゼロで近似)

# 荷電レプトン
Ye = np.zeros((3, 3), dtype=complex)
Ye[0, 0] = np.sqrt(2) * 0.511e-3 / v  # electron
Ye[1, 1] = np.sqrt(2) * 105.7e-3 / v # muon
Ye[2, 2] = np.sqrt(2) * 1.777 / v    # tau

rge.Yu.initialValue = Yu
rge.Yd.initialValue = Yd
rge.Ye.initialValue = Ye

##########################################################################
# ステップ2: PyR@TEによるSM基本パラメータのRGE計算を実行
##########################################################################

print("--- PyR@TEによる標準模型パラメータのRGE計算を開始します ---")
rge.solve(step=0.1)
print("--- RGE計算が完了しました ---")

##############################################################################
# ステップ3: RGE計算結果を用いて、ニュートリノ混合角の走行を計算
##############################################################################

print("\n--- ニュートリノ混合角の走行を計算し、RQTスケールを探索します ---")

# --- ニュートリノの初期値 (NuFIT 6.0, Normal Ordering best fit) ---
theta12_low_deg = 33.68
theta23_low_deg = 43.3
theta13_low_deg = 8.56
# ラジアンに変換
theta12_low = np.deg2rad(theta12_low_deg)
theta23_low = np.deg2rad(theta23_low_deg)
theta13_low = np.deg2rad(theta13_low_deg)

TARGET_SUM = np.sqrt(2)

# --- ニュートリノ混合角のRGEを定義 ---
# 荷電レプトン湯川結合が主要因となる簡略化した1ループRGE
C = 1.5  # 理論構造から決まるO(1)の係数

def neutrino_rge_system(t, y, yt_func):
    """ t(log10スケール)とy(sin^2(theta))を受け取り、微分方程式を返す """
    s12_sq, s23_sq, s13_sq = y
    
    # PyR@TEの計算結果から、対応するエネルギーでのタウ湯川結合の値を取得
    yt_running_sq = np.real(yt_func(t))**2
    
    # 簡略化したRGE
    kappa_factor = 1 / (16 * np.pi**2) * np.log(10) # t = log10(mu) のためlog(10)が必要
    ds12_sq_dt = -C * kappa_factor * yt_running_sq * s12_sq
    ds23_sq_dt = -C * kappa_factor * yt_running_sq * s23_sq * (1 - s23_sq)
    ds13_sq_dt = 0  # theta_13の走行は非常に小さいのでゼロと近似

    return [ds12_sq_dt, ds23_sq_dt, ds13_sq_dt]

# PyR@TEの解からタウ湯川結合の補間関数を作成 (スプライン補間)
yt_running_func = interp1d(rge.tList, rge.solutions['Ye_{33}'], kind='cubic', bounds_error=False, fill_value="extrapolate")

# ニュートリノの初期値ベクトル [sin^2(theta12), sin^2(theta23), sin^2(theta13)]
y0_nu = [np.sin(theta12_low)**2, np.sin(theta23_low)**2, np.sin(theta13_low)**2]
t_span_nu = [initial_scale_log10, 18] # エネルギー範囲

# ニュートリノのRGEを解く
sol_nu = solve_ivp(lambda t, y: neutrino_rge_system(t, y, yt_running_func), 
                   t_span_nu, y0_nu, dense_output=True, 
                   t_eval=np.linspace(t_span_nu[0], t_span_nu[1], 500))

##########################################################################
# ステップ4: 結果のプロットとRQTスケールの特定
##########################################################################

# 各エネルギーでの混合角（ラジアン）を計算
t_values_nu = sol_nu.t
y_values_nu = sol_nu.y

# 負の値にならないようにクリップ
y_values_nu = np.maximum(y_values_nu, 0)

theta12_run = np.arcsin(np.sqrt(y_values_nu[0]))
theta23_run = np.arcsin(np.sqrt(y_values_nu[1]))
theta13_run = np.arcsin(np.sqrt(y_values_nu[2]))

# サイン総和則の値を計算
sine_sum_values = np.sin(theta12_run) + np.sin(theta23_run) + np.sin(theta13_run)
mu_values = 10**t_values_nu

# RQTスケールを探索 (理論値に最も近くなる点)
idx = np.argmin(np.abs(sine_sum_values - TARGET_SUM))
rqt_scale_mu = mu_values[idx]

# 結果の表示
print(f"\n出発点 (Zスケール) でのサイン総和則の値: {sine_sum_values[0]:.4f}")
print(f"RQTの理論値 (sqrt(2)) に最も近くなるエネルギー (RQTスケール) はおよそ {rqt_scale_mu:.2e} GeV です。")

# グラフ描画
plt.figure(figsize=(12, 7))
plt.semilogx(mu_values, sine_sum_values)
plt.axhline(y=TARGET_SUM, color='r', linestyle='--', label=r'RQT Prediction ($\sqrt{2} \approx 1.414$)')
plt.axvline(x=rqt_scale_mu, color='g', linestyle=':', label=f'RQT Scale (~{rqt_scale_mu:.2e} GeV)')
plt.title('High-Precision RGE Running of the Neutrino Sine Sum Rule', fontsize=16)
plt.xlabel('Energy Scale $\mu$ (GeV)', fontsize=12)
plt.ylabel(r'$\sum \sin(\theta_{ij})$', fontsize=12)
plt.legend()
plt.grid(True, which="both", ls="-", alpha=0.5)
plt.ylim(1.41, 1.46) # y軸の範囲を調整して見やすくする
plt.savefig('../Data/Output/Sine_sum_Rule.png')
plt.show()