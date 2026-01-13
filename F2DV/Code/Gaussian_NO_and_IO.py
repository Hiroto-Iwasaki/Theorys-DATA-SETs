# ✅ 日本語フォント設定（Colab用）-----------------------------
#なぜか!aptに#を付けないと保存できなかったので、#を外してから使用してください。
#!apt-get -y install fonts-noto-cjk
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import os

font_path = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
font_prop = fm.FontProperties(fname=font_path)
plt.rcParams["font.family"] = font_prop.get_name()

# =========================================================
# 0) Core definitions
# =========================================================
dK = {"12": 1, "23": 1, "13": 2}
dC = {"12": 4, "23": 2, "13": 6}

def fit_p_r_from_12_23_vec(s12, s23):
    """
    log2(sinθij) = -(pΔK + rΔC)  を 12,23 の2本で解く（ベクトル版）
    """
    b12 = -np.log2(s12)
    b23 = -np.log2(s23)
    A = np.array([[dK["12"], dC["12"]],
                  [dK["23"], dC["23"]]], dtype=float)
    Ainv = np.linalg.inv(A)
    # [p, r]^T = A^{-1} [b12, b23]^T
    p = Ainv[0,0]*b12 + Ainv[0,1]*b23
    r = Ainv[1,0]*b12 + Ainv[1,1]*b23
    return p, r

def s13_base_from_pr_vec(p, r):
    expo13 = p*dK["13"] + r*dC["13"]
    return 2.0**(-expo13)

def epsilon_from_s13_vec(s13_base, s13_obs, x13_principle=np.sqrt(2.0)):
    # epsilon = sqrt2 - log2(s13_base/s13_obs)
    return x13_principle - np.log2(s13_base / s13_obs)

def x13_eff_from_epsilon_vec(eps, x13_principle=np.sqrt(2.0)):
    return x13_principle - eps

# --- geometric input ---
def mnu_geo_from_me(m_e_MeV=0.51099895, K_geo=-48):
    # eV に直してから (sqrt2)^K
    return float((m_e_MeV * 1e6) * (np.sqrt(2.0)**K_geo))

def light_masses_from_IO_vec(m3, dm21=7.5e-5, dm32_abs=2.5e-3):
    # IO: m3 given, m2^2 = m3^2 + |Δm^2_32|,  m1^2 = m2^2 - Δm^2_21
    m2 = np.sqrt(m3**2 + dm32_abs)
    m1 = np.sqrt(np.maximum(0.0, m2**2 - dm21))
    return m1, m2, m3

def light_masses_from_NO_vec(m1, dm21=7.5e-5, dm31=2.5e-3):
    # NO: m1 given, m2^2 = m1^2 + Δm^2_21, m3^2 = m1^2 + Δm^2_31
    m2 = np.sqrt(m1**2 + dm21)
    m3 = np.sqrt(m1**2 + dm31)
    return m1, m2, m3

def ci(x, lo, hi):
    x = np.asarray(x, dtype=float)
    return float(np.quantile(x, lo)), float(np.quantile(x, hi))

# =========================================================
# 1) Switch: IO / NO
# =========================================================
ORDERING = "IO"   # "IO" も可（ここだけ切替）

# =========================================================
# 2) Baseline PMNS targets（あなたのpmns_fits Uabs_eff を使う版に後で置換可）
# =========================================================
# 注意：ここは「オクタントの名前」を low/high に統一（NO/IOは ORDERING 側）
pmns_targets = {
    "low_octant":  {"sin2_12": 0.303, "sin2_23": 0.451, "sin2_13": 0.02225},
    "high_octant": {"sin2_12": 0.303, "sin2_23": 0.572, "sin2_13": 0.02203},
}

# =========================================================
# 3) Gaussian uncertainties (1σ; 編集ポイント)
# =========================================================
sigma_sin2_12 = 0.012
sigma_sin2_23 = 0.020
sigma_sin2_13 = 0.0005

sigma_dm21_frac = 0.03
sigma_dm3x_frac = 0.03     # IOなら |Δm^2_32|, NOなら Δm^2_31
sigma_mgeo_frac = 0.02     # m3_geo (IO) or m1_geo (NO) の揺らぎ

N = 200000
seed = 12345
rng = np.random.default_rng(seed)

x13_principle = np.sqrt(2.0)

# =========================================================
# 4) Baseline geo masses
# =========================================================
m_geo0 = mnu_geo_from_me()  # IOなら m3_geo, NOなら m1_geo として使う
dm21_0 = 7.5e-5
dm3x_0 = 2.5e-3

print("=== Settings ===")
print("ORDERING =", ORDERING)
print("N =", N, "seed =", seed)
print("m_geo0 [eV] =", m_geo0)

# =========================================================
# 5) Monte Carlo (vectorized)
# =========================================================
rows = []

for octant, tgt in pmns_targets.items():
    s212_0 = float(tgt["sin2_12"])
    s223_0 = float(tgt["sin2_23"])
    s213_0 = float(tgt["sin2_13"])

    # angles: Gaussian in sin^2, clipped
    s212 = np.clip(rng.normal(s212_0, sigma_sin2_12, N), 1e-9, 1-1e-9)
    s223 = np.clip(rng.normal(s223_0, sigma_sin2_23, N), 1e-9, 1-1e-9)
    s213 = np.clip(rng.normal(s213_0, sigma_sin2_13, N), 1e-9, 1-1e-9)

    s12 = np.sqrt(s212)
    s23 = np.sqrt(s223)
    s13_obs = np.sqrt(s213)

    # p,r from 12,23 (vector)
    p, r = fit_p_r_from_12_23_vec(s12, s23)

    # base s13 from p,r
    s13_base = s13_base_from_pr_vec(p, r)

    # epsilon and xeff
    eps13 = epsilon_from_s13_vec(s13_base, s13_obs, x13_principle=x13_principle)
    xeff  = x13_eff_from_epsilon_vec(eps13, x13_principle=x13_principle)

    # oscillations
    dm21 = dm21_0 * (1.0 + rng.normal(0.0, sigma_dm21_frac, N))
    dm3x = dm3x_0 * (1.0 + rng.normal(0.0, sigma_dm3x_frac, N))
    dm21 = np.clip(dm21, 1e-12, None)
    dm3x = np.clip(dm3x, 1e-12, None)

    # geo mass
    mgeo = m_geo0 * (1.0 + rng.normal(0.0, sigma_mgeo_frac, N))
    mgeo = np.clip(mgeo, 1e-12, None)

    # masses + sum mnu (ORDERING switch)
    if ORDERING.upper() == "IO":
        m1, m2, m3 = light_masses_from_IO_vec(mgeo, dm21=dm21, dm32_abs=dm3x)
    elif ORDERING.upper() == "NO":
        m1, m2, m3 = light_masses_from_NO_vec(mgeo, dm21=dm21, dm31=dm3x)
    else:
        raise ValueError("ORDERING must be 'IO' or 'NO'.")

    sum_mnu = m1 + m2 + m3

    # store
    df_oct = pd.DataFrame({
        "ordering": ORDERING.upper(),
        "octant": octant,
        "sin2_12": s212,
        "sin2_23": s223,
        "sin2_13": s213,
        "p": p, "r": r,
        "s13_base": s13_base,
        "s13_obs": s13_obs,
        "epsilon13_l": eps13,
        "x13_eff_l": xeff,
        "dm21": dm21,
        "dm3x": dm3x,
        "mgeo": mgeo,
        "m1": m1, "m2": m2, "m3": m3,
        "sum_mnu": sum_mnu,
    })
    rows.append(df_oct)

df = pd.concat(rows, ignore_index=True)

# =========================================================
# 6) Summary print
# =========================================================
print("\n=== Gaussian MC summary ===")
for octant in pmns_targets.keys():
    sub = df[df["octant"]==octant]
    eps = sub["epsilon13_l"].values
    xef = sub["x13_eff_l"].values
    smn = sub["sum_mnu"].values
    corr = np.corrcoef(eps, smn)[0,1]
    print(f"\n--- {ORDERING.upper()} {octant} ---")
    print("epsilon13_l baseline ≈", float(np.median(eps)))
    print("epsilon13_l mean/std =", (float(np.mean(eps)), float(np.std(eps))))
    print("epsilon13_l 68% CI   =", ci(eps, 0.16, 0.84))
    print("epsilon13_l 95% CI   =", ci(eps, 0.025, 0.975))
    print("x13_eff_l   68% CI   =", ci(xef, 0.16, 0.84))
    print("x13_eff_l   95% CI   =", ci(xef, 0.025, 0.975))
    print("sum mnu mean/std     =", (float(np.mean(smn)), float(np.std(smn))))
    print("corr(eps, sum mnu)   =", float(corr))

# =========================================================
# 7) Plots (overlay low+high in ONE figure)
# =========================================================
outdir = "GAV_outputs"
os.makedirs(outdir, exist_ok=True)

# 7-1) epsilon hist (overlay)
plt.figure(figsize=(10,5))
for octant in pmns_targets.keys():
    sub = df[df["octant"]==octant]
    plt.hist(sub["epsilon13_l"], bins=140, alpha=0.55, density=True, label=f"{ORDERING.upper()}_{octant}")
plt.xlabel(r"$\epsilon^{(\ell)}_{13}$")
plt.ylabel("pdf")
plt.title(rf"{ORDERING.upper()} Gaussian MC: distribution of $\epsilon^{{(\ell)}}_{{13}}$ (low+high overlay)")
plt.legend()
hist_path = f"{outdir}/gaussian_mc_{ORDERING.upper()}_dist_epsilon_l_13.png"
plt.savefig(hist_path, dpi=200, bbox_inches="tight")
plt.show()
print("Saved:", hist_path)

# 7-2) eps vs sum mnu (overlay)
plt.figure(figsize=(10,6))
for octant in pmns_targets.keys():
    sub = df[df["octant"]==octant]
    plt.scatter(sub["epsilon13_l"], sub["sum_mnu"], s=1.5, alpha=0.15, label=f"{ORDERING.upper()}_{octant}")
plt.xlabel(r"$\epsilon^{(\ell)}_{13}$")
plt.ylabel(r"$\sum m_\nu$ [eV]")
plt.title(rf"{ORDERING.upper()} Gaussian MC: $\epsilon^{{(\ell)}}_{{13}}$ vs $\sum m_\nu$ (low+high overlay)")
plt.legend()
scat_path = f"{outdir}/gaussian_mc_{ORDERING.upper()}_epsilon_l_13_vs_sum_mnu.png"
plt.savefig(scat_path, dpi=200, bbox_inches="tight")
plt.show()
print("Saved:", scat_path)

# =========================================================
# 8) Save CSV
# =========================================================
csv_path = f"{outdir}/pmns_epsilon13_gaussian_mc_{ORDERING.upper()}.csv"
df.to_csv(csv_path, index=False)
print("\nSaved:", csv_path)

# =========================================================
# （数式の項の意味と役割）※このセルの主役だけまとめ
# =========================================================
print("\n(数式の項の意味と役割)")
print("・log2(sinθij) = -(pΔK + rΔC)：12,23 から p,r を決め、13 の基準値 s13_base を生成する骨格。")
print("・s13_base：生成則だけで予言される sinθ13（13補正なしの基準）。")
print("・s13_obs：入力フィットの sinθ13（観測側）。")
print("・epsilon13^(l) = √2 - log2(s13_base/s13_obs)：13だけの残差（原理 √2 からのズレとして保存）。")
print("・x13_eff^(l) = √2 - epsilon13^(l)：実効指数（論文の x13=√2-ε と同じ形）。")
print("・m_geo：幾何学入力で固定する基準質量（IOなら m3、NOなら m1）。")
print("・sum mnu：m1+m2+m3（宇宙論制限と直結、位相非依存）。")