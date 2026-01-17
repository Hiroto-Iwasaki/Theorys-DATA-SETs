# =========================================================
# 【日本語】NO/IO の Gaussian MC：ε13^(ℓ)・Σmν・mβ を推定し、
#          代表質量を用いた位相スキャンで mββ（0νββ）の予言帯を要約する
#          （※このファイルは Casas–Ibarra による Yν 生成/直交性検証は行いません）
# ---------------------------------------------------------
# 【English】Gaussian MC for NO/IO: estimate ε13^(ℓ), Σmν, and mβ,
#           and summarize the predicted mββ (0νββ) band via a Majorana-phase scan
#           using representative masses
#           (*This script does NOT perform Casas–Ibarra Yν construction / R-orthogonality checks.)
# =========================================================
#日本語フォント設定（Colab用）-----------------------------
#なぜか!aptに#を付けないと保存できなかったので、#を外してから使用してください。
#!apt-get -y install fonts-noto-cjk
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import os

font_path = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
fm.fontManager.addfont(font_path)  # ← これで findfont 警告が出にくくなります
font_prop = fm.FontProperties(fname=font_path)
plt.rcParams["font.family"] = font_prop.get_name()
plt.rcParams["axes.unicode_minus"] = False

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
    return float((m_e_MeV * 1e6) * (np.sqrt(2.0)**K_geo))  # eV

def light_masses_from_IO_vec(m3, dm21=7.5e-5, dm32_abs=2.5e-3):
    m2 = np.sqrt(m3**2 + dm32_abs)
    m1 = np.sqrt(np.maximum(0.0, m2**2 - dm21))
    return m1, m2, m3

def light_masses_from_NO_vec(m1, dm21=7.5e-5, dm31=2.5e-3):
    m2 = np.sqrt(m1**2 + dm21)
    m3 = np.sqrt(m1**2 + dm31)
    return m1, m2, m3

def ci(x, lo, hi):
    x = np.asarray(x, dtype=float)
    return float(np.quantile(x, lo)), float(np.quantile(x, hi))

# =========================================================
# 1) Switch: IO / NO
# =========================================================
ORDERING = "NO"   # "IO" も可

# =========================================================
# 2) Baseline PMNS targets（sin^2）
#    ※ここは「角度フィット入力」(ガウス揺らぎ元)として使います
# =========================================================
pmns_targets = {
    "low_octant":  {"sin2_12": 0.303, "sin2_23": 0.451, "sin2_13": 0.02225},
    "high_octant": {"sin2_12": 0.303, "sin2_23": 0.572, "sin2_13": 0.02203},
}

# =========================================================
# 2.5) (ADD) e行 |U_ei|（m_beta, m_bb 用）
#    ※あなたの pmns_fits[PMNS_KEY]['Uabs_eff'] の e行と一致させる前提
# =========================================================
Ue_abs_targets = {
    "low_octant":  np.array([0.825525, 0.544296, 0.149164], dtype=float),  # NO_low_octant
    "high_octant": np.array([0.825618, 0.544357, 0.148425], dtype=float),  # NO_high_octant
}

# =========================================================
# 3) Gaussian uncertainties (1σ)
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
# 4) Baseline geo masses & oscillation central values
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

    # (ADD) m_beta per-row (phase independent) using fixed |U_ei| per octant
    Ue_abs = Ue_abs_targets[octant]
    w = (Ue_abs**2).astype(float)  # |U_ei|^2
    m_beta = np.sqrt(w[0]*(m1**2) + w[1]*(m2**2) + w[2]*(m3**2))

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
        "m_beta": m_beta,
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
    mb  = sub["m_beta"].values
    corr = np.corrcoef(eps, smn)[0,1]
    print(f"\n--- {ORDERING.upper()} {octant} ---")
    print("epsilon13_l baseline ≈", float(np.median(eps)))
    print("epsilon13_l mean/std =", (float(np.mean(eps)), float(np.std(eps))))
    print("epsilon13_l 68% CI   =", ci(eps, 0.16, 0.84))
    print("epsilon13_l 95% CI   =", ci(eps, 0.025, 0.975))
    print("x13_eff_l   68% CI   =", ci(xef, 0.16, 0.84))
    print("x13_eff_l   95% CI   =", ci(xef, 0.025, 0.975))
    print("sum mnu mean/std     =", (float(np.mean(smn)), float(np.std(smn))))
    print("m_beta mean/std      =", (float(np.mean(mb)),  float(np.std(mb))))
    print("corr(eps, sum mnu)   =", float(corr))

# =========================================================
# 6.5) (ADD) NO: m_bb phase scan using representative masses per octant
#    ※あなたの LaTeX の「位相スキャンで予言帯」に相当（軽量で論文向き）
# =========================================================
NPHASE = 300_000
seed_phase = 20260117
rng_phase = np.random.default_rng(seed_phase)

def mbb_phase_scan_from_rep_masses(m_rep, Ue_abs, NPHASE, rng):
    """
    m_bb = | m1|Ue1|^2 + m2|Ue2|^2 e^{i alpha21} + m3|Ue3|^2 e^{i alpha31} |
    """
    w = (Ue_abs**2).astype(float)

    alpha21 = rng.uniform(0.0, 2*np.pi, NPHASE)
    alpha31 = rng.uniform(0.0, 2*np.pi, NPHASE)

    mbb = np.abs(
        m_rep[0]*w[0]
        + m_rep[1]*w[1]*np.exp(1j*alpha21)
        + m_rep[2]*w[2]*np.exp(1j*alpha31)
    )
    return mbb

print("\n=== m_bb phase scan (using representative masses; NPHASE=%d) ===" % NPHASE)
mbb_summary = {}

for octant in pmns_targets.keys():
    sub = df[df["octant"]==octant]

    # 代表値（論文と整合させやすいのは mean; 必要なら median に変更可）
    m1_rep = float(np.mean(sub["m1"]))
    m2_rep = float(np.mean(sub["m2"]))
    m3_rep = float(np.mean(sub["m3"]))
    m_rep = np.array([m1_rep, m2_rep, m3_rep], dtype=float)

    Ue_abs = Ue_abs_targets[octant]
    w = (Ue_abs**2).astype(float)

    # 位相非依存：m_beta（代表質量で算出）
    m_beta_rep = float(np.sqrt(np.sum(w*(m_rep**2))))

    # 位相スキャン：m_bb
    mbb = mbb_phase_scan_from_rep_masses(m_rep, Ue_abs, NPHASE, rng_phase)

    out = {
        "m_rep": m_rep,
        "sum_mnu_rep": float(np.sum(m_rep)),
        "Ue_abs": Ue_abs,
        "m_beta_rep": m_beta_rep,
        "mbb_min": float(np.min(mbb)),
        "mbb_max": float(np.max(mbb)),
        "mbb_68": ci(mbb, 0.16, 0.84),
        "mbb_95": ci(mbb, 0.025, 0.975),
    }
    mbb_summary[octant] = out

    print(f"\n--- {ORDERING.upper()} {octant} ---")
    print("|Ue| =", tuple(Ue_abs.tolist()))
    print("rep masses [eV] =", tuple(m_rep.tolist()), " sum=", out["sum_mnu_rep"])
    print("m_beta [eV] =", out["m_beta_rep"])
    print("m_bb [eV]: min/max =", (out["mbb_min"], out["mbb_max"]))
    print("m_bb 68%CI =", out["mbb_68"])
    print("m_bb 95%CI =", out["mbb_95"])

# LaTeX へ貼りやすい形（NO版の値）
print("\n=== LaTeX-ready: e-row |U| values ===")
print(r"\begin{align}")
print(r"(|U_{e1}|,|U_{e2}|,|U_{e3}|)")
print(r"\simeq")
print(r"\begin{cases}")
print(r"(%.6f,\ %.6f,\ %.6f) & ({\rm NO\_low\_octant})\\" % tuple(Ue_abs_targets["low_octant"]))
print(r"(%.6f,\ %.6f,\ %.6f) & ({\rm NO\_high\_octant})" % tuple(Ue_abs_targets["high_octant"]))
print(r"\end{cases}")
print(r"\end{align}")

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

# 7-3) (ADD) m_beta histogram (overlay)
plt.figure(figsize=(10,5))
for octant in pmns_targets.keys():
    sub = df[df["octant"]==octant]
    plt.hist(sub["m_beta"], bins=140, alpha=0.55, density=True, label=f"{ORDERING.upper()}_{octant}")
plt.xlabel(r"$m_\beta$ [eV]")
plt.ylabel("pdf")
plt.title(rf"{ORDERING.upper()} Gaussian MC: distribution of $m_\beta$ (low+high overlay)")
plt.legend()
mbeta_path = f"{outdir}/gaussian_mc_{ORDERING.upper()}_dist_mbeta.png"
plt.savefig(mbeta_path, dpi=200, bbox_inches="tight")
plt.show()
print("Saved:", mbeta_path)

# =========================================================
# 8) Save CSV
# =========================================================
csv_path = f"{outdir}/pmns_epsilon13_gaussian_mc_{ORDERING.upper()}_with_mbeta.csv"
df.to_csv(csv_path, index=False)
print("\nSaved:", csv_path)

# =========================================================
# 9) Save m_bb summary (CSV)
# =========================================================
mbb_rows = []
for octant, s in mbb_summary.items():
    mbb_rows.append({
        "ordering": ORDERING.upper(),
        "octant": octant,
        "Ue1": s["Ue_abs"][0], "Ue2": s["Ue_abs"][1], "Ue3": s["Ue_abs"][2],
        "m1_rep": s["m_rep"][0], "m2_rep": s["m_rep"][1], "m3_rep": s["m_rep"][2],
        "sum_mnu_rep": s["sum_mnu_rep"],
        "m_beta_rep": s["m_beta_rep"],
        "mbb_min": s["mbb_min"],
        "mbb_max": s["mbb_max"],
        "mbb_68_lo": s["mbb_68"][0], "mbb_68_hi": s["mbb_68"][1],
        "mbb_95_lo": s["mbb_95"][0], "mbb_95_hi": s["mbb_95"][1],
        "NPHASE": NPHASE,
        "seed_phase": seed_phase,
    })
df_mbb = pd.DataFrame(mbb_rows)
mbb_csv_path = f"{outdir}/mbb_phase_scan_summary_{ORDERING.upper()}.csv"
df_mbb.to_csv(mbb_csv_path, index=False)
print("Saved:", mbb_csv_path)

# =========================================================
# （数式の項の意味と役割）※主役だけまとめ（論文対応）
# =========================================================
print("\n(数式の項の意味と役割)")
print("・log2(sinθij)=-(pΔK+rΔC)：12,23 から p,r を決め、13の基準 s13_base を生成する骨格です。")
print("・epsilon13^(ℓ)=√2-log2(s13_base/s13_obs)：13だけの残差で、原理 √2 からのズレを保存する量です。")
print("・x13_eff^(ℓ)=√2-epsilon13^(ℓ)：実効指数で、論文の x13=√2-ε の形に一致させるための変数です。")
print("・(m1,m2,m3)：ORDERING(=NO/IO) と (Δm21^2,Δm3x^2) と m_geo から得る軽いニュートリノ質量です。")
print("・sum mν=m1+m2+m3：宇宙論制限と直結する位相非依存量です。")
print("・m_beta^2=Σ|U_ei|^2 m_i^2：単一β崩壊有効質量で、位相に依存しません。")
print("・m_bb=|Σ m_i |U_ei|^2 e^{iα_i}|：0νββ有効質量で、Majorana位相(α21,α31)の干渉で帯になります。")
print("・|U_ei|：pmns_fits の e行（絶対値）を入力として固定し、m_beta と m_bb の重みを決めます。")