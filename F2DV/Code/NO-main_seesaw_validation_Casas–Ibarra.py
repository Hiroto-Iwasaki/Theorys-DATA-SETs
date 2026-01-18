# =========================================================
# NO主分岐（Casas–Ibarra整合）：Yν生成→mν一致・R直交性・η/LFVを検証（任意でreal Rスキャン）
# NO-main seesaw validation (Casas–Ibarra): build Yν, check mν match, R-orthogonality, and η/LFV (optional real-R scan)
# =========================================================

# 日本語フォント設定（Colab用）-----------------------------
#なぜか!aptに#を付けないと保存できなかったので、#を外してから使用してください。
#!apt-get -y install fonts-noto-cjk
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
font_path = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
font_prop = fm.FontProperties(fname=font_path)
plt.rcParams["font.family"] = font_prop.get_name()

# =========================================================
# NOメイン：Casas–Ibarra（規約整合）で Yν を生成し、
# mν一致・R直交性・η/LFV を1セルで検証 +（任意）スキャンで最適点選択
# =========================================================

# -----------------------------
# 0) Settings (触るのは基本ここだけ)
# -----------------------------
ORDERING = "NO"
OCTANT_KEY = "NO_low_octant"   # "NO_low_octant" or "NO_high_octant"

# PMNS入力（例：あなたのセルの値）
pmns_targets = {
    "NO_low_octant":  {"sin2_12": 0.303, "sin2_23": 0.451, "sin2_13": 0.02225, "delta_deg": 232.0},
    "NO_high_octant": {"sin2_12": 0.303, "sin2_23": 0.572, "sin2_13": 0.02203, "delta_deg": 197.0},
}
# Majorana位相（必要なら後でスキャン；まずは0固定）
alpha21_deg = 0.0
alpha31_deg = 0.0

# あなたの論文で固定した NO 代表値（FREEZE）
m_light_eV = np.array([0.030457910895347492, 0.03166519, 0.05854643], dtype=float)

# 参考：振動差分（FREEZE表示用）
dm21 = 7.5e-5
dm31 = 2.5e-3

# Yukawa固有値仮説（free_geom）
# ここは「IOで合わせた y0=1e-5, r=30」をそのまま使い、NO側で MR を再推定します
y0 = 1e-5
r_geom = 30.0
y_nu_diag = np.array([y0, y0*r_geom, y0*(r_geom**2)], dtype=float)

# スキャンするか（ターゲットに近い点を探すなら True）
DO_SCAN = True
Nscan = 20000
SEED_SCAN = 20260117

# ターゲット（あなたの論文記載に合わせて：例）
# ※無ければ None にしてください（スキャン最適化をしません）
target_eta_em = 5.5924e-17

# 物理定数・参照限界
v = 246.0  # GeV
alpha_em = 1/137.035999084
Br_limit = 1.5e-13  # MEG II (rough reference)

# -----------------------------
# 1) Utilities
# -----------------------------
def U_pdg_complex(s12, s23, s13, delta_rad, alpha21_rad=0.0, alpha31_rad=0.0):
    """PDG parameterization (Dirac + Majorana phases). Returns complex U_PMNS."""
    s12 = float(s12); s23 = float(s23); s13 = float(s13)
    c12 = np.sqrt(1 - s12**2); c23 = np.sqrt(1 - s23**2); c13 = np.sqrt(1 - s13**2)

    e_minus_id = np.exp(-1j*float(delta_rad))

    U = np.zeros((3,3), dtype=complex)
    U[0,0] = c12*c13
    U[0,1] = s12*c13
    U[0,2] = s13*e_minus_id

    U[1,0] = -s12*c23 - c12*s23*s13*np.conj(e_minus_id)
    U[1,1] =  c12*c23 - s12*s23*s13*np.conj(e_minus_id)
    U[1,2] =  s23*c13

    U[2,0] =  s12*s23 - c12*c23*s13*np.conj(e_minus_id)
    U[2,1] = -c12*s23 - s12*c23*s13*np.conj(e_minus_id)
    U[2,2] =  c23*c13

    # Majorana phase matrix: diag(1, e^{i alpha21/2}, e^{i alpha31/2})
    P = np.diag([1.0, np.exp(1j*alpha21_rad/2.0), np.exp(1j*alpha31_rad/2.0)])
    return U @ P

def random_orthogonal_3(rng):
    """Random real orthogonal R with det=+1 (QR)."""
    A = rng.normal(size=(3,3))
    Q, _ = np.linalg.qr(A)
    if np.linalg.det(Q) < 0:
        Q[:,0] *= -1
    return Q

def infer_MR_from_yukawa_and_mnu_diag(y_diag, mnu_eV, v_GeV=246.0, eps_floor_eV=1e-30):
    """
    Diagonal rough matching:
      m_i ≈ (v^2/2) y_i^2 / M_i  =>  M_i ≈ (v^2/2) y_i^2 / m_i
    """
    y = np.asarray(y_diag, float)
    m = np.asarray(mnu_eV, float)

    pref = (float(v_GeV)**2)/2.0  # GeV^2
    m_eff = m.copy()
    flags = []
    for i in range(3):
        if m_eff[i] <= eps_floor_eV:
            m_eff[i] = np.nan
            flags.append(f"M_R[{i}] unconstrained (m_nu[{i}] ~ 0)")
    m_GeV = m_eff * 1e-9

    MR = np.full(3, np.inf, dtype=float)
    mask = ~np.isnan(m_GeV)
    MR[mask] = pref * (y[mask]**2) / m_GeV[mask]
    return MR, flags

def casas_ibarra_Ynu_consistent(U, m_light_eV, M_heavy_GeV, R, v=246.0):
    """
    Convention (consistent set):
      mnu_target = U* Dnu U^†
      seesaw: mnu = - mD MR^{-1} mD^T
      choose: mD = i U* sqrt(Dnu) R sqrt(DR)
      Y = sqrt(2)/v * mD
    """
    m_light_GeV = np.asarray(m_light_eV, float) * 1e-9
    sqrt_m = np.diag(np.sqrt(m_light_GeV))
    sqrt_M = np.diag(np.sqrt(np.asarray(M_heavy_GeV, float)))
    mD = 1j * (U.conjugate() @ sqrt_m @ R @ sqrt_M)  # GeV
    Y  = (np.sqrt(2)/float(v)) * mD
    return Y

def mnu_target_from_U_m(U, m_light_eV):
    D = np.diag(np.asarray(m_light_eV, float) * 1e-9)  # GeV
    return U.conjugate() @ D @ U.conjugate().T  # U* D U^T  (GeV)

def mnu_seesaw_from_Y_M(Y, M_heavy_GeV, v=246.0):
    Y = np.asarray(Y, complex)
    M = np.asarray(M_heavy_GeV, float)
    mD = (float(v)/np.sqrt(2)) * Y  # GeV
    MR_inv = np.diag(1.0/M)
    return -(mD @ MR_inv @ mD.T)    # GeV

def theta_eta_from_Ynu(Y, M_heavy_GeV, v=246.0):
    mD = (float(v)/np.sqrt(2)) * np.asarray(Y, complex)  # GeV
    MR_inv = np.diag(1.0/np.asarray(M_heavy_GeV, float))
    Theta = mD @ MR_inv
    eta = 0.5 * (Theta @ Theta.conjugate().T)
    return Theta, eta

def Br_muegamma_from_eta(eta):
    return (3*alpha_em/(8*np.pi)) * (np.abs(eta[0,1])**2)

def extract_R_from_Ynu_consistent(Y, U, m_light_eV, M_heavy_GeV, v=246.0, eps_floor_GeV=0.0):
    """
    Invert the consistent convention:
      mD = (v/sqrt2) Y
      mD = i U* sqrt(m) R sqrt(M)
    => R = (-i) (sqrt(m))^{-1} (U^T mD) (sqrt(M))^{-1}
    (because U* -> U^T on left when moving across; see indices)
    """
    Y = np.asarray(Y, complex)
    U = np.asarray(U, complex)
    m_light_GeV = np.asarray(m_light_eV, float) * 1e-9
    if eps_floor_GeV > 0:
        m_light_GeV = np.maximum(m_light_GeV, eps_floor_GeV)

    sqrtm_inv = np.diag(1.0/np.sqrt(m_light_GeV))
    sqrtM_inv = np.diag(1.0/np.sqrt(np.asarray(M_heavy_GeV, float)))
    mD = (float(v)/np.sqrt(2)) * Y

    R = (-1j) * (sqrtm_inv @ (U.T @ mD) @ sqrtM_inv)
    return R

def frob_norm(A):
    return float(np.linalg.norm(A, ord="fro"))

# -----------------------------
# 2) FREEZE: NO代表値 + PMNS複素U
# -----------------------------
tgt = pmns_targets[OCTANT_KEY]
s12 = np.sqrt(float(tgt["sin2_12"]))
s23 = np.sqrt(float(tgt["sin2_23"]))
s13 = np.sqrt(float(tgt["sin2_13"]))
delta_rad = np.deg2rad(float(tgt["delta_deg"]))
alpha21_rad = np.deg2rad(float(alpha21_deg))
alpha31_rad = np.deg2rad(float(alpha31_deg))

U = U_pdg_complex(s12, s23, s13, delta_rad, alpha21_rad, alpha31_rad)
unitarity_dev = np.linalg.norm(U.conjugate().T @ U - np.eye(3))

print("=== [FREEZE] NO representative masses (paper-fixed) ===")
print("dm21 =", dm21, " dm31 =", dm31)
print("(m1,m2,m3) [eV] =", m_light_eV)
print("sum mnu [eV] =", float(np.sum(m_light_eV)))

print("\n=== [FREEZE] Complex PMNS U (PDG + Dirac phase) ===")
print("OCTANT_KEY =", OCTANT_KEY)
print("sin^2(12,23,13) =", (tgt["sin2_12"], tgt["sin2_23"], tgt["sin2_13"]))
print("delta [deg] =", tgt["delta_deg"], " alpha21,alpha31 [deg] =", (alpha21_deg, alpha31_deg))
print("Unitary check ||U^†U-I|| =", float(unitarity_dev))
print("U (complex) =\n", U)
print("|U| =\n", np.abs(U))

# -----------------------------
# 3) NO版で M_R を再推定（IO値の流用は禁止、ここで作り直す）
# -----------------------------
M_heavy_GeV, flags = infer_MR_from_yukawa_and_mnu_diag(y_nu_diag, m_light_eV, v_GeV=v, eps_floor_eV=1e-30)

print("\n=== [FREEZE] Heavy masses inferred for NO (from chosen y_nu ansatz) ===")
print("y_nu_diag =", y_nu_diag, f"(y0={y0}, r={r_geom})")
print("M_heavy_GeV =", M_heavy_GeV)
print("log10(M/GeV) =", np.log10(M_heavy_GeV))
if flags:
    print("flags:", flags)

# -----------------------------
# 4) まず1点（再現性のため R=I）で GAV候補を作って検証
# -----------------------------
R0 = np.eye(3)
Y0 = casas_ibarra_Ynu_consistent(U, m_light_eV, M_heavy_GeV, R0, v=v)

mnu_t = mnu_target_from_U_m(U, m_light_eV)         # GeV
mnu_s = mnu_seesaw_from_Y_M(Y0, M_heavy_GeV, v=v)  # GeV

rel_err = frob_norm(mnu_t - mnu_s) / max(1e-300, frob_norm(mnu_t))

Theta0, eta0 = theta_eta_from_Ynu(Y0, M_heavy_GeV, v=v)
br0 = Br_muegamma_from_eta(eta0)

R_ex0 = extract_R_from_Ynu_consistent(Y0, U, m_light_eV, M_heavy_GeV, v=v, eps_floor_GeV=0.0)
RtR0 = R_ex0.T @ R_ex0
orth_dev0 = np.linalg.norm(RtR0 - np.eye(3))

print("\n=== [TEST] One-point (R=I) consistency check ===")
print("||Y||_F =", frob_norm(Y0))
print("Relative Frobenius error ||mnu_target - mnu_seesaw||/||mnu_target|| =", float(rel_err))
print("|eta_eμ| =", float(np.abs(eta0[0,1])))
print("Br(mu->e gamma) ≈", float(br0))
print("Br/MEGII_limit =", float(br0/Br_limit))
print("Extracted R (should ~ I) =\n", R_ex0)
print("||R^T R - I|| =", float(orth_dev0))
print("||Re(R)|| =", float(np.linalg.norm(np.real(R_ex0))), " ||Im(R)|| =", float(np.linalg.norm(np.imag(R_ex0))))

# -----------------------------
# 5) （任意）real R スキャン：ターゲット |eta_eμ| に最も近い点を選ぶ
# -----------------------------
Y_best = Y0
R_best = R0
eta_best = eta0
br_best = br0
score_best = np.inf if (target_eta_em is not None) else None

if DO_SCAN:
    rng = np.random.default_rng(SEED_SCAN)
    eta_abs_e_mu = np.empty(Nscan, dtype=float)
    brs = np.empty(Nscan, dtype=float)

    best_idx = 0
    for i in range(Nscan):
        R = random_orthogonal_3(rng)
        Y = casas_ibarra_Ynu_consistent(U, m_light_eV, M_heavy_GeV, R, v=v)
        _, eta = theta_eta_from_Ynu(Y, M_heavy_GeV, v=v)
        br = Br_muegamma_from_eta(eta)

        val = float(np.abs(eta[0,1]))
        eta_abs_e_mu[i] = val
        brs[i] = float(br)

        if target_eta_em is not None:
            score = abs(val - float(target_eta_em))
            if score < score_best:
                score_best = score
                best_idx = i
                Y_best, R_best, eta_best, br_best = Y, R, eta, br

    print("\n=== [SCAN] real R scan summary ===")
    print("Nscan =", Nscan, " seed =", SEED_SCAN)
    print("Br(mu->e gamma): min/median/max =", float(brs.min()), float(np.median(brs)), float(brs.max()))
    print("frac above MEGII =", float(np.mean(brs > Br_limit)))
    print("|eta_eμ|: min/median/max =", float(eta_abs_e_mu.min()), float(np.median(eta_abs_e_mu)), float(eta_abs_e_mu.max()))

    if target_eta_em is not None:
        print("\n=== [SCAN] BEST (closest to target |eta_eμ|) ===")
        print("target |eta_eμ| =", float(target_eta_em))
        print("best   |eta_eμ| =", float(np.abs(eta_best[0,1])))
        print("score  |Δ|      =", float(score_best))
        print("Br(best)        =", float(br_best))

    # quick plots
    plt.figure()
    plt.title(r"Distribution of $\mathrm{Br}(\mu\to e\gamma)$ (real $R$ scan)")
    plt.hist(brs, bins=80, log=True)
    plt.axvline(Br_limit, linestyle="--")
    plt.xlabel("Br")
    plt.ylabel("counts (log)")
    plt.show()

    plt.figure()
    plt.title(r"Correlation: $|\eta_{e\mu}|$ vs $\mathrm{Br}(\mu\to e\gamma)$")
    plt.scatter(eta_abs_e_mu, brs, s=2, alpha=0.3)
    plt.axhline(Br_limit, linestyle="--")
    plt.xlabel(r"$|\eta_{e\mu}|$")
    plt.ylabel("Br")
    plt.yscale("log")
    plt.show()

# -----------------------------
# 6) BEST点で「mν一致」「R逆抽出直交性」を最終レポート
# -----------------------------
mnu_s_best = mnu_seesaw_from_Y_M(Y_best, M_heavy_GeV, v=v)
rel_err_best = frob_norm(mnu_t - mnu_s_best) / max(1e-300, frob_norm(mnu_t))

R_ex_best = extract_R_from_Ynu_consistent(Y_best, U, m_light_eV, M_heavy_GeV, v=v, eps_floor_GeV=0.0)
RtR_best = R_ex_best.T @ R_ex_best
orth_dev_best = np.linalg.norm(RtR_best - np.eye(3))

print("\n=== [RESULT] NO-main GAV candidate (BEST) ===")
print("v =", v)
print("M_heavy_GeV =", M_heavy_GeV)
print("||Y||_F =", frob_norm(Y_best))
print("Relative Frobenius error ||mnu_target - mnu_seesaw||/||mnu_target|| =", float(rel_err_best))
print("|eta| matrix =\n", np.abs(eta_best))
print("|eta_eμ| =", float(np.abs(eta_best[0,1])))
print("Br(mu->e gamma) ≈", float(br_best), "  Br/MEGII =", float(br_best/Br_limit))
print("\nExtracted R (from BEST Y) =\n", R_ex_best)
print("||R^T R - I|| =", float(orth_dev_best))
print("||Re(R)|| =", float(np.linalg.norm(np.real(R_ex_best))), " ||Im(R)|| =", float(np.linalg.norm(np.imag(R_ex_best))))

# --- export variables (paper-writing / next cells) ---
Y_GAV_NO = Y_best
M_HEAVY_NO = M_heavy_GeV
R_GAV_NO = R_ex_best
ETA_NO = eta_best
BR_NO = br_best

print("\n=== Exported variables ===")
print("Y_GAV_NO shape:", Y_GAV_NO.shape)
print("M_HEAVY_NO:", M_HEAVY_NO)

# =========================================================
# （数式の項の意味と役割）※このセルの主役だけまとめ
# =========================================================
print("\n(数式の項の意味と役割)")
print("・mν_target = U* Dν U^T：採用した複素PMNS U と NO代表(m_i)が意味する低エネルギーのターゲット質量行列。")
print("・mν_seesaw = - mD MR^{-1} mD^T：GAV候補(Y,M)が予言するType-I seesaw側の質量行列（モデル側）。")
print("・mD = i U* sqrt(Dν) R sqrt(DR)：上の2つを整合させるためのCasas–Ibarraの一貫規約（U* と i が重要）。")
print("・Y = (sqrt2/v) mD：Dirac質量行列をYukawa行列へ戻す定義（EWSB後の標準関係）。")
print("・R抽出（逆算）：生成したYが同じ規約のCasas–Ibarra解空間に乗っているかを、R^T R≈Iで検証する。")
print("・η = 1/2 ΘΘ^†,  Br(μ→eγ)∝|η_eμ|^2：非ユニタリ性とLFVが十分小さいか（安全領域）を評価する指標。")