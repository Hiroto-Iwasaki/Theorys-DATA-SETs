# =========================================================
# δの角度固定が他の角度に対してもほぼ普遍であるかどうかの全検定
# delta-check for m_bb: NO_low / NO_high, delta sweep
#  - Fix the same (alpha21, alpha31) random samples for fair comparison
#  - Compare min/max and 68%CI (and 95%CI) of m_bb
# =========================================================
import numpy as np

# -----------------------------
# 0) Inputs (edit here only)
# -----------------------------
# NO-main masses (paper-fixed) [eV]
m_light_eV = np.array([0.030457910895347492, 0.03166519, 0.05854643], dtype=float)

# PMNS inputs (your current targets)
pmns_targets = {
    "NO_low_octant":  {"sin2_12": 0.303, "sin2_23": 0.451, "sin2_13": 0.02225, "delta_deg": 232.0},
    "NO_high_octant": {"sin2_12": 0.303, "sin2_23": 0.572, "sin2_13": 0.02203, "delta_deg": 197.0},
}

# Phase scan settings for m_bb
N_phase = 300000
SEED_phase = 20260111

# delta test list: use [baseline_delta, 0, 180] for each octant
delta_extra = [0.0, 180.0]

# If you also want 95%CI printed
PRINT_95CI = True

# -----------------------------
# 1) PMNS U (PDG Dirac only, no Majorana)
# -----------------------------
def U_pdg_dirac_only(s12, s23, s13, delta_rad):
    """PDG parameterization (Dirac phase only). Returns complex U_PMNS (no Majorana phases)."""
    s12 = float(s12); s23 = float(s23); s13 = float(s13)
    c12 = np.sqrt(1.0 - s12**2); c23 = np.sqrt(1.0 - s23**2); c13 = np.sqrt(1.0 - s13**2)

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
    return U

# -----------------------------
# 2) m_bb samples with Majorana phases scanned
# -----------------------------
def sample_mbb_eV(m_eV, U_e_row, alpha21, alpha31):
    """
    m_bb = | sum_i (U_ei^2 * m_i) | with Majorana phases:
      U_e2 -> U_e2 * exp(i alpha21/2)
      U_e3 -> U_e3 * exp(i alpha31/2)
    Thus: U_e2^2 gets exp(i alpha21), U_e3^2 gets exp(i alpha31)
    """
    m1, m2, m3 = map(float, m_eV)
    Ue1, Ue2, Ue3 = U_e_row

    # apply Majorana phases to squared terms
    term1 = (Ue1**2) * m1
    term2 = (Ue2**2) * np.exp(1j*alpha21) * m2
    term3 = (Ue3**2) * np.exp(1j*alpha31) * m3

    mbb = np.abs(term1 + term2 + term3)
    return mbb

def summarize_samples(x):
    x = np.asarray(x, float)
    out = {
        "min": float(np.min(x)),
        "max": float(np.max(x)),
        "q16": float(np.quantile(x, 0.16)),
        "q84": float(np.quantile(x, 0.84)),
        "q025": float(np.quantile(x, 0.025)),
        "q975": float(np.quantile(x, 0.975)),
    }
    return out

# -----------------------------
# 3) Generate ONE fixed phase sample set (alpha21, alpha31) for all deltas
# -----------------------------
rng = np.random.default_rng(SEED_phase)
alpha21 = rng.uniform(0.0, 2*np.pi, size=N_phase)
alpha31 = rng.uniform(0.0, 2*np.pi, size=N_phase)

print("=== Phase samples fixed for all delta tests ===")
print("N_phase =", N_phase, " seed_phase =", SEED_phase)
print("alpha21 range:", float(alpha21.min()), float(alpha21.max()))
print("alpha31 range:", float(alpha31.min()), float(alpha31.max()))

# -----------------------------
# 4) Run delta checks for both octants
# -----------------------------
results = []  # list of dicts

for oct_key, tgt in pmns_targets.items():
    s12 = np.sqrt(float(tgt["sin2_12"]))
    s23 = np.sqrt(float(tgt["sin2_23"]))
    s13 = np.sqrt(float(tgt["sin2_13"]))

    delta_list = [float(tgt["delta_deg"])] + [float(d) for d in delta_extra]

    for delta_deg in delta_list:
        U = U_pdg_dirac_only(s12, s23, s13, np.deg2rad(delta_deg))
        Ue = U[0, :]  # e-row (complex)

        mbb = sample_mbb_eV(m_light_eV, Ue, alpha21, alpha31)
        stats = summarize_samples(mbb)

        results.append({
            "octant": oct_key,
            "delta_deg": delta_deg,
            **stats
        })

# -----------------------------
# 5) Pretty print
# -----------------------------
def fmt(x):  # compact scientific
    return f"{x:.6e}"

print("\n=== m_bb delta-check summary (eV) ===")
for r in results:
    print(f"\n[{r['octant']}]  delta = {r['delta_deg']:.1f} deg")
    print("  min/max     =", fmt(r["min"]), "/", fmt(r["max"]))
    print("  68% CI      =", f"({fmt(r['q16'])}, {fmt(r['q84'])})")
    if PRINT_95CI:
        print("  95% CI      =", f"({fmt(r['q025'])}, {fmt(r['q975'])})")

# -----------------------------
# 6) (Optional) Quick "almost same?" diagnostic within each octant
#     Compare baseline vs delta=0,180 by relative difference of CI endpoints
# -----------------------------
def rel_diff(a, b):
    a = float(a); b = float(b)
    den = max(1e-300, abs(a))
    return abs(b - a) / den

print("\n=== Relative change vs baseline within each octant (CI endpoints) ===")
for oct_key, tgt in pmns_targets.items():
    base = next(rr for rr in results if rr["octant"]==oct_key and abs(rr["delta_deg"]-float(tgt["delta_deg"]))<1e-12)
    for d in delta_extra:
        rr = next(rr for rr in results if rr["octant"]==oct_key and abs(rr["delta_deg"]-float(d))<1e-12)
        print(f"\n[{oct_key}] delta {tgt['delta_deg']} -> {d}")
        print("  68%CI low  rel.diff =", rel_diff(base["q16"], rr["q16"]))
        print("  68%CI high rel.diff =", rel_diff(base["q84"], rr["q84"]))
        print("  min        rel.diff =", rel_diff(base["min"], rr["min"]))
        print("  max        rel.diff =", rel_diff(base["max"], rr["max"]))

# =========================================================
# （数式の項の意味と役割）※このチェックの主役
# =========================================================
print("\n(数式の項の意味と役割)")
print("・m_{ββ} = |Σ_i U_{ei}^2 m_i|：0νββ の有効質量（位相干渉で帯が決まる）。")
print("・α21, α31：Majorana 位相（m_{ββ} に本質寄与）。本セルでは [0,2π) を一様走査。")
print("・δ：Dirac 位相。ここでは δ を切替え、m_{ββ} の帯がどれだけ動くかを検査。")
print("・同一乱数列(α21,α31)固定：δ 以外の揺らぎを消して、δ 依存だけを公平に比較する役割。")
print("・min/max, 68%CI(16–84%分位)：帯の代表指標。δ を変えてもほぼ不変なら『δ固定で十分』を裏付ける。")