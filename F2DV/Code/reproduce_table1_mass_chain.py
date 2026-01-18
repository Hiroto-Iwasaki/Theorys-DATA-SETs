# 1-5
# 【目的】論文Table 1の質量チェーン予測（μ, τ, t_MS, b など）とRMSをそのまま再現します。
# [Goal] Reproduce the paper’s Table 1 mass-chain predictions (μ, τ, t_MS, b, etc.) and the RMS exactly.

# 【前提】β（lep/up/down）とmass_steps_bestは論文の固定値（式(39)–(41)およびTable 1）を使用します。
# [Assumption] Use the paper-fixed β (lep/up/down) and mass_steps_best (Eqs. (39)–(41) and Table 1).

# 【方法】最軽量（e,u,d）を観測値でアンカーし、各ステップで m_next = m_base*(√2)^K_int*exp(β*ΔC) を伝播します。
# [Method] Anchor the lightest masses (e,u,d) to data and propagate each step via m_next = m_base*(√2)^K_int*exp(β*ΔC).

# 【出力】観測値との誤差（%）と、全粒子を含む相対RMSを表示します。
# [Output] Print the percent errors vs observation and the overall relative RMS across the listed particles.

# 【注意】ここは「探索」ではなく「固定パラメータでの再現チェック」です（最適化は行いません）。
# [Note] This is not a search/fit; it is a fixed-parameter reproducibility check (no optimization).

import numpy as np
import pandas as pd

# ---------------------------------------------------------
# BEST案（あなたの出力そのもの）
# ---------------------------------------------------------
beta2 = {'lep': 0.022439401006563695, 'up': 0.06916966635852716, 'down': -0.06171501948288168}

mass_steps_best = [
    ("e->mu",   "e",  "mu",    "lep",  15, 6),
    ("mu->tau", "mu", "tau",   "lep",   8, 2),
    ("u->c",    "u",  "c",     "up",   18, 2),
    ("c->t",    "c",  TOP_KEY, "up",   14, 0),
    ("d->s",    "d",  "s",     "down",  9, 2),
    ("s->b",    "s",  "b",     "down", 11, 0),
]

# ---------------------------------------------------------
# チェーン伝播で予測質量を作る
#   m_next = m_base * (sqrt2)^K_int * exp(beta_sector * ΔC)
# ---------------------------------------------------------
def step_factor(sector, K_int, dC):
    return (np.sqrt(2.0)**int(K_int)) * np.exp(float(beta2[sector]) * float(dC))

# 予測質量辞書：最軽量（e,u,d）は観測をアンカーにして伝播
m_pred = {k: None for k in ["e","mu","tau","u","c",TOP_KEY,"d","s","b"]}
m_pred["e"] = mass_obs["e"]
m_pred["u"] = mass_obs["u"]
m_pred["d"] = mass_obs["d"]

# 伝播（順に適用）
for (name, base_k, tgt_k, sector, K_int, dC) in mass_steps_best:
    if m_pred[base_k] is None:
        raise RuntimeError(f"Base mass not set for {base_k} in chain.")
    m_pred[tgt_k] = m_pred[base_k] * step_factor(sector, K_int, dC)

# 表：観測と比較
rows = []
for k in ["e","mu","tau","u","c",TOP_KEY,"d","s","b"]:
    obs = mass_obs[k]
    pred = m_pred[k]
    errpct = 100.0*(pred-obs)/obs
    rows.append([k, pred, obs, errpct])

df_chain = pd.DataFrame(rows, columns=["particle","m_pred_chain[MeV]","m_obs[MeV]","err%"])
display(df_chain)

rms = float(np.sqrt(np.mean((df_chain["err%"]/100.0)**2)))
print("RMS(relative, chain) =", rms)

# m_next = m_base * (sqrt2)^K_int * exp(beta_sector * ΔC)
# - m_next: 次の粒子の予測質量 / predicted mass of the next particle
# - m_base: 基準（前段）粒子の質量 / base (previous) particle mass
# - sqrt2^K_int: 主階層の離散ステップ / discrete main-hierarchy step
# - beta_sector: セクター別の微補正係数 / sector-dependent micro-correction coefficient
# - ΔC: Casimir差ラベル（離散入力）/ Casimir-difference label (discrete input)
# - 役割: 「√2離散階層＋微補正」で質量連鎖を再構成 / role: reconstruct the mass chain via √2 discreteness + micro-corrections