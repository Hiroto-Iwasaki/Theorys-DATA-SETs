# 1-2
# 【目的】質量チェーン表(df_mass)の整合性を診断し、どのステップがβ一定仮定から外れているか可視化します。
# [Goal] Diagnose the mass-chain table (df_mass) and visualize which steps deviate from the constant-β assumption.

# 【機能】ステップごとの必要β（lnP/ΔC）、lnPとK差分の一致、sector一定βの残差をまとめて表示します。
# [Functions] Show per-step required β (lnP/ΔC), the lnP–Kdiff consistency check, and residuals under sector-constant β.

# 【用途】再現がズレたときに、原因がK_intかΔCかβかを切り分けるデバッグ用です。
# [Use] Debugging aid to isolate whether mismatches come from K_int, ΔC, or β when reproduction fails.

# 【注意】このセルはモデルを更新しません（診断のみ）。入力df_massとbeta_dictをそのまま検査します。
# [Note] This cell does not update the model (diagnostics only); it inspects the given df_mass and beta_dict as-is.

import numpy as np
import pandas as pd

def diagnose_mass_df(df_mass, beta_dict=None):
    df = df_mass.copy()

    # ステップごとの必要beta
    df["beta_step = lnP/ΔC"] = df["lnP(obs)"] / df["ΔC"]

    # 参考：lnP は K差そのもの（確認）
    df["lnP_from_Kdiff"] = (df["K_obs"] - df["K_int"]) * np.log(2)/2
    df["check_lnP - Kdiff"] = df["lnP(obs)"] - df["lnP_from_Kdiff"]

    # sector一定betaのときの残差 lnP - beta*ΔC
    if beta_dict is not None:
        df["beta(sector)"] = df["sector"].map(beta_dict)
        df["residual = lnP - betaΔC"] = df["lnP(obs)"] - df["beta(sector)"]*df["ΔC"]

    cols = [
        "step","sector","K_obs","K_int","ΔC",
        "lnP(obs)","beta_step = lnP/ΔC",
        "beta*ΔC","m_pred[MeV]","m_obs[MeV]","err%",
        "check_lnP - Kdiff"
    ]
    if beta_dict is not None:
        cols.insert(cols.index("beta*ΔC")+1, "residual = lnP - betaΔC")

    display(df[cols])

# すでにある beta, df_mass を使う想定
diagnose_mass_df(df_mass, beta_dict=beta)

# beta_step = lnP/ΔC
# - beta_step: 各ステップが要求する局所β / local β required by each step
# - lnP: 主階層(√2^K_int)からの対数ズレ / log deviation from the main hierarchy (√2^K_int)
# - ΔC: 微補正ラベル / micro-correction label
# - 役割: 「sector一定β」仮定が妥当かを局所的に検査 / role: test the sector-constant β assumption locally

# lnP_from_Kdiff = (K_obs - K_int) * ln2/2
# - K_obs: 観測から定義した連鎖指数 / chain exponent inferred from data
# - K_int: モデルで採用した整数ステップ / integer step used in the model
# - ln2/2: 変換係数（2*log2 ↔ ln）/ conversion factor between 2*log2 and natural log
# - 役割: lnPが「K差分由来のズレ」と一致しているか確認 / role: verify lnP matches the K-difference-based deviation