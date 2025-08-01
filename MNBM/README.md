# MNBM – Supplementary Materials

## Overview / 概要
This project provides a quantitative test of the Special Quantum Relativity (SQR) framework, a novel mass-generation scenario, using public data from the IceCube DeepCore detector.The analysis focuses on 9 years of low-energy neutrino events (5–100 GeV) to empirically determine the "Information Lock Coefficient," `<R_lock>`.By applying the SQR mass formula, $m = m_0 \langle R_{lock} \rangle$, and inputting cosmological constraints on neutrino mass ($m_{target}$), we derive a value for the bare mass, $m_0 \approx 41.7$ eV.The corresponding sterile-active mixing strength, $|U_{\mu4}|^2 \approx 1.2 \times 10^{-3}$, is found to be consistent with the current exclusion limits from IceCube's 7.5-year sterile neutrino search, leaving the SQR model viable.The results are also compared with the future sensitivities of experiments like DUNE and Hyper-K.

本プロジェクトは、IceCube DeepCore検出器の公開データを用い、新しい質量生成シナリオである特殊量子相対性理論 (Special Quantum Relativity; SQR) の実証的検証を行うものです。解析では、9年間にわたる低エネルギーニュートリノ事象（5–100 GeV）に注目し、「情報ロック係数」`<R_lock>`を経験的に算出します。SQRの質量式 $m = m_0 \langle R_{lock} \rangle$ に、宇宙論的なニュートリノ質量の上限値 ($m_{target}$) を適用することで、裸の質量 $m_0 \approx 41.7$ eV を導出しました 。この結果が示唆するステライル-アクティブ混合強度 $|U_{\mu4}|^2 \approx 1.2 \times 10^{-3}$ は、IceCubeによる7.5年間のステライルニュートリノ探索の90%信頼度レベル排他領域と矛盾せず、SQRモデルが依然として許容されることを示しています。また、本解析の結果はDUNEやHyper-Kといった将来実験の感度とも比較されます 。

---

## Code / コード内容
The `/Code/` directory contains the Jupyter Notebook used for the analysis.

`/Code/` ディレクトリには、解析に使用した Jupyter Notebook が含まれています。

**`mnbm_analysis.ipynb`**
* **Description**: This notebook contains the complete Python code for the analysis. It loads public IceCube data, calculates the average Information Lock Coefficient `<R_lock>` as a function of reconstructed energy, and estimates the bare mass parameter $m_0$ based on the SQR framework. It also includes code for systematic error evaluation and visualization of the results on the sterile-active mixing plane, comparing the SQR prediction with IceCube's exclusion limits and future experimental reaches.
* **Dependencies**: `pandas`, `numpy`, `matplotlib`.
* **How to Run**: Open the notebook in a Jupyter environment (e.g., Jupyter Lab, Google Colab). Before running, download the required data files from the external sources listed in the `/Data/` section and place them in a directory accessible by the notebook. Update the file paths in the notebook accordingly.

***内容**: 解析の全工程を含むPythonコードです 。IceCubeの公開データを読み込み、SQR理論に基づいて再構成エネルギーの関数として情報ロック係数`<R_lock>`を計算し、裸質量パラメータ`m0`を推定します 。また、系統誤差の評価や、SQRの予測値をIceCubeの排他領域および将来実験の感度と比較するための、ステライル-アクティブ混合平面上での可視化コードも含まれています 。
* **依存ライブラリ**: `pandas`, `numpy`, `matplotlib`。
***実行方法**: Jupyter環境（Jupyter Lab, Google Colabなど）でノートブックを開いてください 。実行前に、`/Data/`セクションに記載されている外部リンクから必要なデータファイルをダウンロードし、ノートブックからアクセス可能なディレクトリに配置してください。ノートブック内のファイルパスは適宜修正が必要です。

---

## Data / データ概要
This repository does not contain data files directly. All data used in this analysis is publicly available from the Harvard Dataverse. Please download the necessary files from the links below. The notebook primarily uses `data.csv` (for event data) and `wilks_contour_90pct.csv` (for exclusion limits) from these datasets.

このリポジトリには、データファイルは直接含まれていません。本解析で使用されるデータはすべてHarvard Dataverseから公開されています。以下のリンクから必要なファイルをダウンロードしてください。ノートブックでは主に、これらのデータセットに含まれる`data.csv`（イベントデータ）および`wilks_contour_90pct.csv`（排他領域データ）を使用します。

* **1. Measurement of atmospheric neutrino mixing with improved IceCube DeepCore calibration and data processing**
    * **Description**: This dataset contains 9 years of low-energy muon-neutrino events from IceCube DeepCore, used for the primary `<R_lock>` calculation.
    * **Source**: [https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/B4RITM](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/B4RITM)

* **2. Results of the Search for eV-scale sterile neutrinos with 7.5 years of DeepCore data**
    * **Description**: This dataset provides the 90% C.L. exclusion contour on the sterile-active mixing plane, used for comparison with the SQR prediction.
    * **Source**: [https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/QKL28Z](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/QKL28Z)

---

## Quick-Start Script / クイックスタート
The following Python code snippet reproduces the main result of this study: a plot of the sterile-active mixing plane showing the SQR prediction against the IceCube 90% C.L. exclusion region.

以下のPythonコードは、本研究の主要な結果である、SQRの予測値とIceCubeの90%信頼度レベル排他領域を比較したステライル-アクティブ混合平面のプロットを再現するものです。

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# --- File paths (please update to your local paths)
# --- ファイルパス (ローカルのパスに更新してください)
CONTOUR_CSV = "path/to/your/data/wilks_contour_90pct.csv"

# --- Load IceCube 90% C.L. contour
# --- IceCube 90%CL等高線の読み込み
cont = pd.read_csv(CONTOUR_CSV)
x = cont['Umu4_sq']
y = cont['Utau4_sq']

# --- SQR predicted point (from the analysis)
# --- SQRの予測点 (解析結果より)
m0       = 41.7  # eV
m_target = 0.05  # eV
R_pred   = m_target / m0
u_pred   = R_pred
v_pred   = R_pred

# --- Plotting
# --- 描画
plt.figure(figsize=(5, 5))
plt.fill(x, y, alpha=0.15, label="IceCube 90% CL")
plt.plot(u_pred, v_pred, 'r*', ms=12, label=f"SQR prediction (m₀={m0:.1f} eV)")
plt.xlabel(r"$|U_{\mu4}|^{2}$")
plt.ylabel(r"$|U_{\tau4}|^{2}$")
plt.title("Sterile–Active mixing plane")
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e-4, 1)
plt.ylim(1e-4, 1)
plt.legend()
plt.grid(True, which='both', ls=':')
plt.tight_layout()
plt.show()