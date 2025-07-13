Theorys-DATA-SETs

Overview / 概要

This repository contains datasets, code, and supplementary materials related to multiple theoretical physics research papers by Hiroto Iwasaki.

本リポジトリは、Hiroto Iwasaki による複数の理論物理学研究論文に関連する データセット・コード・補足資料 を含んでいます。

⸻

Repository Structure / リポジトリ構成

Each top-level directory is named after a research project or paper, and contains:
	•	/Code/ – Scripts and computational models used in the study (e.g., Python, LaTeX)
	•	/Data/ – Datasets (e.g., CSV files, parameter tables, numerical results)
	•	README.md – Description of the paper, its objectives, data, and code

Included Papers / 含まれる研究論文（略称）
	•	DIRT
	•	Extension
	•	DDSG
	•	SCIO
	•	CIOPP
	•	IRU
	•	IPDC

Example Directory Layout

/DIRT/
  ├── Code/
  │   └── simulation_model.py
  ├── Data/
  │   └── results.csv
  └── README.md

/SCIO/
  ├── Code/
  ├── Data/
  └── README.md

Each paper directory includes its own README.md with the following structure:

# [Paper Short Title] – Supplementary Materials

## Overview / 概要
A brief explanation of the research paper, its objective, and theoretical background.

## Code / コード内容
Description of the scripts included, software dependencies, and how to run them.

## Data / データ概要
Explanation of the datasets, format (e.g., CSV), and how they relate to the paper.

## Citation / 引用情報
If you use this material, please cite the related paper:
"[Full Paper Title], Hiroto Iwasaki, [Journal/Conference], [Year]"


⸻

📦 1. Overview of Distributed Files (English)

File name	Purpose	Typical shape / dtype	Sample load
energy_bins.npy	Energy-bin edges for the low range (≈ 10–30 GeV)	(N + 1,), float64	E_lo = np.load("energy_bins.npy")
energy_bins_hi.npy	Energy-bin edges for the high range (≈ 30–100 GeV)	(M + 1,), float64	E_hi = np.load("energy_bins_hi.npy")
hi_ref.npy	Model spectrum on the high-energy side (expected counts μ)	(M,), float64	mu_hi = np.load("hi_ref.npy")
counts_obs.npy	Observed event counts (low + high concatenated)	(N + M,), int32	n_obs = np.load("counts_obs.npy")
mcmc_chain.npy	MCMC chain → n_step × n_walk × n_param (0 = τ_det, 1 = k, 2 = κ, …)	(2500, 60, 3), float64	chain = np.load("mcmc_chain.npy")
stream_vectors.npy	3-step stream vectors	(M, 3), float64	vec = np.load("stream_vectors.npy")  # each row = [Δκ, ΔΔt, Δρ]
stream_cluster_labels.npy	DBSCAN labels for the stream vectors	(M,), int32	lab = np.load("stream_cluster_labels.npy")

Notes
	•	Bin-edge arrays list the left and right edges of each energy bin (elements 0 & 1 form the first bin).
	•	hi_ref.npy stores only the model expectation, computed as
\mu_{\text{hi}} = \Phi(E_{\text{mid}})\,A_{\text{eff}}\,(4\pi)\,\Delta T.
	•	Use stream_vectors.npy together with stream_cluster_labels.npy for quiver plots or cluster-wise mean-vector analysis.
	•	lab == -1 marks noise/outliers.

⸻

📦 1. 配布ファイル概要 (Japanese)

ファイル名	役割	典型 shape / dtype	サンプルロード
energy_bins.npy	低エネルギー帯（10–30 GeV 想定）のエネルギー境界配列	(N+1,), float64	E_lo = np.load('energy_bins.npy')
energy_bins_hi.npy	高エネルギー帯（30–100 GeV）のエネルギー境界配列	(M+1,), float64	E_hi = np.load('energy_bins_hi.npy')
hi_ref.npy	高エネルギー側モデルスペクトル（期待値 μ）	(M,), float64	mu_hi = np.load('hi_ref.npy')
counts_obs.npy	実観測カウント（低 + 高を連結済み）	(N+M,), int32	n_obs = np.load('counts_obs.npy')
mcmc_chain.npy	MCMC 連鎖（shape = n_step × n_walk × n_param）0: τ_det, 1: k, 2: κ など	(2500, 60, 3), float64	chain = np.load('mcmc_chain.npy')
stream_vectors.npy	ストリームベクトル（3-step 差分）	(M, 3), float64	vec = np.load('stream_vectors.npy')  # each row = [Δκ, ΔΔt, Δρ]
stream_cluster_labels.npy	上記ベクトルの DBSCAN ラベル	(M,), int32	lab = np.load('stream_cluster_labels.npy')

メモ
	•	境界配列は「エネルギー bin の左右端」を列挙する形式です（例: 0 番目と 1 番目で 1 つの bin）。
	•	hi_ref.npy は式  \mu_{\text{hi}} = \Phi(E_{\text{mid}}) \cdot A_{\text{eff}} \cdot 4\pi\,\Delta T で計算したモデル期待値のみを格納しています。
	•	stream_vectors.npy と stream_cluster_labels.npy は、Quiver plot や「クラスタ別平均ベクトル」解析に ペアで 使用します。
	•	lab == -1 はノイズ点。

⸻

🔧 2. Derived Variables Used in the Analyses / 解析で用いる派生変数

The five core variables below match the sample script; combining x_data, y_data, dc_dt, primes, and kappa (κ) reproduces the phase analysis and DIRT parameter estimation.

Variable	Meaning (JP)	Meaning (EN)	Definition / code snippet
x_data	零点間隔 Δt 系列	Zero–spacing intervals Δt	x_data = np.load('x_data.npy')[:9999]
y_data	原系列 (情報ポテンシャル c)	Original time series (info-potential c)	y_data = np.load('y_data.npy')[:9999]
dc_dt	情報降下の時間微分 dc/dt	Time derivative of info-descent	dc_dt = np.load('dc_dt.npy')[:9999]
primes	先頭 9999 個の素数列	First 9999 prime numbers	primes = np.load('primes.npy')[:9999]
kappa	構造強度 κ(p)	Structural intensity κ(p)	kappa = 1 - np.log(np.log(primes)) / np.log(primes)  or kappa = np.load('kappa_list.npy')
prime_density	素数密度 ρ(p)	Prime density ρ(p)	prime_density = 1 / np.log(primes)
log_primes	log p 系列	Natural logs of primes	log_primes = np.log(primes)[:10000]
info_potential	情報ポテンシャル c	Information potential c	info_potential = x_data * prime_density

Usage Guidelines / 使い分けの目安

Category	Static / Dynamic	Typical use
y_data	static	FFT, derived metrics, baseline plots
dc_dt	dynamic	Detect structural breakpoints, track preservation
x_data	static	Common x-axis for periodicity analysis
primes, kappa	quasi-static	Evaluate prime-order structure, clustering features
info_potential	static	Main series for phase extraction / FFT
prime_density, log_primes	quasi-static	Visualize prime regularity, generate info_potential
stream_vectors, stream_cluster_labels	dynamic (3-step streams)	Vector-field direction (“information flow”) & cluster analysis

info_potential sometimes appears as an alias for y_data in scripts.
Prime density adopts the Riemann prime-number theorem approximation ρ(p) ≈ 1/ln p.

⸻

🚀 3. Quick-Start Mini Script / すぐに動かすミニコード例

import numpy as np
import matplotlib.pyplot as plt

# --- Load core data
E_lo  = np.load('energy_bins.npy')
E_hi  = np.load('energy_bins_hi.npy')
mu_hi = np.load('hi_ref.npy')
n_obs = np.load('counts_obs.npy')

# --- Mid-bin energies & residuals
E_mid_hi = 0.5 * (E_hi[:-1] + E_hi[1:])
residual = (n_obs[-len(mu_hi):] - mu_hi) / np.sqrt(mu_hi)

plt.loglog(E_mid_hi, residual, 'o')
plt.axhline(0, ls=':')
plt.xlabel('Energy [GeV]')
plt.ylabel('(obs − μ) / √μ')
plt.tight_layout()
plt.show()


⸻

Usage / 使い方

You can clone or download this repository to explore and reuse the research materials.

以下のコマンドでクローンできます：

git clone https://github.com/Hiroto-Iwasaki/Theorys-DATA-SETs.git

All materials are provided for academic and non-commercial use under the license stated below.

⸻

License / ライセンス

All contents in this repository are licensed under the
Creative Commons Zero v1.0 Universal (CC0).

このリポジトリ内のすべてのコンテンツは、パブリックドメインとして自由に利用可能です。出典の明示や再利用はご自由に行ってください。

⸻

Contact / 連絡先

For questions, citation details, or collaboration inquiries, please contact:
	•	Name  : Hiroto Iwasaki
	•	Email : hirotoiwasaki25@gmail.com
	•	GitHub: Hiroto-Iwasaki

⸻
