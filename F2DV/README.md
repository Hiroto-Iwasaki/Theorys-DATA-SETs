⸻

F2DV — √2 離散真空からのフレーバー（F2DV）

F2DV — Flavor from a √2-Discrete Vacuum (F2DV)

本リポジトリは、研究論文（本稿）に付随する解析コード・生成図・補助データをまとめたものです。
This repository collects the analysis code, generated figures, and auxiliary data associated with the paper.

対象論文（作業タイトル）： “Flavor from a √2-Discrete Vacuum: Hybrid Mass Chains, GAV Mixing, and EFT-Linked Residuals”
Target paper (working title): “Flavor from a √2-Discrete Vacuum: Hybrid Mass Chains, GAV Mixing, and EFT-Linked Residuals”

含む要素：質量連鎖（主階層＋微補正）＋混合則＋\epsilon_{13} 残差保存＋dim-6 EFT 接続＋幾何学シーソー（NO/IO）
Includes: hybrid mass chains (main hierarchy + micro-correction) + mixing rule + residual \epsilon_{13} preservation + dim-6 EFT link + geometric seesaw (NO/IO)

⸻

Overview / 概要

このディレクトリには、NO/IO 分岐を含む Gaussian MC、図の生成、統計量（平均・分散・CI・相関）の算出スクリプトが含まれます。
This directory contains scripts for Gaussian Monte Carlo with NO/IO branches, figure generation, and computation of statistics (mean, variance, confidence intervals, correlations).

主な生成物は、\epsilon_{13}^{(\ell)} の分布ヒストグラム、\epsilon_{13}^{(\ell)} vs \sum m_\nu の散布図（low/high octant overlay）です。
The main outputs are histograms of \epsilon_{13}^{(\ell)} and scatter plots of \epsilon_{13}^{(\ell)} vs \sum m_\nu (low/high octant overlay).

加えて、m_\beta の分布図（low/high overlay）と、代表質量を用いた m_{\beta\beta} 位相スキャン要約（CSV）も生成します。
In addition, it generates the m_\beta distribution (low/high overlay) and a representative-mass m_{\beta\beta} phase-scan summary (CSV).

⸻

Research Context / 研究背景（今回の論文の主張）

本稿の狙いは、SM Yukawa セクターの多数パラメータを「少数の生成則＋残差保存」に分解し、観測可能量へ直接接続することです。
The goal is to decompose the many SM Yukawa parameters into “a small set of generation rules + residual preservation,” and connect them directly to observables.

具体的には、(i) 質量階層、(ii) 混合角階層、(iii) 13チャネル残差 \epsilon_{13}、(iv) dim-6 EFT、(v) 幾何学シーソーによる絶対質量固定（NO/IO）を一本化します。
Concretely, we unify (i) mass hierarchies, (ii) mixing-angle hierarchies, (iii) the 13-channel residual \epsilon_{13}, (iv) dim-6 EFT, and (v) geometric seesaw anchoring of absolute masses (NO/IO).

⸻

Core Ingredients / 中核要素（式 → 直後に「意味」と「役割」）

1) 質量連鎖：主階層＋微補正

1) Mass chain: main hierarchy + micro-correction

m_{\rm target}=m_{\rm base}(\sqrt2)^{K_{\rm int}}\exp(\beta_{\rm sector}\Delta C)

意味（各項）
Meaning (each term)
	•	m_{\rm base}, m_{\rm target}：連鎖でつなぐ基準質量・目標質量です。
m_{\rm base}, m_{\rm target}: reference and target masses linked by the chain.
	•	K_{\rm int}：\sqrt2 ステップの整数回数（主階層）です。
K_{\rm int}: an integer count of \sqrt2 steps (main hierarchy).
	•	\Delta C：Casimir差に対応する離散入力（微補正）です。
\Delta C: a discrete input corresponding to Casimir differences (micro-correction).
	•	\beta_{\rm sector}：セクター別（up/down/lep）の微補正係数です。
\beta_{\rm sector}: a sector-dependent micro-correction coefficient (up/down/lep).

役割
Role
質量階層を「\sqrt2 の整数ステップ」と「Casimir差の微補正」に分解し、少数自由度で再構成します。
This decomposes mass hierarchies into “integer \sqrt2 steps” plus “Casimir-difference micro-corrections,” enabling reconstruction with few degrees of freedom.

⸻

2) 混合則（GAV型）：CKM/PMNS を同一形式で生成

2) Mixing rule (GAV-type): generating CKM/PMNS hierarchies in a unified form

\log_2(\sin\theta_{ij})\simeq-(p\Delta K_{ij}+r\Delta C_{ij})

意味（各項）
Meaning (each term)
	•	\theta_{ij}：混合角（CKM/PMNS）です。
\theta_{ij}: a mixing angle (CKM/PMNS).
	•	\Delta K_{ij}：主階層側の離散差分ラベルです。
\Delta K_{ij}: a discrete difference label on the main-hierarchy side.
	•	\Delta C_{ij}：微補正側の離散差分ラベルです。
\Delta C_{ij}: a discrete difference label on the micro-correction side.
	•	p,r：(12),(23) 等から決まる係数です。
p,r: coefficients determined (e.g., from 12 and 23 channels).

役割
Role
混合角の階層を統一的に生成し、(12)(23) の情報から 13 のベース予測へ接続します。
This generates mixing-angle hierarchies uniformly and connects information from (12)(23) to a baseline prediction for the 13 channel.

⸻

3) 13チャネル残差の保存：\epsilon_{13}^{(\ell)}（レプトン側）

3) 13-channel residual preservation: \epsilon_{13}^{(\ell)} (lepton sector)

\epsilon_{13}^{(\ell)}=\sqrt2-\log_2\!\left(\frac{s_{13,\rm base}}{s_{13,\rm obs}}\right),\qquad
x_{13,\rm eff}^{(\ell)}=\sqrt2-\epsilon_{13}^{(\ell)}

意味（各項）
Meaning (each term)
	•	s_{13,\rm base}：(12,23) から得た (p,r) により生成される 13 の基準値です。
s_{13,\rm base}: a 13-channel baseline generated from (p,r) inferred by (12,23).
	•	s_{13,\rm obs}：観測（入力）の s_{13} です。
s_{13,\rm obs}: the observed/input s_{13}.
	•	\log_2(s_{13,\rm base}/s_{13,\rm obs})：13成分の実効指数（ズレの強さ）です。
\log_2(s_{13,\rm base}/s_{13,\rm obs}): an effective exponent quantifying the deviation.
	•	\epsilon_{13}^{(\ell)}：\sqrt2 原理値からの残差（保存されるズレ）です。
\epsilon_{13}^{(\ell)}: the residual from the \sqrt2 principle value (the preserved deviation).
	•	x_{13,\rm eff}^{(\ell)}：論文の x=\sqrt2-\epsilon に一致させる実効変数です。
x_{13,\rm eff}^{(\ell)}: an effective variable to match the paper form x=\sqrt2-\epsilon.

役割
Role
原理値 \sqrt2 を維持しつつ、現実の歪みを \epsilon_{13}^{(\ell)} として分離・保存し、EFT/UV解釈へ渡します。
While keeping \sqrt2 as the principle value, this isolates real-world distortions into \epsilon_{13}^{(\ell)} and passes them to EFT/UV interpretations.

⸻

4) dim-6 EFT 接続：13専用 Yukawa 演算子（概念枠組み）

4) dim-6 EFT link: a 13-only Yukawa operator (conceptual bridge)

\mathcal O^{(6)}_{13}=\frac{c_{13}}{\Lambda^2}(\bar Q_L\Phi d_R)(\Phi^\dagger\Phi)+{\rm h.c.}

意味（各項）
Meaning (each term)
	•	c_{13}：Wilson係数（無次元）です。
c_{13}: a Wilson coefficient (dimensionless).
	•	\Lambda：新物理スケールです。
\Lambda: the new-physics scale.
	•	\bar Q_L, d_R：左手クォークダブレット／右手ダウン型です。
\bar Q_L, d_R: left-handed quark doublet / right-handed down-type quark.
	•	\Phi：Higgs ダブレットです。
\Phi: the Higgs doublet.
	•	(\Phi^\dagger\Phi)：ゲージ不変な Higgs 二次です。
(\Phi^\dagger\Phi): a gauge-invariant Higgs bilinear.
	•	{\rm h.c.}：エルミート共役です。
{\rm h.c.}: Hermitian conjugate.

役割
Role
\epsilon_{13} を「13チャネルに局在した有効演算子の寄与」として翻訳し、\Lambda を含む検証可能な予言に接続します。
This translates \epsilon_{13} into a contribution from an effective operator localized to the 13 channel, connecting it to testable predictions involving \Lambda.

⸻

5) 幾何学シーソー入力：絶対質量スケールの固定（NO/IO）

5) Geometric seesaw input: anchoring the absolute mass scale (NO/IO)

m_\nu^{\rm geo}=m_e(\sqrt2)^{K_{\rm geo}},\qquad K_{\rm geo}=-48

意味（各項）
Meaning (each term)
	•	m_\nu^{\rm geo}：幾何学入力から与える基準ニュートリノ質量スケールです。
m_\nu^{\rm geo}: a reference neutrino mass scale given by geometric input.
	•	m_e：電子質量です。
m_e: the electron mass.
	•	K_{\rm geo}：幾何学指数（固定）です。
K_{\rm geo}: a fixed geometric index.

役割
Role
NO/IO のいずれでも絶対質量スケールを固定し、\sum m_\nu、m_\beta、m_{\beta\beta} を予言帯として出力可能にします。
This anchors the absolute mass scale for both NO and IO, enabling prediction bands for \sum m_\nu, m_\beta, and m_{\beta\beta}.

⸻

6) 観測量：\sum m_\nu, m_\beta, m_{\beta\beta}

6) Observables: \sum m_\nu, m_\beta, m_{\beta\beta}

\sum m_\nu=m_1+m_2+m_3,\qquad
m_\beta^2=\sum_i |U_{ei}|^2 m_i^2,\qquad
m_{\beta\beta}=\left|\sum_i U_{ei}^2 m_i\right|

意味（各項）
Meaning (each term)
	•	m_i：軽いニュートリノ質量固有値です。
m_i: light neutrino mass eigenvalues.
	•	\sum m_\nu：ニュートリノ質量和（宇宙論量）です。
\sum m_\nu: the neutrino mass sum (cosmological quantity).
	•	|U_{ei}|：PMNS の e 行の絶対値（重み）です。
|U_{ei}|: absolute values of the PMNS e-row (weights).
	•	m_\beta：単一β崩壊の有効質量（位相非依存）です。
m_\beta: effective mass in single beta decay (phase-independent).
	•	m_{\beta\beta}：0νββの有効質量（位相干渉で帯）です。
m_{\beta\beta}: effective mass for 0νββ (forms a band via phase interference).

役割
Role
Gaussian MC により \epsilon_{13}^{(\ell)} と同時に観測量を統計評価し、論文の「検証可能な予言帯」を直接生成します。
Gaussian MC evaluates these observables alongside \epsilon_{13}^{(\ell)}, directly producing the paper’s “testable prediction bands.”

⸻

Directory Structure / ディレクトリ構造

F2DV/
 ├── Code/
 │    ├── Gaussian_NO_and_IO.py
 │    └── delta-angle_check_for_m_bb.py
 ├── Data/
 │    ├── gaussian_mc_NO_dist_epsilon_l_13.png
 │    ├── gaussian_mc_NO_epsilon_l_13_vs_sum_mnu.png
 │    ├── gaussian_mc_IO_dist_epsilon_l_13.png
 │    ├── gaussian_mc_IO_epsilon_l_13_vs_sum_mnu.png
 │    ├── gaussian_mc_NO_dist_mbeta.png            (optional)
 │    ├── gaussian_mc_IO_dist_mbeta.png            (optional)
 │    ├── pmns_epsilon13_gaussian_mc_NO_with_mbeta.csv  (optional)
 │    ├── pmns_epsilon13_gaussian_mc_IO_with_mbeta.csv  (optional)
 │    └── mbb_phase_scan_summary_NO.csv            (optional)
 │        mbb_phase_scan_summary_IO.csv            (optional)
 └── README.md

Code/ には解析スクリプト、Data/ には論文に貼る最終図・要約CSVを置きます。
Code/ stores analysis scripts, and Data/ stores final figures and summary CSVs for the paper.

⸻

Code / コード内容と役割

Gaussian_NO_and_IO.py

NO/IO（および low/high octant）で Gaussian MC を走らせ、\epsilon_{13}^{(\ell)}、\sum m_\nu、m_\beta の統計量と図を生成します。
Runs Gaussian MC for NO/IO (and low/high octants) to generate statistics and plots for \epsilon_{13}^{(\ell)}, \sum m_\nu, and m_\beta.

代表質量（各オクタントの代表値）を用いて Majorana 位相走査を行い、m_{\beta\beta} の予言帯要約（min/max, 68%CI, 95%CI）を CSV に保存します。
Using representative masses per octant, it performs a Majorana-phase scan and saves an m_{\beta\beta} band summary (min/max, 68%CI, 95%CI) to CSV.

delta-angle_check_for_m_bb.py

同一の Majorana 位相サンプルを固定した上で δ を切り替え、m_{\beta\beta} の帯がどれだけ動くか（ほぼ普遍か）を検定します。
Fixes the same Majorana-phase samples and switches δ to test how much the m_{\beta\beta} band moves (i.e., whether it is nearly universal).

⸻

How to Use / 使用方法

1) 環境構築

1) Environment setup

pip install numpy pandas matplotlib

2) 実行（図とCSVを生成）

2) Run (generate figures and CSVs)

python Code/Gaussian_NO_and_IO.py

生成物はスクリプト内の outdir（例：GAV_outputs/）に出力されます。必要に応じて Data/ に移動してください。
Outputs are written to the script’s outdir (e.g., GAV_outputs/). Move the final files into Data/ as needed.

3) δチェック（任意）

3) δ-check (optional)

python Code/delta-angle_check_for_m_bb.py


⸻

Notes / 注意

\epsilon_{13}^{(\ell)} は「log2 比そのもの」ではなく、論文の x=\sqrt2-\epsilon 形式に一致するように 残差として定義しています。
\epsilon_{13}^{(\ell)} is defined as a residual (not simply the log2 ratio itself) to match the paper’s form x=\sqrt2-\epsilon.

（この定義により、本文の x_{13}=\sqrt2-\epsilon+\sin\theta_{13} へのマッチの整合性が保たれます。）
(This preserves consistency with the paper’s matching to x_{13}=\sqrt2-\epsilon+\sin\theta_{13}.)

⸻

Citation / 引用

本リポジトリの内容を使用する場合、以下を引用してください。
If you use materials from this repository, please cite:

Iwasaki, H. (2026).
“Flavor from a √2-Discrete Vacuum: Hybrid Mass Chains, GAV Mixing, and EFT-Linked Residuals.”
Preprint / Draft.

⸻