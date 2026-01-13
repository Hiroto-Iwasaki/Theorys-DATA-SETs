F2DV — √2 離散真空からのフレーバー（F2DV）

F2DV — Flavor from a √2-Discrete Vacuum (F2DV)

質量連鎖（主階層＋微補正）＋混合則＋\epsilon_{13} 残差保存＋dim-6 EFT 接続＋幾何学シーソー（NO/IO）
Hybrid mass chains (main hierarchy + micro-correction) + mixing rule + residual \epsilon_{13} preservation + dim-6 EFT link + geometric seesaw (NO/IO)

本リポジトリは、研究論文（本稿）に付随する解析コード・生成図・補助データをまとめたものです。
This repository collects the analysis code, generated figures, and auxiliary data associated with the paper.

対象論文（作業タイトル）：“Flavor from a √2-Discrete Vacuum: Hybrid Mass Chains, GAV Mixing, and EFT-Linked Residuals”
Target paper (working title): “Flavor from a √2-Discrete Vacuum: Hybrid Mass Chains, GAV Mixing, and EFT-Linked Residuals”

⸻

Overview / 概要

Overview / Overview

このディレクトリには、NO/IO 分岐を含む Gaussian MC、図の生成、統計量（平均・分散・CI・相関）の算出に必要なスクリプトが含まれます。
This directory contains scripts for Gaussian Monte Carlo with NO/IO branches, figure generation, and computation of statistics (mean, variance, confidence intervals, correlations).

主な生成物は、\epsilon_{13}^{(\ell)} の分布ヒストグラムと、\epsilon_{13}^{(\ell)} vs \sum m_\nu の散布図（low/high octant overlay）です。
The main outputs are histograms of \epsilon_{13}^{(\ell)} and scatter plots of \epsilon_{13}^{(\ell)} vs \sum m_\nu (low/high octant overlay).

⸻

Research Context / 研究背景（今回の論文の主張）

Research Context / Research context (claims of this paper)

本稿の狙いは、SM Yukawa セクターの多数パラメータを「少数の生成則＋残差保存」に分解し、観測可能量へ直接接続することです。
The goal of this work is to decompose the many SM Yukawa parameters into “a small set of generation rules + residual preservation,” and connect them directly to observable quantities.

具体的には、(i) 質量階層、(ii) 混合角階層、(iii) 13チャネル残差 \epsilon_{13}、(iv) dim-6 EFT、(v) 幾何学シーソーによる絶対質量固定（NO/IO）を一本化します。
Concretely, we unify (i) mass hierarchies, (ii) mixing-angle hierarchies, (iii) the 13-channel residual \epsilon_{13}, (iv) dim-6 EFT, and (v) geometric seesaw anchoring of absolute masses (NO/IO).

⸻

Core Ingredients / 中核要素（式 → 直後に「意味」と「役割」）

Core Ingredients / Core ingredients (equation → immediate “meaning” and “role”)

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

2) 混合則（GAV型）：CKM/PMNS の階層を同一形式で生成

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

3) 13チャネル残差の保存：\epsilon_{13}

3) 13-channel residual preservation: \epsilon_{13}

x_{13}=\sqrt2-\epsilon_{13}

意味（各項）
Meaning (each term)
	•	x_{13}：13抑圧の実効指数（クォーク／レプトン側で定義）です。
x_{13}: an effective exponent for 13 suppression (defined for quark/lepton sectors).
	•	\sqrt2：原理値として固定したい基準です。
\sqrt2: the reference “principle value” we aim to keep fixed.
	•	\epsilon_{13}：原理値からのズレ（残差）です。
\epsilon_{13}: deviation from the principle value (residual).

役割
Role
原理値 \sqrt2 を維持しつつ、現実の歪みを \epsilon_{13} に分離して保存し、EFT/UV解釈へ渡します。
While keeping \sqrt2 as the principle value, this isolates real-world distortions into \epsilon_{13} and passes them to EFT/UV interpretations.

⸻

4) dim-6 EFT 接続：13専用 Yukawa 演算子

4) dim-6 EFT link: a 13-only Yukawa operator

\mathcal O^{(6)}_{13}=\frac{c_{13}}{\Lambda^2}(\bar Q_L\Phi d_R)(\Phi^\dagger\Phi)+{\rm h.c.}

意味（各項）
Meaning (each term)
	•	c_{13}：Wilson係数（無次元）です。
c_{13}: a Wilson coefficient (dimensionless).
	•	\Lambda：新物理スケールです。
\Lambda: the new-physics scale.
	•	\bar Q_L, d_R：左手クォークダブレット／右手ダウン型です。
\bar Q_L, d_R: left-handed quark doublet / right-handed down-type quark.
	•	\Phi：Higgsダブレットです。
\Phi: the Higgs doublet.
	•	(\Phi^\dagger\Phi)：ゲージ不変なHiggs二次です。
(\Phi^\dagger\Phi): a gauge-invariant Higgs bilinear.
	•	h.c.：エルミート共役です。
h.c.: Hermitian conjugate.

役割
Role
\epsilon_{13} を「13チャネルに局在した有効演算子の寄与」として表現し、\Lambda を含む検証可能な予言に翻訳します。
This expresses \epsilon_{13} as a contribution from an effective operator localized to the 13 channel, translating it into testable predictions involving \Lambda.

⸻

5) 幾何学シーソー入力：絶対質量スケールの固定

5) Geometric seesaw input: anchoring the absolute mass scale

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
NO/IO のいずれでも絶対質量スケールを固定し、\sum m_\nu や m_{\beta\beta} を予言帯として出力可能にします。
This anchors the absolute mass scale for both NO and IO, enabling prediction bands for \sum m_\nu and m_{\beta\beta}.

⸻

6) 観測量：\sum m_\nu（図の縦軸）

6) Observable: \sum m_\nu (the y-axis in the plots)

\sum m_\nu=m_1+m_2+m_3

意味（各項）
Meaning (each term)
	•	m_i：質量固有値です。
m_i: mass eigenvalues.
	•	\sum m_\nu：ニュートリノ質量和（宇宙論量）です。
\sum m_\nu: the neutrino mass sum (cosmological quantity).

役割
Role
Gaussian MC により \epsilon_{13}^{(\ell)} と同時に \sum m_\nu の分布や相関を評価し、観測制限と直接比較します。
With Gaussian MC, we evaluate distributions/correlations of \sum m_\nu alongside \epsilon_{13}^{(\ell)}, enabling direct comparison with observational bounds.

⸻

What the Gaussian MC produces / Gaussian MC が生成するもの

What the Gaussian MC produces / What the Gaussian MC produces

図（NO/IO それぞれ）として、\epsilon_{13}^{(\ell)} の分布（low/high overlay）を生成します。
As figures (for both NO and IO), it generates the distribution of \epsilon_{13}^{(\ell)} with low/high octant overlay.

さらに、\epsilon_{13}^{(\ell)} と \sum m_\nu の散布図（low/high overlay）を生成します。
It also generates scatter plots of \epsilon_{13}^{(\ell)} versus \sum m_\nu with low/high overlay.

典型的に、\epsilon_{13}^{(\ell)} はオクタントで符号・中心が分離し、\sum m_\nu との相関は非常に小さい（|\mathrm{corr}|\sim 10^{-3}）ことを確認します。
Typically, \epsilon_{13}^{(\ell)} separates in sign/central value by octant, and its correlation with \sum m_\nu is extremely small (|\mathrm{corr}|\sim 10^{-3}).

⸻

Directory Structure / ディレクトリ構造

F2DV/
 ├── Code/	
 │    └── Gaussian_NO_and_IO.py
 │
 └─── Data/
      ├── gaussian_mc_NO_dist_epsilon_l_13.png
      ├── gaussian_mc_NO_epsilon_l_13_vs_sum_mnu.png
      ├── gaussian_mc_IO_dist_epsilon_l_13.png
      └── gaussian_mc_IO_epsilon_l_13_vs_sum_mnu.png

Code/ には NO/IO の Gaussian MC を実行するスクリプトを置きます。
Place the Gaussian MC script for NO/IO in Code/.

Data/ には論文に貼り付ける最終図を置きます。
Store final figures to be used in the paper in Data/.

⸻

Code / コード内容と役割

Code / Code contents and roles

Gaussian_NO_and_IO.py

このスクリプトは、NO/IO それぞれについて Gaussian MC を実行し、\epsilon_{13}^{(\ell)} と \sum m_\nu の統計と図を生成します。
This script runs Gaussian MC for NO and IO, generating statistics and plots for \epsilon_{13}^{(\ell)} and \sum m_\nu.

出力として、分布ヒスト、散布図、平均・標準偏差・68%/95% CI・相関係数を得ます。
It outputs histograms, scatter plots, and mean/std, 68%/95% CIs, and correlation coefficients.

⸻

How to Use / 使用方法

How to Use / How to use

1) 環境構築

1) Environment setup

pip install numpy matplotlib

2) 実行（NO/IO をまとめて生成）

2) Run (generate NO/IO in one go)

python Code/Gaussian_NO_and_IO.py

3) 生成物の確認

3) Check outputs

Figures/ に PNG が生成されます。
PNG figures will be generated in Figures/.

⸻

Citation / 引用

Citation / Citation

このディレクトリの内容を使用する場合、以下を引用してください。
If you use materials from this directory, please cite the following.

Iwasaki, H. (2026).
“Flavor from a √2-Discrete Vacuum: Hybrid Mass Chains, GAV Mixing, and EFT-Linked Residuals.”
Preprint / Draft.

⸻